#!/usr/bin/env Rscript
# haploblock_pathway_characterization_plots.R
#
# Combined characterization + plotting script for haploblock functional and
# pathway analysis pipeline: annotate haploblock clusters to genes, then 
# characterize the biological content of those genes via pathway/GO enrichment 
# and standard enrichplot visualizations. 
#
# Outputs (all under --out-dir)
# ------------------------------
#   per_haploblock_gene_annotation.csv   descriptive, no statistics
#   panel_pathway_enrichment_go_bp.csv   \
#   panel_pathway_enrichment_go_cc.csv    > whole-genome background (no
#   panel_pathway_enrichment_go_mf.csv   /  `universe` argument passed)
#   dotplot_go_bp.pdf / _cc.pdf / _mf.pdf   GO dotplots, run separately
#   manhattanplot_go_all.pdf                GO BP+CC+MF landscape together
#   panel_pathway_enrichment_kegg.csv        + cnetplot_kegg.pdf, emapplot_kegg.pdf       (--run-kegg)
#   panel_pathway_enrichment_reactome.csv    + cnetplot_reactome.pdf, emapplot_reactome.pdf (--run-reactome)
#
# Install
# -------
#   install.packages(c("optparse", "dplyr", "ggplot2"))
#   BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot"))
#
# Usage
# -----
#   Rscript haploblock_pathway_characterization_plots.R \
#     --hap2gene out_dir/functional_annotation/haploblock_to_gene.tsv \
#     --run-kegg --run-reactome \
#     --out-dir out_dir/pathway_characterization
#
#   # restrict to specific haploblocks of interest:
#   Rscript haploblock_pathway_characterization_plots.R \
#     --hap2gene out_dir/functional_annotation/haploblock_to_gene.tsv \
#     --haploblock-ids chr6_31480875-31598421,chr6_32000000-32100000 \
#     --out-dir out_dir/pathway_characterization

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ReactomePA)
})

option_list <- list(
  make_option("--hap2gene", type = "character", dest = "hap2gene",
              help = "haploblock_to_gene.tsv from the Python align+intersect step"),
  make_option("--haploblock-ids", type = "character", default = NULL, dest = "haploblock_ids",
              help = "comma-separated haploblock_ids to restrict to (default: entire panel)"),
  make_option("--run-kegg", action = "store_true", default = FALSE, dest = "run_kegg"),
  make_option("--run-reactome", action = "store_true", default = FALSE, dest = "run_reactome"),
  make_option("--show-category", type = "integer", default = 15, dest = "show_category"),
  make_option("--out-dir", type = "character", dest = "out_dir")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(name) file.path(opt$out_dir, name)
has_results <- function(res) !is.null(res) && nrow(as.data.frame(res)) > 0
save_plot <- function(p, name, width = 9, height = 7) {
  ggsave(out(name), p, width = width, height = height, device = "pdf")
  message(sprintf("Wrote %s", out(name)))
}

# --------------------------------------------------------------------------- #
# Load and (optionally) restrict the haploblock -> gene mapping
# --------------------------------------------------------------------------- #

hap2gene <- read.delim(opt$hap2gene, stringsAsFactors = FALSE)

if (!is.null(opt$haploblock_ids)) {
  ids <- strsplit(opt$haploblock_ids, ",")[[1]]
  before <- n_distinct(hap2gene$haploblock_id)
  hap2gene <- hap2gene %>% filter(haploblock_id %in% ids)
  message(sprintf(
    "Restricted to %d/%d requested haploblocks (of %d in panel)",
    n_distinct(hap2gene$haploblock_id), length(ids), before
  ))
}

gene_list <- unique(na.omit(hap2gene$gene_name))
message(sprintf("Panel covers %d genes across %d haploblocks", length(gene_list), n_distinct(hap2gene$haploblock_id)))

# --------------------------------------------------------------------------- #
# Output 1: per-haploblock descriptive annotation, no statistics
# --------------------------------------------------------------------------- #

annotation <- hap2gene %>%
  group_by(haploblock_id) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    genes = paste(sort(unique(na.omit(gene_name))), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genes))

write.csv(annotation, out("per_haploblock_gene_annotation.csv"), row.names = FALSE)

n_zero <- sum(annotation$n_genes == 0)
if (n_zero > 0) message(sprintf("NOTE: %d/%d haploblocks overlap zero genes", n_zero, nrow(annotation)))

# --------------------------------------------------------------------------- #
# Output 2: GO enrichment, run ONCE per ontology -- each result feeds both
# its CSV and its dotplot, and one extra ont="ALL" run feeds the manhattan
# plot. Whole-genome background throughout (no `universe` argument).
# --------------------------------------------------------------------------- #

go_ontologies <- c("BP", "CC", "MF")
for (ont in go_ontologies) {
  ego <- enrichGO(
    gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
    ont = ont, pAdjustMethod = "BH"
  )
  csv_name <- sprintf("panel_pathway_enrichment_go_%s.csv", tolower(ont))
  write.csv(as.data.frame(ego), out(csv_name), row.names = FALSE)

  if (has_results(ego)) {
    p <- dotplot(ego, showCategory = opt$show_category, label_format = NULL) +
      ggtitle(sprintf("GO %s dotplot", ont))
    save_plot(p, sprintf("dotplot_go_%s.pdf", tolower(ont)))
  } else {
    message(sprintf("No significant GO %s terms -- skipping dotplot", ont))
  }
}

ego_all <- enrichGO(
  gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = "ALL", pAdjustMethod = "BH"
)
if (has_results(ego_all)) {
  p_manhattan <- manhattanplot(
    ego_all,
    color = "p.adjust",
    showCategory = 12,
    size = "Count",
    split = "ONTOLOGY",
    title = "GO enrichment landscape (BP, CC, MF)"
  )
  save_plot(p_manhattan, "manhattanplot_go_all.pdf", width = 11, height = 6)
} else {
  message("No significant GO terms across any ontology -- skipping manhattan plot")
}

# --------------------------------------------------------------------------- #
# Output 3: KEGG and Reactome, each computed ONCE -- CSV + cnetplot + emapplot
# from the same object.
# --------------------------------------------------------------------------- #

plot_cnet_and_emap <- function(res, label) {
  if (!has_results(res)) {
    message(sprintf("No significant %s terms -- skipping cnetplot/emapplot", label))
    return(invisible(NULL))
  }
  res_readable <- if (isTRUE(res@readable)) res else setReadable(res, "org.Hs.eg.db", "ENTREZID")

  p_cnet <- cnetplot(res_readable, showCategory = min(opt$show_category, 10)) +
    ggtitle(sprintf("%s gene-concept network", label))
  save_plot(p_cnet, sprintf("cnetplot_%s.pdf", tolower(label)), width = 10, height = 9)

  res_sim <- tryCatch(pairwise_termsim(res_readable), error = function(e) NULL)
  if (!is.null(res_sim)) {
    p_emap <- emapplot(res_sim, showCategory = opt$show_category) +
      ggtitle(sprintf("%s enrichment map", label))
    save_plot(p_emap, sprintf("emapplot_%s.pdf", tolower(label)), width = 10, height = 9)
  } else {
    message(sprintf("Could not compute pairwise term similarity for %s -- skipping emapplot", label))
  }
}

if (opt$run_kegg || opt$run_reactome) {
  map_gene <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  dropped <- length(gene_list) - nrow(map_gene)
  if (dropped > 0) {
    message(sprintf(
      "NOTE: %d/%d genes did not map to an Entrez ID and were dropped from KEGG/Reactome",
      dropped, length(gene_list)
    ))
  }

  if (opt$run_kegg) {
    kegg_res <- enrichKEGG(gene = map_gene$ENTREZID, organism = "hsa", pAdjustMethod = "BH")
    write.csv(as.data.frame(kegg_res), out("panel_pathway_enrichment_kegg.csv"), row.names = FALSE)
    plot_cnet_and_emap(kegg_res, "KEGG")
  }

  if (opt$run_reactome) {
    reactome_res <- enrichPathway(
      gene = map_gene$ENTREZID, organism = "human",
      pAdjustMethod = "BH", readable = TRUE
    )
    write.csv(as.data.frame(reactome_res), out("panel_pathway_enrichment_reactome.csv"), row.names = FALSE)
    plot_cnet_and_emap(reactome_res, "Reactome")
  }
}

