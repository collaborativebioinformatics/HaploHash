# What we're trying to do with eQTLs

**Goal:** check whether our AVI-based haploblock ranking (how "functionally important/disruptive" a block is, from AlphaGenome) actually means something biologically — using a completely independent dataset as a sanity check, instead of just trusting the model.

**The check:** eQTLs are variants that are *known, from real population data (GTEx)*, to change how much a nearby gene gets expressed. If our AVI ranking is capturing something real about selection/constraint, we expected: the more constrained a block is, the fewer eQTLs it should have — because purifying selection should be actively removing disruptive regulatory variants before they can spread through the population.

**What we did first (block level):** ranked all 39,074 haploblocks by AVI score, counted GTEx eQTLs (pooled across all 49 tissues) landing in each block, checked the correlation.

**Result:** opposite of what we expected. More AVI-flagged blocks have *more* eQTLs, not fewer (moderate positive correlation, very statistically solid given 39k blocks). Same pattern also showed up independently in the cluster/re-identification work — more AVI-constrained blocks also have more haplotype clusters and more singleton haplotypes.

**Why, probably:** AVI score at the block level mostly reflects how much regulatory/gene machinery is packed into that stretch of DNA (a "how much biology is happening here" measure), not how depleted of actual variation the region is. Gene-dense regions naturally have both more positions that *could* be damaging (pushing AVI up) and more genes to be eQTL targets for (pushing eQTL count up) — same underlying cause, not a contradiction of selection theory. It just means block-level AVI density isn't the right lens for testing selection directly.

**The better test (in progress):** instead of "how much potential damage exists in this block," ask "of all the positions AVI flags as potentially damaging, how many actually show up as real variants in people" — a proper obs/exp ratio, same idea as gnomAD's constraint scores. Low ratio = selection is actively suppressing those sites. We piloted this on the 12 chr22 test blocks we already had data for, and the direction flipped to match the original hypothesis (fewer eQTLs where obs/exp is lower) — but 12 blocks is way too small a sample to call it real (not statistically significant yet).

**What's needed to finish it:** the same obs/exp calculation done genome-wide instead of just chr22. That requires phased 1000G VCFs for every chromosome — right now we only have chr22's downloaded. If anyone has already pulled the other chromosomes' VCFs, or knows where they live in our shared data folder, that would unblock this immediately. Otherwise it's a straightforward download + rerun once the cluster's reachable again.
