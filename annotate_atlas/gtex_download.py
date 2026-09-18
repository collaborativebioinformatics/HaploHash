#!/usr/bin/env python3
"""Download original GTEx v8 results from GTEx and the eQTL Catalogue mirror."""

import argparse
import csv
import subprocess
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "gtex/data"
COMMON = {
    "GTEx_Analysis_v8_eQTL.tar": "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar",
    "GTEx_v8_finemapping_DAPG.tar": "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/fine-mapping-cis-eqtl/GTEx_v8_finemapping_DAPG.tar",
    "README_eQTL_v8.txt": "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/README_eQTL_v8.txt",
}


def download(url, path, connections=1):
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".part")
    print(f"Downloading {path.name}", flush=True)
    if connections > 1 or Path(str(partial) + ".aria2").exists():
        subprocess.run(["aria2c", "--continue=true", "--auto-file-renaming=false",
                        "--file-allocation=none", f"--split={connections}",
                        f"--max-connection-per-server={connections}", "--max-tries=5",
                        "--retry-wait=5", "--connect-timeout=60", "--timeout=120",
                        "--summary-interval=60", "--console-log-level=warn", "--show-console-readout=false",
                        "--dir", str(partial.parent), "--out", partial.name, url], check=True)
    else:
        subprocess.run(["curl", "--silent", "--show-error", "--fail", "--location", "--retry", "5", "--retry-delay", "5",
                        "--connect-timeout", "60", "--speed-time", "120", "--speed-limit", "1000",
                        "--continue-at", "-", "--output", str(partial), url], check=True)
    partial.rename(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", help="tissue row number, or 'common'")
    parser.add_argument("--connections", type=int, choices=range(1, 17), default=1,
                        help="parallel connections using aria2c")
    args = parser.parse_args()
    if args.index == "common":
        for name, url in COMMON.items():
            download(url, DATA / name, args.connections)
        with tarfile.open(DATA / "GTEx_Analysis_v8_eQTL.tar") as archive:
            for member in archive:
                if member.isfile() and member.name.endswith(".egenes.txt.gz"):
                    target = DATA / "egenes" / Path(member.name).name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.extractfile(member) as src, target.open("wb") as dst:
                        import shutil
                        shutil.copyfileobj(src, dst)
        print("Common GTEx files ready", flush=True)
    else:
        with (HERE / "gtex/tissues.tsv").open() as stream:
            tissue = list(csv.DictReader(stream, delimiter="\t"))[int(args.index)]
        download(tissue["url"], DATA / "allpairs" / (tissue["tissue"] + ".tsv.gz"), args.connections)
        print(tissue["tissue"], "downloaded", flush=True)


if __name__ == "__main__":
    main()
