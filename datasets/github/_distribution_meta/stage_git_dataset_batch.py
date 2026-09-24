#!/usr/bin/env python3
"""Stage exactly one generated GitHub dataset batch; does not commit or push."""
from __future__ import annotations
import argparse, json, subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
GITHUB_ROOT = HERE.parent
ROOT = GITHUB_ROOT.parent.parent
BATCHES = json.loads((HERE / "GIT_UPLOAD_BATCHES.json").read_text(encoding="utf-8"))
META = [
    ".gitignore",
    "datasets/github/_distribution_meta/README.md",
    "datasets/github/_distribution_meta/dataset_manifest.json.gz",
    "datasets/github/_distribution_meta/GITHUB_DATASET_INVENTORY.csv",
    "datasets/github/_distribution_meta/LARGE_FILE_REPORT.csv",
    "datasets/github/_distribution_meta/BUNDLED_DIRECTORY_REPORT.csv",
    "datasets/github/_distribution_meta/LOCAL_ONLY_EXCLUSIONS.csv",
    "datasets/github/_distribution_meta/GIT_UPLOAD_BATCHES.json",
    "datasets/github/_distribution_meta/PUSH_BATCH_PLAN.md",
    "datasets/github/_distribution_meta/recombine_dataset.py",
    "datasets/github/_distribution_meta/stage_git_dataset_batch.py",
]

def add_paths(paths):
    for i in range(0, len(paths), 300):
        subprocess.run(["git", "-C", str(ROOT), "add", "--", *paths[i:i+300]], check=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("batch", type=int)
    args = ap.parse_args()
    if args.batch < 1 or args.batch > len(BATCHES["batches"]):
        raise SystemExit(f"batch must be 1..{len(BATCHES['batches'])}")
    batch = BATCHES["batches"][args.batch - 1]
    paths = ["datasets/github/" + p for p in batch["github_relative_paths"]]
    if args.batch == 1: add_paths(META)
    add_paths(paths)
    print(f"STAGED batch {args.batch}/{len(BATCHES['batches'])}: {len(paths)} payload paths")
    print("Review with: git -C", ROOT, "status --short")
    print("This helper did NOT commit or push anything.")

if __name__ == "__main__": main()
