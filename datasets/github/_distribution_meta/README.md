# GitHub dataset distribution

This directory is a GitHub-safe representation of the project's authoritative local `data/` dataset.

## Representation

- `data/` is the authoritative **local dataset** and is excluded from Git.
- `datasets/github/` is the **GitHub distribution dataset**.
- Ordinary files <=24 MiB are copied byte-for-byte.
- Larger ordinary files are split into fixed 24 MiB binary parts.
- Source directories wider than 3000 immediate entries are represented by deterministic `tar.gz` directory bundles; bundles are also split at 24 MiB when necessary.
- Bundle reconstruction restores the original directory members and verifies every member SHA-256.
- The machine-readable manifest is stored as deterministic `dataset_manifest.json.gz` to keep metadata comfortably below GitHub file-size limits.
- Nested VCS metadata (`.git`, `.hg`, `.svn`) is intentionally local-only and listed in `LOCAL_ONLY_EXCLUSIONS.csv`.
- `datasets/reconstructed/` is ignored by Git.

Distributable source entries: **100784**  
Directory bundles: **10**  
Full local `data/` footprint: **4.26 GiB**

## Reconstruct / verify

Reconstruct safely into the default local-only `datasets/reconstructed/`:

```bash
python datasets/github/_distribution_meta/recombine_dataset.py
```

Verify the authoritative local dataset without modifying it:

```bash
python datasets/github/_distribution_meta/recombine_dataset.py --output data --verify-only
```

Non-verify reconstruction into authoritative `data/` is explicitly refused.

## Initial GitHub upload

The generated payload requires **3** staged Git batches. Each batch is capped by both payload size and path count. Stage, inspect, commit, and push one batch before proceeding to the next. The helper never commits or pushes automatically.

See `_distribution_meta/PUSH_BATCH_PLAN.md`.

## Integrity digests

Full local tree digest: `940711f6cb205fb0f02dd2a38f8242e8b927e4c81b5f7a027398e553c576bd88`  
Distributable scientific tree digest: `0ef2b3e2770fd65f1d7ffdfb2bb6bf8dafd71bcbfef1b038be8245d2814655da`
