# GitHub dataset push batch plan

Each batch must be staged, committed, and pushed before moving to the next batch.
Batch constraints: <= 1.25 GiB payload and <= 10000 payload paths.
The helper stages only; it never commits or pushes.

## Batch 1/3
- Payload: 1.24 GiB
- Files/parts: 53

```bash
python datasets/github/stage_git_dataset_batch.py 1
git status --short
git commit -m "Add dataset batch 1/3"
git push origin HEAD
```

## Batch 2/3
- Payload: 1.24 GiB
- Files/parts: 53

```bash
python datasets/github/stage_git_dataset_batch.py 2
git status --short
git commit -m "Add dataset batch 2/3"
git push origin HEAD
```

## Batch 3/3
- Payload: 570.44 MiB
- Files/parts: 7975

```bash
python datasets/github/stage_git_dataset_batch.py 3
git status --short
git commit -m "Add dataset batch 3/3"
git push origin HEAD
```
