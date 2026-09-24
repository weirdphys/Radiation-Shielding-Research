#!/usr/bin/env python3
"""Reconstruct and verify the distributable scientific data tree."""
from __future__ import annotations
import argparse, gzip, hashlib, json, os, shutil, sys, tarfile, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
GITHUB_ROOT = HERE.parent
ROOT = GITHUB_ROOT.parent.parent
MANIFEST = HERE / "dataset_manifest.json.gz"


def sha256_file(path: Path, chunk: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def safe_tar(tf: tarfile.TarFile):
    for m in tf.getmembers():
        p = Path(m.name)
        if p.is_absolute() or ".." in p.parts or not (m.isfile() or m.issym()):
            raise RuntimeError(f"unsafe/unexpected archive member: {m.name}")


def materialize_packaged(obj, temp_path: Path):
    mode = obj["packaging_mode"]
    if mode in ("direct", "bundle_direct"):
        shutil.copy2(GITHUB_ROOT / obj["github_relative_path"], temp_path)
    elif mode in ("split", "bundle_split"):
        with temp_path.open("wb") as out:
            for part in sorted(obj["parts"], key=lambda x: x["part_index"]):
                pp = GITHUB_ROOT / part["github_relative_path"]
                if sha256_file(pp) != part["sha256"]:
                    raise RuntimeError(f"part hash mismatch: {pp}")
                with pp.open("rb") as src:
                    shutil.copyfileobj(src, out, length=8 * 1024 * 1024)
    else:
        raise RuntimeError(f"unsupported packaged-object mode: {mode}")


def verify_existing_member(out_root: Path, e):
    dst = out_root / e["relative_path"]
    if e["entry_type"] == "symlink":
        return dst.is_symlink() and os.readlink(dst) == e["symlink_target"]
    return dst.is_file() and sha256_file(dst) == e["sha256"]


def install_member(src: Path, dst: Path, e, overwrite: bool):
    if dst.exists() or dst.is_symlink():
        if e["entry_type"] == "symlink" and dst.is_symlink() and os.readlink(dst) == e["symlink_target"]:
            return False
        if e["entry_type"] == "file" and dst.is_file() and sha256_file(dst) == e["sha256"]:
            return False
        if not overwrite:
            raise RuntimeError(f"exists with different content/type: {e['relative_path']}")
        if dst.is_dir() and not dst.is_symlink():
            shutil.rmtree(dst)
        else:
            dst.unlink(missing_ok=True)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if e["entry_type"] == "symlink":
        dst.symlink_to(e["symlink_target"])
    else:
        tmp = dst.with_name(dst.name + ".reconstructing")
        tmp.unlink(missing_ok=True)
        shutil.copy2(src, tmp)
        if sha256_file(tmp) != e["sha256"]:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(f"member SHA mismatch before install: {e['relative_path']}")
        os.replace(tmp, dst)
        try:
            os.chmod(dst, int(e.get("mode", 0o644)))
        except Exception:
            pass
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=ROOT / "datasets" / "reconstructed")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--verify-only", action="store_true")
    args = ap.parse_args()

    with gzip.open(MANIFEST, "rt", encoding="utf-8") as f:
        manifest = json.load(f)
    out_root = args.output.expanduser().resolve()
    authoritative = (ROOT / "data").resolve()
    if out_root == authoritative and not args.verify_only:
        raise SystemExit("REFUSED: non-verify reconstruction into authoritative data/ is not allowed")

    failures = []
    written = 0

    # Unbundled entries.
    for e in manifest["entries"]:
        dst = out_root / e["relative_path"]
        if e["entry_type"] == "symlink":
            ok = dst.is_symlink() and os.readlink(dst) == e["symlink_target"]
            if args.verify_only:
                if not ok: failures.append(f"symlink mismatch: {e['relative_path']}")
                continue
            if ok: continue
            if dst.exists() or dst.is_symlink():
                if not args.overwrite:
                    failures.append(f"exists with different content/type: {e['relative_path']}")
                    continue
                if dst.is_dir() and not dst.is_symlink(): shutil.rmtree(dst)
                else: dst.unlink(missing_ok=True)
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.symlink_to(e["symlink_target"])
            written += 1
            continue

        expected = e["sha256"]
        if dst.is_file() and sha256_file(dst) == expected:
            continue
        if args.verify_only:
            failures.append(f"missing/hash mismatch: {e['relative_path']}")
            continue
        if (dst.exists() or dst.is_symlink()) and not args.overwrite:
            failures.append(f"exists with different hash/type: {e['relative_path']}")
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="rr_file_reconstruct_") as td:
            tmp = Path(td) / "payload"
            materialize_packaged(e, tmp)
            if sha256_file(tmp) != expected:
                raise RuntimeError(f"reconstructed SHA-256 mismatch: {e['relative_path']}")
            if dst.exists() or dst.is_symlink():
                if dst.is_dir() and not dst.is_symlink(): shutil.rmtree(dst)
                else: dst.unlink(missing_ok=True)
            shutil.copy2(tmp, dst)
            try: os.chmod(dst, int(e.get("mode", 0o644)))
            except Exception: pass
            written += 1

    # Bundled directories.
    for b in manifest.get("directory_bundles", []):
        members = b["members"]
        if args.verify_only:
            for e in members:
                if not verify_existing_member(out_root, e):
                    failures.append(f"bundle member mismatch: {e['relative_path']}")
            continue

        with tempfile.TemporaryDirectory(prefix="rr_bundle_reconstruct_") as td:
            td = Path(td)
            archive = td / "bundle.tar.gz"
            materialize_packaged(b, archive)
            if sha256_file(archive) != b["archive_sha256"]:
                raise RuntimeError(f"bundle archive hash mismatch: {b['bundle_directory']}")
            extracted = td / "extracted"
            extracted.mkdir()
            with tarfile.open(archive, mode="r:gz") as tf:
                safe_tar(tf)
                tf.extractall(extracted, filter="fully_trusted")
            bd = Path(b["bundle_directory"])
            for e in members:
                inner = Path(e["relative_path"]).relative_to(bd)
                src = extracted / inner
                if e["entry_type"] == "symlink":
                    if not src.is_symlink() or os.readlink(src) != e["symlink_target"]:
                        raise RuntimeError(f"extracted bundle symlink mismatch: {e['relative_path']}")
                else:
                    if not src.is_file() or sha256_file(src) != e["sha256"]:
                        raise RuntimeError(f"extracted bundle file mismatch: {e['relative_path']}")
                dst = out_root / e["relative_path"]
                try:
                    if install_member(src, dst, e, args.overwrite): written += 1
                except RuntimeError as exc:
                    failures.append(str(exc))

    if failures:
        print("RECONSTRUCTION/VERIFICATION FAILED", file=sys.stderr)
        for x in failures[:200]: print(" -", x, file=sys.stderr)
        if len(failures) > 200: print(f" ... {len(failures)-200} additional failures", file=sys.stderr)
        raise SystemExit(2)

    total = manifest["distribution_entry_count"]
    print(f"PASS: {total} distributable entries verified/restorable; {written} entries written.")
    print("Output data root:", out_root)
    print("Distribution tree digest:", manifest["distribution_tree_digest_sha256"])
    if manifest.get("local_only_excluded_entry_count", 0):
        print("NOTE: local-only VCS metadata is intentionally absent; see _distribution_meta/LOCAL_ONLY_EXCLUSIONS.csv")

if __name__ == "__main__":
    main()
