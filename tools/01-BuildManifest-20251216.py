#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
from pathlib import Path


def _hash_file(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Copy files listed in tt into repo, and generate MANIFEST + checksums."
    )
    parser.add_argument("--from-tt", default="tt", help="Path to tt file (default: tt).")
    parser.add_argument(
        "--write",
        action="store_true",
        help="Actually copy files and write MANIFEST/checksums (default: dry-run).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite destination files if already exist.",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    tt_path = repo_root / args.from_tt
    lines = [l.strip() for l in tt_path.read_text(encoding="utf-8").splitlines() if l.strip()]

    records = []
    for i, src_rel in enumerate(lines, start=1):
        if not src_rel.startswith("../"):
            raise SystemExit(f"Unsupported path in tt (expected ../...): {src_rel}")
        dst_rel = src_rel[3:]
        src = (repo_root / src_rel).resolve()
        dst = repo_root / dst_rel
        if not src.exists():
            raise SystemExit(f"Missing source file: {src_rel} -> {src}")

        if args.write:
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists() and not args.overwrite:
                pass
            else:
                shutil.copy2(src, dst)

        size = dst.stat().st_size if dst.exists() else src.stat().st_size
        md5 = _hash_file(dst if dst.exists() else src, "md5")
        sha256 = _hash_file(dst if dst.exists() else src, "sha256")
        records.append((i, src_rel, dst_rel, size, md5, sha256))

    if args.write:
        manifest_path = repo_root / "MANIFEST.tsv"
        with manifest_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["index", "source_rel", "repo_rel", "bytes", "md5", "sha256"])
            w.writerows(records)

        (repo_root / "tt.repo").write_text(
            "\n".join([r[2] for r in records]) + "\n", encoding="utf-8"
        )

        (repo_root / "CHECKSUMS.md5").write_text(
            "\n".join([f"{r[4]}  {r[2]}" for r in records]) + "\n", encoding="utf-8"
        )
        (repo_root / "CHECKSUMS.sha256").write_text(
            "\n".join([f"{r[5]}  {r[2]}" for r in records]) + "\n", encoding="utf-8"
        )

        print(f"Wrote {len(records)} records to MANIFEST.tsv / checksums.")
    else:
        print("Dry-run OK. Use --write to copy and write MANIFEST/checksums.")
        print(f"Would process {len(records)} files.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
