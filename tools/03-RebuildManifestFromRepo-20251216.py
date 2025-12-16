#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def _hash_file(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _posix_relpath(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rebuild MANIFEST.tsv + checksums by scanning current repo for .R files."
    )
    parser.add_argument("--root", default=".", help="Repo root to scan (default: .).")
    parser.add_argument(
        "--write",
        action="store_true",
        help="Actually write MANIFEST.tsv/CHECKSUMS.* (default: dry-run).",
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    if not root.exists():
        raise SystemExit(f"Root not found: {root}")

    r_files = []
    for p in root.rglob("*.R"):
        if not p.is_file():
            continue
        if "/.git/" in p.as_posix():
            continue
        r_files.append(p)

    r_files = sorted(set(r_files))
    if not r_files:
        print("No .R files found.")
        return 1

    records: list[tuple[int, str, str, int, str, str]] = []
    for i, path in enumerate(r_files, start=1):
        rel = _posix_relpath(path, root)
        size = path.stat().st_size
        md5 = _hash_file(path, "md5")
        sha256 = _hash_file(path, "sha256")
        records.append((i, rel, rel, size, md5, sha256))

    if not args.write:
        print("Dry-run OK. Use --write to write MANIFEST.tsv/checksums.")
        print(f"Would write {len(records)} records.")
        return 0

    manifest_path = root / "MANIFEST.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["index", "source_rel", "repo_rel", "bytes", "md5", "sha256"])
        w.writerows(records)

    (root / "CHECKSUMS.md5").write_text(
        "\n".join([f"{r[4]}  {r[2]}" for r in records]) + "\n", encoding="utf-8"
    )
    (root / "CHECKSUMS.sha256").write_text(
        "\n".join([f"{r[5]}  {r[2]}" for r in records]) + "\n", encoding="utf-8"
    )

    print(f"Wrote {len(records)} records to MANIFEST.tsv / checksums.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

