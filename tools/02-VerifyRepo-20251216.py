#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


def _hash_file(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    repo_root = Path.cwd()
    manifest = repo_root / "MANIFEST.tsv"
    if not manifest.exists():
        raise SystemExit("MANIFEST.tsv not found. Run tools/01-BuildManifest-20251216.py --write first.")

    bad = 0
    with manifest.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rel = row["repo_rel"]
            path = repo_root / rel
            if not path.exists():
                print(f"[MISSING] {rel}")
                bad += 1
                continue
            md5 = _hash_file(path, "md5")
            sha256 = _hash_file(path, "sha256")
            if md5 != row["md5"] or sha256 != row["sha256"]:
                print(f"[HASH MISMATCH] {rel}")
                bad += 1

    if bad:
        print(f"FAILED: {bad} file(s) missing or mismatched.")
        return 2
    print("OK: all files match MANIFEST.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

