#!/usr/bin/env bash
set -euo pipefail

python tools/02-VerifyRepo-20251216.py
Rscript tests/package-test-function-20251216/01-ParseAllR-20251216.R

echo "Preflight OK."
