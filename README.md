# Data availability

The raw transcriptomic data for this study were retrieved from the EMBL-EBI ArrayExpress database under accession number E-MTAB-8142. To focus on the immune landscape, subsequent downstream analyses were specifically restricted to the CD45⁺ cell population. Detailed metadata regarding the samples utilized in this research are available in the [sample_information
](https://github.com/boiscat/ILC_single_cell/blob/main/sample_information)

# Analysis Scripts Pack (R)

This repository is a packaged snapshot of a set of R analysis scripts (Seurat / Harmony / CellChat / enrichment), organized for GitHub sharing and reproducible execution.

## What's included

- R scripts copied from the original working directory using `tt` as the source list
- A traceable manifest (`MANIFEST.tsv`) and checksums (`CHECKSUMS.md5`, `CHECKSUMS.sha256`)
- Step-by-step Chinese docs under `docs/`
- Environment helper files under `env/`

## Quick start

1) Create an environment (recommended: mamba/conda):

```bash
mamba env create -f env/environment.yml
mamba activate pso-analysis
```

2) Install R packages (CRAN + Bioconductor + GitHub):

```bash
Rscript env/install.R
```

3) Basic syntax check (does not run analysis):

```bash
Rscript tests/package-test-function-20251216/01-ParseAllR-20251216.R
```

## How to run scripts

- Many scripts use **relative paths**, so the **working directory matters**.
- For scripts under `figures/fig1/` that read paths like `fig1/...`, run from `figures/`:

```bash
cd figures
Rscript fig1/01-patientILC-20251204.R
```

- For pipeline scripts under `scripts/`, run from the repository root:

```bash
Rscript scripts/prep/01-sample-site-extract-20251108.R
```

## Reproducibility / traceability

- `tt` is the original source list (paths outside the repo).
- `tt.repo` is the in-repo file list (paths inside the repo).
- `MANIFEST.tsv` records source path → repo path and file hashes.

See `docs/REPRODUCIBILITY.md` for details.

## Data

Input data files are **not included**. See `docs/DATA.md` for where to put data and expected file names.

## Push to GitHub

See `docs/GITHUB.md` (step-by-step).

---
