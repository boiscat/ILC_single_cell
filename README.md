# Study Overview
Innate lymphoid cells (ILCs) are critical regulators of immune homeostasis in barrier tissues, yet the mechanisms underlying this balance are not fully elucidated. This study identifies a subset of CD73⁺ ILCs that plays a protective role in skin and colonic inflammation.

# Data availability

The raw transcriptomic data for this study were retrieved from the EMBL-EBI ArrayExpress database under accession number E-MTAB-8142. To focus on the immune landscape, subsequent downstream analyses were specifically restricted to the CD45⁺ cell population. Detailed metadata regarding the samples utilized in this research are available in the [sample_information
](https://github.com/boiscat/ILC_single_cell/blob/main/sample_information)

# Analysis Scripts Pack (R)

This repository is a packaged snapshot of a set of R analysis scripts (Seurat / Harmony / CellChat / enrichment), organized for GitHub sharing and reproducible execution.

## What's included

- R scripts copied from the original working directory
- A traceable manifest (`MANIFEST.tsv`) 
- Environment helper files under `env/`

## Environment

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


```

- For pipeline scripts under `scripts/`, run from the repository root:

```bash
Rscript scripts/prep/01-sample-site-extract-20251108.R
```
