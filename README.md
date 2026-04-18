# Sorghum defense architecture

This repository contains the input data and analysis pipeline used to generate the supplementary data workbook for the manuscript:

**Canopy leaf-angle architecture and leaf blade–midrib susceptibility are associated with anthracnose severity in sorghum and johnsongrass**

## Contents

### Input data

- `8leaf_final.xlsx`
- `susceptibility score_Midrib.xlsx`
- `qpcr_leaf_midrib_allocation.csv`

### Analysis script

- `run_pipeline.py`

The pipeline reads the input files and writes one Excel workbook containing the derived tables, model summaries, audit sheets, and manuscript-ready table outputs.

## Environment

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate sorghum-defense-architecture
```

## Run

From the repository root:

```bash
python run_pipeline.py --data-dir . --out-dir outputs --overwrite
```

The output workbook will be written to:

```text
outputs/Supplementary_Data_S1.xlsx
```

## Missing-value handling for tissue probabilities

Severe-infection probabilities for the FSP35 detached-leaf assay are calculated from scored observations only. Rows with missing `FSP35_score` values are excluded from the denominator for the corresponding genotype and tissue. The output workbook includes an audit sheet reporting total rows, non-missing score rows, and missing score rows by genotype and tissue.

## Output workbook

The workbook includes:

- input audit
- leaf-angle summaries
- canopy trait summaries
- canopy PCA scores and explained variance
- canopy model data and model summary
- tissue severe-infection probability audit
- tissue severe-infection probabilities
- genotype-level ordinal tissue comparisons
- pooled cumulative-logit model tests
- qRT-PCR defense-allocation summaries
- corrected attack-defense alignment table
- Table 2
- Table S1

The pipeline does not generate manuscript figures.
