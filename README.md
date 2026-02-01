# Sorghum defense architecture

This repository provides a reproducible Python pipeline for the multi-scale defense architecture analyses in sorghum:
1) canopy geometry and greenhouse severity (FSP35),
2) tissue-level susceptibility contrasts (leaf vs midrib; FSP35 excised-leaf assays),
3) information-theoretic alignment (Jensen–Shannon divergence) between normalized attack patterns and measured defense allocation (FSP53 qPCR).

The pipeline focuses on **derived tables and Excel supplementary outputs** and does **not** generate manuscript figures.

## Repository layout

- `scripts/run_pipeline.py`: end-to-end pipeline (raw inputs → derived tables → `Supplementary_Data_S1.xlsx`)
- `scripts/compute_funneling_r2.py`: quick reproducibility check for the canopy funneling R² values
- `data/raw/`: place the three primary input files here (not committed by default)
- `outputs/`: generated outputs (ignored by git)

## Inputs

Place the following three files into `data/raw/`:

- `8leaf_final.xlsx`
- `susceptibility score_Midrib.xlsx`
- `qpcr_leaf_midrib_allocation.csv`

If you cannot or do not want to commit raw data to GitHub, keep `data/raw/` untracked and provide the files locally.

## Environment

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate gwas_env
```

## Quick check: canopy funneling R² (n=14 complete-case)

From the repository root:

```bash
python scripts/compute_funneling_r2.py --data-dir data/raw
```

Expected output includes:

- `R2_slope_quad` ≈ 0.838  
- `R2_full_quad`  ≈ 0.870  
- `n` = 14 complete-case cultivars

## Run the full pipeline

```bash
python scripts/run_pipeline.py --data-dir data/raw --out-dir outputs
```

Outputs:

- `outputs/Supplementary_Data_S1.xlsx` (primary deliverable)
- `outputs/derived_*.csv` (intermediate tables for auditing)

## Notes on probabilities vs distributions

- `p_leaf` and `p_midrib` are **absolute** probabilities of severe infection under the FSP35 excised-leaf assay and do **not** sum to 1.
- For Jensen–Shannon divergence (JSD), they are normalized to a two-element distribution:

`P_attack = [p_leaf/(p_leaf+p_midrib), p_midrib/(p_leaf+p_midrib)]`.


MIT License (see `LICENSE`).
