# Sorghum defense architecture

This repository provides a reproducible Python pipeline for the multi-scale defense architecture analyses in sorghum:
1) canopy geometry and greenhouse severity (FSP35),
2) tissue-level susceptibility contrasts (leaf vs midrib; FSP35 excised-leaf assays),
3) information-theoretic alignment (Jensen–Shannon divergence) between normalized attack patterns and measured defense allocation (FSP53 qPCR).

This repository focuses on derived tables and Excel supplementary outputs and does not generate manuscript figures.

## Files in this repository

Primary inputs (raw data):
- `8leaf_final.xlsx`
- `susceptibility score_Midrib.xlsx`
- `qpcr_leaf_midrib_allocation.csv`

Scripts:
- `run_pipeline.py` (end-to-end pipeline; raw inputs → derived tables → `outputs/Supplementary_Data_S1.xlsx`)
- `compute_funneling_r2.py` (quick reproducibility check for canopy funneling R² values)

## Environment

Recommended conda environment:

```bash
conda env create -f environment.yml
conda activate gwas_env



## Quick check: canopy funneling R² (n=14 complete-case)

From the repository root:

python compute_funneling_r2.py --data-dir .

Run the full pipeline

python run_pipeline.py --data-dir . --out-dir outputs --overwrite

Expected output includes:
Outputs:

outputs/Supplementary_Data_S1.xlsx (primary deliverable)

outputs/derived_*.csv (intermediate tables for auditing)

Notes on probabilities vs distributions

p_leaf and p_midrib are absolute probabilities of severe infection under the FSP35 excised-leaf assay and do not sum to 1.

For JSD, they are normalized to:
P_attack = [p_leaf/(p_leaf+p_midrib), p_midrib/(p_leaf+p_midrib)].

JSD is computed between normalized P_attack and P_defense using log base 2.
