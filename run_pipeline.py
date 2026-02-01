#!/usr/bin/env python3
"""
End-to-end pipeline to reproduce derived tables and Supplementary Data S1.

Inputs (expected in --data-dir):
  - 8leaf_final.xlsx
  - susceptibility score_Midrib.xlsx
  - qpcr_leaf_midrib_allocation.csv

Outputs (written to --out-dir):
  - Supplementary_Data_S1.xlsx
  - derived_canopy_features.csv
  - derived_midrib_probs.csv
  - derived_qpcr_allocation.csv
  - derived_alignment_jsd.csv

This pipeline does not generate manuscript figures.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Dict, Any

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

import openpyxl
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter


# ---------------------------
# Utilities
# ---------------------------

def _require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input file: {path}")


def _safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def jsd_base2(p: np.ndarray, q: np.ndarray) -> float:
    """
    Jensen–Shannon divergence (base 2) for discrete distributions p and q.
    p and q must be 1D arrays that sum to 1.
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)

    if p.ndim != 1 or q.ndim != 1 or p.shape != q.shape:
        raise ValueError("p and q must be 1D arrays of the same shape.")

    # Normalize defensively
    p_sum = p.sum()
    q_sum = q.sum()
    if p_sum <= 0 or q_sum <= 0:
        return float("nan")
    p = p / p_sum
    q = q / q_sum

    m = 0.5 * (p + q)

    def kl(a: np.ndarray, b: np.ndarray) -> float:
        # 0 * log(0/.) treated as 0
        mask = a > 0
        return float(np.sum(a[mask] * (np.log2(a[mask]) - np.log2(b[mask]))))

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


# ---------------------------
# Canopy features (8leaf_final.xlsx)
# ---------------------------

@dataclass
class CanopyModelResults:
    r2_lambda_linear: float
    r2_funneling_full: float
    r2_funneling_slope_only: float
    n_complete_case: int


def build_canopy_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build genotype-level canopy features from leaf angle measurements.

    Required columns:
      - cultivar
      - leaf_no
      - angle_deg
      - GH_FSP35 (may contain NaN)

    Returns a genotype-level table with:
      mean_profile_angle, profile_sd, min_angle, max_angle, angle_range,
      leaf_slope, intercept,
      mean_basal_angle (mean of all angles for the genotype),
      runoff_index, retention_index,
      GH_FSP35_mean (mean across rows for genotype, ignoring NaN)
    """
    required = {"cultivar", "leaf_no", "angle_deg", "GH_FSP35"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"8leaf_final.xlsx missing required columns: {sorted(missing)}")

    rows = []
    for cultivar, g in df.groupby("cultivar", dropna=True):
        g = g.copy()
        g["leaf_no"] = pd.to_numeric(g["leaf_no"], errors="coerce")
        g["angle_deg"] = pd.to_numeric(g["angle_deg"], errors="coerce")
        g["GH_FSP35"] = pd.to_numeric(g["GH_FSP35"], errors="coerce")

        g = g.dropna(subset=["leaf_no", "angle_deg"])
        if g.empty:
            continue

        # per-leaf mean angle profile (leaf 1..8)
        leaf_means = g.groupby("leaf_no")["angle_deg"].mean().sort_index()

        mean_profile_angle = float(leaf_means.mean())
        profile_sd = float(leaf_means.std(ddof=1))
        min_angle = float(leaf_means.min())
        max_angle = float(leaf_means.max())
        angle_range = float(max_angle - min_angle)

        # slope of profile means vs leaf_no
        if leaf_means.shape[0] >= 2:
            x = leaf_means.index.to_numpy(dtype=float)
            y = leaf_means.to_numpy(dtype=float)
            slope, intercept = np.polyfit(x, y, 1)
        else:
            slope, intercept = np.nan, np.nan

        mean_basal_angle = float(g["angle_deg"].mean())

        # physics indices based on mean basal angle
        theta = np.deg2rad(mean_basal_angle)
        runoff_index = float(np.cos(theta))
        retention_index = float(np.sin(theta))

        gh_vals = g["GH_FSP35"].dropna()
        gh_mean = float(gh_vals.mean()) if len(gh_vals) > 0 else np.nan

        rows.append(
            {
                "cultivar": str(cultivar),
                "mean_profile_angle": mean_profile_angle,
                "profile_sd": profile_sd,
                "min_angle": min_angle,
                "max_angle": max_angle,
                "angle_range": angle_range,
                "leaf_slope": float(slope),
                "intercept": float(intercept),
                "mean_basal_angle": mean_basal_angle,
                "runoff_index": runoff_index,
                "retention_index": retention_index,
                "GH_FSP35_mean": gh_mean,
            }
        )

    return pd.DataFrame(rows)


def fit_canopy_models(feats: pd.DataFrame, eps: float = 1e-6) -> Tuple[pd.DataFrame, CanopyModelResults]:
    """
    Fit:
      - linear model: lambda_arch ~ mean_basal_angle + leaf_slope + profile_sd
      - funneling full: GH_FSP35_mean ~ mean_profile_angle + leaf_slope + profile_sd + leaf_slope^2
      - funneling slope-only: GH_FSP35_mean ~ leaf_slope + leaf_slope^2

    Note: lambda_arch uses p_severe_GH normalized from GH_FSP35_mean to [0,1] using max scaling.
    """
    df = feats.copy()

    # complete-case for GH models
    cc = df.dropna(subset=["GH_FSP35_mean", "mean_profile_angle", "leaf_slope", "profile_sd", "mean_basal_angle"]).reset_index(drop=True)
    n = cc.shape[0]
    if n < 3:
        raise ValueError(f"Not enough complete-case genotypes for canopy models (n={n}).")

    # Normalize GH severity to [0,1] via max scaling (matches earlier manuscript usage in this project)
    gh = cc["GH_FSP35_mean"].to_numpy(dtype=float)
    gh_max = float(np.nanmax(gh))
    if gh_max <= 0:
        raise ValueError("GH_FSP35_mean max is non-positive; cannot normalize.")
    p_severe = gh / gh_max
    lambda_arch = -np.log(p_severe + float(eps))
    cc["p_severe_GH"] = p_severe
    cc["lambda_arch"] = lambda_arch

    # linear lambda model
    X_l = cc[["mean_basal_angle", "leaf_slope", "profile_sd"]].to_numpy(dtype=float)
    y_l = cc["lambda_arch"].to_numpy(dtype=float)
    m_l = LinearRegression().fit(X_l, y_l)
    r2_lambda = r2_score(y_l, m_l.predict(X_l))

    # funneling full (predict GH directly)
    slope = cc["leaf_slope"].to_numpy(dtype=float)
    X_f = np.column_stack([
        cc["mean_profile_angle"].to_numpy(dtype=float),
        slope,
        cc["profile_sd"].to_numpy(dtype=float),
        slope ** 2,
    ])
    y_f = cc["GH_FSP35_mean"].to_numpy(dtype=float)
    m_f = LinearRegression().fit(X_f, y_f)
    r2_full = r2_score(y_f, m_f.predict(X_f))

    # slope-only quadratic
    X_s = np.column_stack([slope, slope ** 2])
    m_s = LinearRegression().fit(X_s, y_f)
    r2_slope_only = r2_score(y_f, m_s.predict(X_s))

    results = CanopyModelResults(
        r2_lambda_linear=float(r2_lambda),
        r2_funneling_full=float(r2_full),
        r2_funneling_slope_only=float(r2_slope_only),
        n_complete_case=int(n),
    )

    # Add predictions back for auditing
    cc["lambda_arch_pred"] = m_l.predict(X_l)
    cc["GH_pred_full"] = m_f.predict(X_f)
    cc["GH_pred_slope_only"] = m_s.predict(X_s)

    return cc, results


# ---------------------------
# Tissue susceptibility (susceptibility score_Midrib.xlsx)
# ---------------------------

def compute_midrib_probs(df: pd.DataFrame, severe_threshold: float = 3.0) -> pd.DataFrame:
    """
    Compute p_leaf and p_midrib as absolute probabilities of severe infection for FSP35 in excised-leaf assays.

    Severe is defined as FSP35_score >= severe_threshold.

    Important: Rows with missing FSP35_score are retained and treated as non-severe in the denominator.
    This matches the convention used in the manuscript tables where unscored sites contribute to the total
    number of evaluated sites for a cultivar/tissue.
    """
    required = {"cultivar", "FSP35_score", "tissue"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"susceptibility score_Midrib.xlsx missing required columns: {sorted(missing)}")

    d = df.copy()
    d["FSP35_score"] = pd.to_numeric(d["FSP35_score"], errors="coerce")
    d = d.dropna(subset=["cultivar", "tissue"])

    # Comparison with NaN yields False, which is the desired behavior here.
    d["is_severe"] = d["FSP35_score"] >= float(severe_threshold)

    out_rows = []
    for cultivar, g in d.groupby("cultivar", dropna=True):
        for tissue in ["Leaf", "Midrib"]:
            gg = g[g["tissue"].astype(str).str.strip().str.lower() == tissue.lower()]
            if gg.empty:
                continue
            p = float(gg["is_severe"].mean())  # denominator includes NaN scores (treated as non-severe)
            out_rows.append({"cultivar": str(cultivar), "tissue": tissue, "p_severe": p, "n_sites": int(gg.shape[0])})

    out = pd.DataFrame(out_rows)

    # Pivot to wide
    wide = out.pivot_table(index="cultivar", columns="tissue", values="p_severe", aggfunc="first")
    wide = wide.rename(columns={"Leaf": "p_leaf", "Midrib": "p_midrib"}).reset_index()
    wide["delta_p"] = wide["p_leaf"] - wide["p_midrib"]
    return wide


# ---------------------------
# qPCR allocation (qpcr_leaf_midrib_allocation.csv)
# ---------------------------

def compute_defense_allocation(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert fold-changes to a two-element defense allocation distribution per cultivar:
      p_def_leaf = leaf_fc / (leaf_fc + mid_fc)
      p_def_midrib = 1 - p_def_leaf

    Then average across genes and time points for each cultivar.
    """
    required = {"cultivar", "leaf_fc", "mid_fc"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"qpcr_leaf_midrib_allocation.csv missing required columns: {sorted(missing)}")

    d = df.copy()
    d["leaf_fc"] = pd.to_numeric(d["leaf_fc"], errors="coerce")
    d["mid_fc"] = pd.to_numeric(d["mid_fc"], errors="coerce")
    d = d.dropna(subset=["cultivar", "leaf_fc", "mid_fc"])

    denom = d["leaf_fc"] + d["mid_fc"]
    d = d[denom > 0].copy()
    d["p_def_leaf"] = d["leaf_fc"] / denom
    d["p_def_midrib"] = 1.0 - d["p_def_leaf"]

    # Average across gene/time per cultivar
    summary = d.groupby("cultivar", dropna=True)[["p_def_leaf", "p_def_midrib"]].mean().reset_index()
    summary["cultivar"] = summary["cultivar"].astype(str)
    return summary


# ---------------------------
# Alignment table (Table 2) and portfolio
# ---------------------------

def build_alignment_table(midrib: pd.DataFrame, defense: pd.DataFrame) -> pd.DataFrame:
    """
    Join midrib probabilities (p_leaf, p_midrib) with defense allocation (p_def_leaf, p_def_midrib),
    normalize attack, and compute JSD (base 2).
    """
    m = midrib.copy()
    d = defense.copy()

    merged = m.merge(d, on="cultivar", how="inner")
    merged = merged.dropna(subset=["p_leaf", "p_midrib", "p_def_leaf", "p_def_midrib"]).reset_index(drop=True)

    p_sum = merged["p_leaf"] + merged["p_midrib"]
    merged["P_attack_leaf_norm"] = merged["p_leaf"] / p_sum
    merged["P_attack_midrib_norm"] = merged["p_midrib"] / p_sum

    jsd_vals = []
    for _, r in merged.iterrows():
        p_attack = np.array([r["P_attack_leaf_norm"], r["P_attack_midrib_norm"]], dtype=float)
        p_def = np.array([r["p_def_leaf"], r["p_def_midrib"]], dtype=float)
        jsd_vals.append(jsd_base2(p_attack, p_def))
    merged["JSD_attack_defense"] = jsd_vals
    merged["risk_JSD"] = merged["JSD_attack_defense"]
    merged["return_defense"] = 1.0 - (merged["p_leaf"] + merged["p_midrib"]) / 2.0

    # Standardize column naming for manuscript-friendly exports
    out = merged[[
        "cultivar",
        "p_leaf", "p_midrib",
        "P_attack_leaf_norm", "P_attack_midrib_norm",
        "p_def_leaf", "p_def_midrib",
        "JSD_attack_defense", "risk_JSD", "return_defense",
        "delta_p",
    ]].copy()

    return out


# ---------------------------
# Excel writer
# ---------------------------

def write_supplementary_excel(
    out_path: Path,
    canopy_raw: pd.DataFrame,
    canopy_features: pd.DataFrame,
    canopy_cc: pd.DataFrame,
    midrib_probs: pd.DataFrame,
    qpcr_summary: pd.DataFrame,
    alignment: pd.DataFrame,
    model_results: CanopyModelResults,
) -> None:
    """
    Write a submission-friendly Supplementary_Data_S1.xlsx with a README sheet and key derived tables.
    """
    wb = openpyxl.Workbook()
    # Remove default sheet
    wb.remove(wb.active)

    bold = Font(bold=True)

    def add_sheet_from_df(name: str, df: pd.DataFrame) -> None:
        ws = wb.create_sheet(name)
        # Header
        for j, col in enumerate(df.columns, start=1):
            cell = ws.cell(row=1, column=j, value=str(col))
            cell.font = bold
            ws.column_dimensions[get_column_letter(j)].width = max(14, min(32, len(str(col)) + 2))
        # Rows
        for i, (_, r) in enumerate(df.iterrows(), start=2):
            for j, col in enumerate(df.columns, start=1):
                val = r[col]
                ws.cell(row=i, column=j, value=_safe_float(val) if isinstance(val, (np.floating, float, int, np.integer)) else val)
        ws.freeze_panes = "A2"

    # README
    ws0 = wb.create_sheet("00_README", 0)
    ws0["A1"] = "Supplementary Data S1"
    ws0["A1"].font = Font(bold=True, size=14)
    ws0["A3"] = "Raw inputs"
    ws0["A3"].font = bold
    ws0["B3"] = "8leaf_final.xlsx; susceptibility score_Midrib.xlsx; qpcr_leaf_midrib_allocation.csv"
    ws0["A5"] = "Canopy model (complete-case)"
    ws0["A5"].font = bold
    ws0["B5"] = f"n = {model_results.n_complete_case}; R2_full_quadratic = {model_results.r2_funneling_full:.3f}; R2_slope_only_quadratic = {model_results.r2_funneling_slope_only:.3f}; R2_lambda_linear = {model_results.r2_lambda_linear:.3f}"
    ws0["A7"] = "Probability vs distribution (important)"
    ws0["A7"].font = bold
    ws0["B7"] = ("p_leaf and p_midrib are absolute probabilities of severe infection (FSP35 excised-leaf assays) "
                 "and do not sum to 1. For JSD, they are normalized to P_attack = [p_leaf/(p_leaf+p_midrib), p_midrib/(p_leaf+p_midrib)].")
    ws0["A9"] = "JSD"
    ws0["A9"].font = bold
    ws0["B9"] = "JSD is computed between normalized P_attack and P_defense using log base 2."
    ws0.column_dimensions["A"].width = 32
    ws0.column_dimensions["B"].width = 120
    ws0["B3"].alignment = Alignment(wrap_text=True)
    ws0["B7"].alignment = Alignment(wrap_text=True)

    # Sheets
    add_sheet_from_df("Angle_by_leaf", canopy_raw[["cultivar", "leaf_no", "angle_deg", "GH_FSP35"]].copy())
    add_sheet_from_df("Arch_summary", canopy_features.copy())
    add_sheet_from_df("Physics_LambdaArch_cc", canopy_cc.copy())
    add_sheet_from_df("Midrib_probs", midrib_probs.copy())
    add_sheet_from_df("qPCR_defense_alloc", qpcr_summary.copy())
    add_sheet_from_df("Table2_alignment_clean", alignment.copy())

    wb.save(out_path)


# ---------------------------
# Main
# ---------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=str, required=True, help="Directory containing the three raw input files.")
    ap.add_argument("--out-dir", type=str, required=True, help="Output directory.")
    ap.add_argument("--eps", type=float, default=1e-6, help="Small constant epsilon for lambda_arch.")
    ap.add_argument("--severe-threshold", type=float, default=3.0, help="Threshold for severe infection in excised-leaf scores (FSP35_score).")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs if present.")
    args = ap.parse_args()

    data_dir = Path(args.data_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    arch_path = data_dir / "8leaf_final.xlsx"
    midrib_path = data_dir / "susceptibility score_Midrib.xlsx"
    qpcr_path = data_dir / "qpcr_leaf_midrib_allocation.csv"

    _require_file(arch_path)
    _require_file(midrib_path)
    _require_file(qpcr_path)

    # Load inputs
    canopy_raw = pd.read_excel(arch_path)
    midrib_raw = pd.read_excel(midrib_path)
    qpcr_raw = pd.read_csv(qpcr_path)

    # Derived tables
    canopy_feats = build_canopy_features(canopy_raw)
    canopy_cc, model_results = fit_canopy_models(canopy_feats, eps=float(args.eps))

    midrib_probs = compute_midrib_probs(midrib_raw, severe_threshold=float(args.severe_threshold))
    qpcr_summary = compute_defense_allocation(qpcr_raw)

    alignment = build_alignment_table(midrib_probs, qpcr_summary)

    # Write intermediate CSVs
    csv_map = {
        "derived_canopy_features.csv": canopy_feats,
        "derived_canopy_complete_case.csv": canopy_cc,
        "derived_midrib_probs.csv": midrib_probs,
        "derived_qpcr_allocation.csv": qpcr_summary,
        "derived_alignment_jsd.csv": alignment,
    }
    for name, df in csv_map.items():
        path = out_dir / name
        if path.exists() and not args.overwrite:
            continue
        df.to_csv(path, index=False)

    # Write Excel
    xlsx_path = out_dir / "Supplementary_Data_S1.xlsx"
    if xlsx_path.exists() and not args.overwrite:
        raise FileExistsError(f"{xlsx_path} already exists. Use --overwrite to replace it.")
    write_supplementary_excel(
        out_path=xlsx_path,
        canopy_raw=canopy_raw,
        canopy_features=canopy_feats,
        canopy_cc=canopy_cc,
        midrib_probs=midrib_probs,
        qpcr_summary=qpcr_summary,
        alignment=alignment,
        model_results=model_results,
    )

    print("Wrote outputs to:", out_dir)
    print("Primary deliverable:", xlsx_path)
    print(f"Canopy models (complete-case n={model_results.n_complete_case}): "
          f"R2_full={model_results.r2_funneling_full:.3f}, "
          f"R2_slope_only={model_results.r2_funneling_slope_only:.3f}, "
          f"R2_lambda_linear={model_results.r2_lambda_linear:.3f}")


if __name__ == "__main__":
    main()
