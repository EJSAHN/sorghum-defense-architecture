#!/usr/bin/env python3
"""
Compute funneling-model R^2 values from 8leaf_final.xlsx.

This script reports:
  - R2_slope_quad: GH_FSP35 ~ slope + slope^2
  - R2_full_quad : GH_FSP35 ~ mean_angle_profile + slope + sd_profile + slope^2

Both are computed on complete-case cultivars (cultivars with non-missing GH_FSP35).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


REQUIRED_COLUMNS = ["cultivar", "leaf_no", "angle_deg", "GH_FSP35"]


def build_architecture_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Per cultivar:
      1) Compute per-leaf mean angle: mean(angle_deg) for each leaf_no
      2) mean_angle_profile = mean(per-leaf means)
      3) angle_sd_profile   = sd(per-leaf means)
      4) angle_slope_leafno = slope of per-leaf means vs leaf_no
      5) GH_FSP35           = mean GH_FSP35 across rows for that cultivar
    """
    rows = []
    for cultivar, g in df.groupby("cultivar", dropna=True):
        g = g.copy()
        g["leaf_no"] = pd.to_numeric(g["leaf_no"], errors="coerce")
        g["angle_deg"] = pd.to_numeric(g["angle_deg"], errors="coerce")
        g["GH_FSP35"] = pd.to_numeric(g["GH_FSP35"], errors="coerce")

        g = g.dropna(subset=["leaf_no", "angle_deg"])
        if g.empty:
            continue

        leaf_means = g.groupby("leaf_no")["angle_deg"].mean().sort_index()

        if leaf_means.shape[0] < 2:
            slope = np.nan
        else:
            x = leaf_means.index.to_numpy(dtype=float)
            y = leaf_means.to_numpy(dtype=float)
            slope = np.polyfit(x, y, 1)[0]

        gh_vals = g["GH_FSP35"].dropna()
        gh_mean = float(gh_vals.mean()) if len(gh_vals) > 0 else np.nan

        rows.append(
            {
                "cultivar": str(cultivar),
                "mean_angle_profile": float(leaf_means.mean()),
                "angle_sd_profile": float(leaf_means.std(ddof=1)),
                "angle_slope_leafno": float(slope),
                "GH_FSP35": gh_mean,
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--data-dir",
        type=str,
        default=".",
        help="Directory containing 8leaf_final.xlsx",
    )
    ap.add_argument(
        "--arch-file",
        type=str,
        default="8leaf_final.xlsx",
        help="Filename of the architecture/GH dataset",
    )
    args = ap.parse_args()

    data_dir = Path(args.data_dir).expanduser().resolve()
    arch_path = data_dir / args.arch_file

    if not arch_path.exists():
        raise FileNotFoundError(f"Missing input file: {arch_path}")

    df = pd.read_excel(arch_path)

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {arch_path.name}: {missing}")

    feats = build_architecture_features(df)

    # Complete-case: cultivars with GH_FSP35 and all predictors non-missing
    base_cols = ["mean_angle_profile", "angle_sd_profile", "angle_slope_leafno", "GH_FSP35"]
    feats_cc = feats.dropna(subset=base_cols).reset_index(drop=True)

    n = feats_cc.shape[0]
    if n < 3:
        raise ValueError(f"Not enough complete-case cultivars to fit models (n={n}).")

    y = feats_cc["GH_FSP35"].to_numpy(dtype=float)

    # Model A: slope + slope^2
    X_a = feats_cc[["angle_slope_leafno"]].to_numpy(dtype=float)
    X_a = np.column_stack([X_a[:, 0], X_a[:, 0] ** 2])
    m_a = LinearRegression().fit(X_a, y)
    r2_a = r2_score(y, m_a.predict(X_a))

    # Model B: mean + slope + sd + slope^2
    mean_ = feats_cc["mean_angle_profile"].to_numpy(dtype=float)
    slope_ = feats_cc["angle_slope_leafno"].to_numpy(dtype=float)
    sd_ = feats_cc["angle_sd_profile"].to_numpy(dtype=float)
    X_b = np.column_stack([mean_, slope_, sd_, slope_ ** 2])
    m_b = LinearRegression().fit(X_b, y)
    r2_b = r2_score(y, m_b.predict(X_b))

    print(f"Input file: {arch_path}")
    print(f"Complete-case cultivars (n): {n}")
    print(f"R2_slope_quad (slope + slope^2): {r2_a:.3f}")
    print(f"R2_full_quad  (mean + slope + sd + slope^2): {r2_b:.3f}")

    # Optional: show which cultivars were used
    used = ", ".join(feats_cc["cultivar"].astype(str).tolist())
    print(f"Cultivars used: {used}")


if __name__ == "__main__":
    main()
