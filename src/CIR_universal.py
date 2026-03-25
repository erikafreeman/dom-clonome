"""
Universal CIR Pipeline
======================
Computes Constrained Isomeric Richness for ANY FT-ICR-MS formula list.
Input: CSV with at minimum columns C, H, O (and optionally N, S, P, sample intensities).
Output: formula_level_CIR.csv, summary stats, diagnostic plots.

Usage:
  python CIR_universal.py <input.csv> [--formula-col MolecularFormula] [--output-dir ./output]

The script auto-detects:
  - Whether element columns (C, H, O, N, S, P) are present or need parsing from a formula string
  - Whether sample intensity columns are present (wide format) or it's just a formula list
  - Column names for common formats (WHONDRS, CoreMS, Formularity, custom)
"""
import os, sys, re, argparse, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")
NL = chr(10)


# ═══════════════════════════════════════════════════════════════
# CORE FUNCTIONS (identical to main pipeline)
# ═══════════════════════════════════════════════════════════════

def parse_formula_string(formula_str):
    """Parse 'C15H20O7N2' into element dict."""
    elements = {"C": 0, "H": 0, "O": 0, "N": 0, "S": 0, "P": 0}
    for match in re.finditer(r"([A-Z][a-z]?)(\d*)", str(formula_str)):
        el, count = match.groups()
        if el in elements:
            elements[el] = int(count) if count else 1
    return elements


def calculate_indices(C, H, O, N, S, P):
    """Vectorised index calculation."""
    DBE = 1 + C - (H / 2) + (N / 2)
    HC = np.where(C > 0, H / C, np.nan)
    OC = np.where(C > 0, O / C, np.nan)
    NC = np.where(C > 0, N / C, np.nan)
    ai_num = 1 + C - (O / 2) - S - (H / 2)
    ai_den = C - (O / 2) - S - N - P
    AI_mod = np.where(ai_den > 0, np.maximum(0, ai_num / ai_den), 0)
    NOSC = np.where(C > 0, -((4*C + H - 3*N - 2*O + 5*P - 2*S) / C) + 4, np.nan)
    return DBE, HC, OC, NC, AI_mod, NOSC


def assign_classes(HC, OC, AI_mod, N):
    """Vectorised compound class assignment."""
    cc = np.full(len(HC), "other", dtype=object)
    cc = np.where(AI_mod >= 0.66, "condensed_aromatic", cc)
    cc = np.where((AI_mod >= 0.50) & (AI_mod < 0.66), "aromatic", cc)
    mask_cram = (HC < 1.5) & (OC < 0.5) & (AI_mod < 0.50)
    cc = np.where(mask_cram, "CRAM", cc)
    mask_carb = (HC >= 1.5) & (OC >= 0.5) & (AI_mod < 0.50)
    cc = np.where(mask_carb, "carbohydrate_like", cc)
    mask_lip = (HC >= 1.5) & (OC < 0.5) & (AI_mod < 0.50)
    cc = np.where(mask_lip, "lipid_like", cc)
    mask_prot = (N > 0) & (HC > 1.3) & (AI_mod < 0.50) & (HC < 1.5 if False else True)
    # protein_like: N>0, HC>1.3, not already assigned
    already = np.isin(cc, ["condensed_aromatic", "aromatic", "CRAM", "carbohydrate_like", "lipid_like"])
    mask_prot = (N > 0) & (HC > 1.3) & ~already
    cc = np.where(mask_prot, "protein_like", cc)
    mask_ncon = (N > 0) & ~np.isin(cc, ["condensed_aromatic", "aromatic", "CRAM",
                 "carbohydrate_like", "lipid_like", "protein_like"])
    cc = np.where(mask_ncon, "N_containing", cc)
    return cc


def get_constraint_counts(cc, DBE, OC, O, N, AI_mod):
    """Vectorised constraint penalty calculation."""
    n_forbidden = np.full(len(cc), 4, dtype=float)  # universal: 4 forbidden
    n_required = np.zeros(len(cc))
    max_rings = np.full(len(cc), np.nan)
    allow_aromatic = np.zeros(len(cc), dtype=bool)

    # condensed_aromatic
    m = cc == "condensed_aromatic"
    allow_aromatic[m] = True
    n_required[m] = 1
    max_rings[m] = np.floor(DBE[m] * 0.8)

    # aromatic
    m = cc == "aromatic"
    allow_aromatic[m] = True
    n_required[m] = 1
    n_required[m & (OC > 0.3)] += 1
    max_rings[m] = np.maximum(1, np.floor(DBE[m] * 0.6))

    # CRAM
    m = cc == "CRAM"
    n_required[m & (O >= 2)] = 1
    max_rings[m] = np.minimum(3, np.floor(DBE[m] * 0.7))

    # carbohydrate_like
    m = cc == "carbohydrate_like"
    n_required[m] = 1
    n_required[m & (DBE >= 1)] += 1
    n_forbidden[m] += 1  # double_bond_CC
    max_rings[m] = 2

    # lipid_like
    m = cc == "lipid_like"
    n_required[m & (OC < 0.2)] = 1
    max_rings[m] = np.maximum(0, np.floor(DBE[m] * 0.5))

    # protein_like
    m = cc == "protein_like"
    n_required[m & (N > 0)] = 1
    max_rings[m] = 1

    # N_containing
    m = cc == "N_containing"
    n_required[m] = 1

    return n_forbidden, n_required, max_rings, allow_aromatic


def estimate_CIR_vectorised(C, O, N, DBE, AI_mod, cc,
                             n_forbidden, n_required, max_rings, allow_aromatic):
    """Fully vectorised CIR estimation."""
    carbon_complexity = np.log10(np.maximum(1, C - 3)) * 0.8
    dbe_contribution = np.maximum(0.1, -0.15 * (DBE - 4)**2 + 1.2)
    oxygen_contribution = np.where(
        cc == "carbohydrate_like",
        np.log10(np.maximum(1, O)) * 0.4,
        np.log10(np.maximum(1, O)) * 0.7
    )
    nitrogen_contribution = np.log10(np.maximum(1, N + 1)) * 0.5
    required_penalty = n_required * 0.12
    forbidden_penalty = n_forbidden * 0.05
    aromaticity_penalty = np.where(allow_aromatic & (AI_mod >= 0.5), AI_mod * 0.8, 0)
    ring_penalty = np.where(
        ~np.isnan(max_rings),
        np.maximum(0, (DBE - max_rings)) * 0.1,
        0
    )
    log_CIR = (carbon_complexity + dbe_contribution + oxygen_contribution
               + nitrogen_contribution - required_penalty - forbidden_penalty
               - aromaticity_penalty - ring_penalty)
    log_CIR_cal = 0.78 + (log_CIR / 3.0) * 0.82
    return np.clip(log_CIR_cal, 0.3, 2.5)


# ═══════════════════════════════════════════════════════════════
# AUTO-DETECT AND LOAD
# ═══════════════════════════════════════════════════════════════

def load_and_detect(path, formula_col=None):
    """Load CSV and auto-detect format."""
    # Try different encodings
    for enc in ["utf-8", "latin-1", "cp1252"]:
        try:
            df = pd.read_csv(path, encoding=enc)
            break
        except UnicodeDecodeError:
            continue

    print(f"Loaded: {len(df)} rows x {len(df.columns)} columns")
    print(f"Columns: {list(df.columns[:20])}{'...' if len(df.columns) > 20 else ''}")

    # Auto-detect element columns
    el_cols = {}
    for el in ["C", "H", "O", "N", "S", "P"]:
        candidates = [c for c in df.columns if c.strip().upper() == el or c.strip() == el]
        if candidates:
            el_cols[el] = candidates[0]

    has_elements = len(el_cols) >= 3  # at least C, H, O

    if has_elements:
        print(f"Element columns found: {el_cols}")
        for el in ["C", "H", "O", "N", "S", "P"]:
            col = el_cols.get(el)
            if col:
                df[el] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
            else:
                df[el] = 0
    else:
        # Try to find a formula string column
        fc = formula_col
        if fc is None:
            for candidate in ["MolecularFormula", "Molecular Formula", "Formula",
                              "MF", "formula", "molecular_formula", "Mol_Formula"]:
                if candidate in df.columns:
                    fc = candidate
                    break
        if fc and fc in df.columns:
            print(f"Parsing formulae from column: {fc}")
            parsed = df[fc].apply(parse_formula_string)
            for el in ["C", "H", "O", "N", "S", "P"]:
                df[el] = parsed.apply(lambda x: x[el])
        else:
            print("ERROR: Cannot find element columns or formula column.")
            print(f"Available columns: {list(df.columns)}")
            sys.exit(1)

    # Filter to assigned formulae
    df = df[df["C"] > 0].copy()
    print(f"Assigned formulae (C > 0): {len(df)}")

    # Build formula string
    def make_formula(row):
        parts = []
        for el in ["C", "H", "O", "N", "S", "P"]:
            if row[el] > 0:
                parts.append(f"{el}{row[el]}" if row[el] > 1 else el)
        return "".join(parts)
    df["MolecularFormula"] = df.apply(make_formula, axis=1)

    # Detect sample columns (numeric columns not in the element/metadata set)
    meta_names = {"Mass", "C", "H", "O", "N", "S", "P", "C13", "Na",
                  "El_comp", "NeutralMass", "Error_ppm", "Candidates",
                  "MolecularFormula", "Formula", "MF", "Molecular Formula",
                  "DBE", "HC", "OC", "AI_mod", "NOSC", "compound_class",
                  "log_CIR", "m/z", "mz", "Intensity", "Class"}
    sample_cols = []
    for c in df.columns:
        if c in meta_names or c in ["C", "H", "O", "N", "S", "P"]:
            continue
        if df[c].dtype in [np.float64, np.int64, float, int]:
            # Check if it looks like intensity data (non-negative, some zeros)
            if df[c].min() >= 0 and (df[c] == 0).mean() > 0.1:
                sample_cols.append(c)

    print(f"Sample intensity columns detected: {len(sample_cols)}")
    return df, sample_cols


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def run(input_path, formula_col=None, output_dir=None, label=None):
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(input_path), "CIR_output")
    os.makedirs(output_dir, exist_ok=True)

    if label is None:
        label = os.path.splitext(os.path.basename(input_path))[0]

    print("=" * 70)
    print(f"CIR UNIVERSAL PIPELINE: {label}")
    print("=" * 70)

    # Load
    df, sample_cols = load_and_detect(input_path, formula_col)

    # Compute indices (vectorised)
    print("\nComputing indices...")
    C = df["C"].values.astype(float)
    H = df["H"].values.astype(float)
    O = df["O"].values.astype(float)
    N = df["N"].values.astype(float)
    S = df["S"].values.astype(float)
    P = df["P"].values.astype(float)

    DBE, HC, OC, NC, AI_mod, NOSC = calculate_indices(C, H, O, N, S, P)
    df["DBE"] = DBE
    df["HC"] = HC
    df["OC"] = OC
    df["NC"] = NC
    df["AI_mod"] = AI_mod
    df["NOSC"] = NOSC

    # Assign classes
    cc = assign_classes(HC, OC, AI_mod, N)
    df["compound_class"] = cc

    print("Compound classes:")
    for cls, n in pd.Series(cc).value_counts().items():
        print(f"  {cls:25s} {n:6d} ({100*n/len(df):.1f}%)")

    # Compute CIR (vectorised)
    print("\nEstimating CIR...")
    n_forbidden, n_required, max_rings, allow_aromatic = get_constraint_counts(
        cc, DBE, OC, O, N, AI_mod)

    log_CIR = estimate_CIR_vectorised(C, O, N, DBE, AI_mod, cc,
                                       n_forbidden, n_required, max_rings, allow_aromatic)
    # Set NaN for C=0
    log_CIR[C == 0] = np.nan
    df["log_CIR"] = log_CIR

    # Summary
    print("\nCIR by class:")
    summary = df.groupby("compound_class")["log_CIR"].describe().round(3)
    print(summary.to_string())

    # Validation
    cir_med = df.groupby("compound_class")["log_CIR"].median()
    if "condensed_aromatic" in cir_med and "CRAM" in cir_med:
        ok = cir_med["condensed_aromatic"] < cir_med["CRAM"]
        print(f"\nValidation: condensed_aromatic < CRAM: {'PASS' if ok else 'WARN'}")

    rng = df["log_CIR"].dropna()
    print(f"Range: {rng.min():.2f} - {rng.max():.2f}")
    print(f"NaN: {df['log_CIR'].isna().sum()}")

    # Save formula-level output
    out_cols = [c for c in ["MolecularFormula", "Mass", "C", "H", "O", "N", "S", "P",
                             "DBE", "HC", "OC", "NC", "AI_mod", "NOSC",
                             "compound_class", "log_CIR"] if c in df.columns]
    out_path = os.path.join(output_dir, f"{label}_formula_CIR.csv")
    df[out_cols].to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")

    # Sample-level aggregation if samples exist
    if len(sample_cols) > 0:
        print(f"\nAggregating {len(sample_cols)} samples...")
        sample_rows = []
        for col in sample_cols:
            intensities = df[col]
            detected = intensities > 0
            sf = df.loc[detected]
            si = intensities[detected]
            if len(sf) == 0:
                continue
            valid = sf["log_CIR"].notna()
            if valid.sum() == 0:
                continue
            w = si[valid] / si[valid].sum()
            row = {
                "sample": col,
                "n_formulas": len(sf),
                "mean_CIR": np.average(sf.loc[valid, "log_CIR"].values, weights=w.values),
                "sd_CIR": sf["log_CIR"].std(),
                "richness": len(sf),
            }
            for k, v in sf.groupby("compound_class")["log_CIR"].mean().items():
                row[f"CIR_{k}"] = v
            sample_rows.append(row)
        df_s = pd.DataFrame(sample_rows)
        s_path = os.path.join(output_dir, f"{label}_sample_CIR.csv")
        df_s.to_csv(s_path, index=False)
        print(f"Saved: {s_path} ({len(df_s)} samples)")
    else:
        df_s = None

    # Plot
    print("\nPlotting...")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    order = df.groupby("compound_class")["log_CIR"].median().sort_values().index
    sns.boxplot(data=df, x="compound_class", y="log_CIR", order=order, ax=ax, palette="Set2")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("log10(CIR)")
    ax.set_title(f"{label}: CIR by Class")

    ax = axes[1]
    v = df.dropna(subset=["OC", "HC", "log_CIR"])
    sc = ax.scatter(v["OC"], v["HC"], c=v["log_CIR"], cmap="plasma", alpha=0.3, s=3, rasterized=True)
    plt.colorbar(sc, ax=ax, label="log10(CIR)", shrink=0.8)
    ax.set_xlabel("O/C"); ax.set_ylabel("H/C")
    ax.set_xlim(0, 1.2); ax.set_ylim(0, 2.5)
    ax.set_title(f"{label}: Van Krevelen")

    ax = axes[2]
    if df_s is not None and len(df_s) > 1:
        ax.hist(df_s["mean_CIR"].dropna(), bins=min(25, len(df_s)), edgecolor="black",
                color="steelblue", alpha=0.8)
        ax.set_xlabel("Mean log10(CIR)")
        ax.set_ylabel("Samples")
        ax.set_title(f"{label}: Sample CIR")
    else:
        ax.hist(df["log_CIR"].dropna(), bins=50, edgecolor="black", color="steelblue", alpha=0.8)
        ax.set_xlabel("log10(CIR)")
        ax.set_ylabel("Formulae")
        ax.set_title(f"{label}: Formula CIR Distribution")

    plt.tight_layout()
    fig_path = os.path.join(output_dir, f"{label}_CIR.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig_path}")

    # Summary text
    txt_path = os.path.join(output_dir, f"{label}_CIR_summary.txt")
    with open(txt_path, "w") as f:
        f.write(f"CIR SUMMARY: {label}" + NL + "=" * 60 + NL)
        f.write(f"Formulae: {len(df)}" + NL)
        f.write(f"Samples: {len(sample_cols)}" + NL + NL)
        f.write(summary.to_string() + NL)
    print(f"Saved: {txt_path}")

    print("\n" + "=" * 70)
    print(f"DONE: {label}")
    print("=" * 70)

    return df, df_s


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Universal CIR Pipeline")
    parser.add_argument("input", help="Path to FT-ICR-MS CSV")
    parser.add_argument("--formula-col", default=None, help="Column name for formula strings")
    parser.add_argument("--output-dir", default=None, help="Output directory")
    parser.add_argument("--label", default=None, help="Dataset label for output files")
    args = parser.parse_args()
    run(args.input, args.formula_col, args.output_dir, args.label)
