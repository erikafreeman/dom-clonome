"""
Robustness checks requested by Gemini review.
Q7:  Does CIR-richness anti-correlation survive replicate filtering?
Q9:  Does class ordering persist under null (uniform) constraints?
Q10: Does CIR_CRAM outperform simpler metrics (mean DBE, mean MW)?
+    Outlier leverage check on CRAM-latitude correlation.
"""
import os, warnings
import pandas as pd
import numpy as np
from scipy import stats
warnings.filterwarnings("ignore")

BASE = os.path.expanduser(os.path.join("~", "OneDrive", "Desktop", "ISOMERs"))
OUT = os.path.join(BASE, "output")

# Load pipeline outputs
df_f = pd.read_csv(os.path.join(OUT, "formula_level_CIR.csv"))
df_s = pd.read_csv(os.path.join(OUT, "sample_level_CIR.csv"))

# Load raw FTICR data for replicate filtering
DATA_DIR = os.path.join(BASE, "WHONDRS_SWMB_DataPackage_2018-11_Rev1", "data")
FTICR_PATH = os.path.join(DATA_DIR, "WHONDRS_SW_Kits_FTICR", "WHONDRS_SW_Kits_FTICR_Aug2018.csv")
df_raw = pd.read_csv(FTICR_PATH)

meta_cols = ["Mass","C","H","O","N","C13","S","P","Na","El_comp","NeutralMass","Error_ppm","Candidates"]
sample_cols = [c for c in df_raw.columns if c not in meta_cols]

# Filter to assigned formulae
df_raw = df_raw[df_raw["C"] > 0].copy()

NL = "\n"
results = []
results.append("=" * 70)
results.append("ROBUSTNESS CHECKS — Response to Gemini Review")
results.append("=" * 70)

# ═══════════════════════════════════════════════════════════════
# Q7: Replicate-filtered CIR-richness correlation
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("Q7: Does CIR-richness anti-correlation survive replicate filtering?")
results.append("-" * 70)

# Group sample columns by site
site_samples = {}
for col in sample_cols:
    import re
    m = re.match(r"(S\d+)\.(\d+)", col)
    if m:
        site = m.group(1)
        if site not in site_samples:
            site_samples[site] = []
        site_samples[site].append(col)

# For each site, find formulae present in ALL replicates
# Then recompute mean_CIR and richness using only those formulae
# Use MolecularFormula as the join key (indices differ between raw and formula CSVs)

# Build a formula-to-CIR lookup from the formula-level output
formula_cir = df_f.set_index("MolecularFormula")["log_CIR"].to_dict()

# Build formula strings for raw data (same logic as pipeline)
def make_formula(row):
    parts = []
    for el in ["C","H","O","N","S","P"]:
        v = int(row.get(el, 0)) if pd.notna(row.get(el, 0)) else 0
        if v > 0:
            parts.append(f"{el}{v}" if v > 1 else el)
    return "".join(parts)

df_raw["MolecularFormula"] = df_raw.apply(make_formula, axis=1)

robust_rows = []
for site, cols in site_samples.items():
    if len(cols) < 2:
        continue
    # Formula is "replicate-confirmed" if intensity > 0 in ALL replicates
    present_in_all = (df_raw[cols] > 0).all(axis=1)
    confirmed = df_raw.loc[present_in_all]

    for col in cols:
        intensities = confirmed[col]
        detected = intensities > 0
        if detected.sum() == 0:
            continue
        det_formulas = confirmed.loc[detected, "MolecularFormula"]
        det_intensities = intensities[detected]
        # Look up CIR for each formula
        cir_vals = det_formulas.map(formula_cir)
        valid = cir_vals.notna()
        if valid.sum() == 0:
            continue
        w = det_intensities[valid] / det_intensities[valid].sum()
        mean_cir = np.average(cir_vals[valid].values, weights=w.values)
        robust_rows.append({
            "sample": col,
            "site_id": site,
            "mean_CIR_robust": mean_cir,
            "richness_robust": detected.sum(),
        })

df_robust = pd.DataFrame(robust_rows)

if len(df_robust) > 5:
    r_orig, p_orig = stats.pearsonr(df_s["richness"], df_s["mean_CIR"])
    r_robust, p_robust = stats.pearsonr(df_robust["richness_robust"], df_robust["mean_CIR_robust"])
    results.append(f"  Original (all formulae):          r = {r_orig:.3f}, p = {p_orig:.4f}, n = {len(df_s)}")
    results.append(f"  Replicate-filtered (all-rep only): r = {r_robust:.3f}, p = {p_robust:.4f}, n = {len(df_robust)}")
    results.append(f"  Mean formulae per sample (original):  {df_s['richness'].mean():.0f}")
    results.append(f"  Mean formulae per sample (filtered):  {df_robust['richness_robust'].mean():.0f}")
    if abs(r_robust) > 0.2 and p_robust < 0.05:
        results.append("  VERDICT: Anti-correlation SURVIVES replicate filtering. Not an artifact of noise.")
    elif p_robust >= 0.05:
        results.append("  VERDICT: Anti-correlation WEAKENS and loses significance after filtering.")
        results.append("           This suggests it may be partly driven by noisy low-intensity peaks.")
    else:
        results.append(f"  VERDICT: Correlation changes from {r_orig:.3f} to {r_robust:.3f}.")

# ═══════════════════════════════════════════════════════════════
# Q9: Null penalty model — uniform constraints across all classes
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("Q9: Does class ordering persist under null (uniform) constraints?")
results.append("-" * 70)

def estimate_CIR_null(row):
    """CIR with ZERO class-specific constraints — only universal forbidden list."""
    C, O, N, DBE = row["C"], row["O"], row["N"], row["DBE"]
    if C == 0 or np.isnan(DBE):
        return np.nan
    carbon_complexity = np.log10(max(1, C - 3)) * 0.8
    dbe_contribution = max(0.1, -0.15 * (DBE - 4)**2 + 1.2)
    oxygen_contribution = np.log10(max(1, O)) * 0.7  # same for all classes
    nitrogen_contribution = np.log10(max(1, N + 1)) * 0.5
    # Only universal penalties: 4 forbidden, 0 required, no aromatic/ring penalties
    forbidden_penalty = 4 * 0.05
    log_CIR = (carbon_complexity + dbe_contribution + oxygen_contribution
               + nitrogen_contribution - forbidden_penalty)
    return np.clip(0.78 + (log_CIR / 3.0) * 0.82, 0.3, 2.5)

df_f["log_CIR_null"] = df_f.apply(estimate_CIR_null, axis=1)

results.append("  Median CIR by class — Original vs. Null (no class-specific constraints):")
results.append(f"  {'Class':<25s} {'Original':>10s} {'Null':>10s} {'Diff':>10s}")
results.append(f"  {'-'*25} {'-'*10} {'-'*10} {'-'*10}")

class_order_orig = df_f.groupby("compound_class")["log_CIR"].median().sort_values()
class_order_null = df_f.groupby("compound_class")["log_CIR_null"].median().sort_values()

for cc in class_order_orig.index:
    orig = class_order_orig[cc]
    null = class_order_null.get(cc, np.nan)
    diff = null - orig
    results.append(f"  {cc:<25s} {orig:>10.3f} {null:>10.3f} {diff:>+10.3f}")

# Check if rank order is preserved
orig_rank = list(class_order_orig.index)
null_rank = list(class_order_null.index)
rank_preserved = orig_rank == null_rank

results.append(f"")
results.append(f"  Original rank order:  {' < '.join(orig_rank)}")
results.append(f"  Null rank order:      {' < '.join(null_rank)}")
results.append(f"  Rank order preserved: {'YES' if rank_preserved else 'NO — ordering changed'}")

if rank_preserved:
    results.append("  VERDICT: Class ordering is driven by UNDERLYING CHEMISTRY (C, O, DBE distributions),")
    results.append("           not by asymmetric constraint penalties. Constraints modulate magnitude,")
    results.append("           not the fundamental ranking.")
else:
    results.append("  VERDICT: Constraint penalties CHANGE the class ordering.")
    results.append("           The ranking is partly driven by asymmetric penalty counts.")

# ═══════════════════════════════════════════════════════════════
# Q10: Does CIR_CRAM outperform simpler metrics?
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("Q10: Does CIR_CRAM outperform simpler metrics for latitude prediction?")
results.append("-" * 70)

# Compute simpler CRAM metrics per sample
# Use MolecularFormula join since indices differ
df_f_indexed = df_f.set_index("MolecularFormula")
simple_rows = []
for _, srow in df_s.iterrows():
    sample = srow["sample"]
    if sample not in df_raw.columns:
        continue
    intensities = df_raw[sample]
    detected_formulas = df_raw.loc[intensities > 0, "MolecularFormula"]
    sf = df_f_indexed.loc[df_f_indexed.index.isin(detected_formulas)]
    cram_mask = sf["compound_class"] == "CRAM"
    cram = sf[cram_mask]
    if len(cram) < 5:
        continue
    simple_rows.append({
        "sample": sample,
        "site_id": srow.get("site_id", ""),
        "latitude": srow.get("latitude", np.nan),
        "mean_DBE_CRAM": cram["DBE"].mean(),
        "mean_HC_CRAM": cram["HC"].mean(),
        "mean_OC_CRAM": cram["OC"].mean(),
        "mean_AI_CRAM": cram["AI_mod"].mean(),
        "mean_Mass_CRAM": cram["Mass"].mean() if "Mass" in cram.columns else np.nan,
        "mean_NOSC_CRAM": cram["NOSC"].mean(),
        "n_CRAM": len(cram),
        "CIR_CRAM": srow.get("CIR_CRAM", np.nan),
    })

df_simple = pd.DataFrame(simple_rows)

results.append(f"  {'Metric':<25s} {'r vs Lat':>10s} {'p':>12s} {'Conclusion':>20s}")
results.append(f"  {'-'*25} {'-'*10} {'-'*12} {'-'*20}")

metrics_to_test = [
    ("CIR_CRAM", "CIR_CRAM"),
    ("mean_DBE_CRAM", "Mean DBE (CRAM)"),
    ("mean_HC_CRAM", "Mean H/C (CRAM)"),
    ("mean_OC_CRAM", "Mean O/C (CRAM)"),
    ("mean_AI_CRAM", "Mean AI_mod (CRAM)"),
    ("mean_Mass_CRAM", "Mean Mass (CRAM)"),
    ("mean_NOSC_CRAM", "Mean NOSC (CRAM)"),
    ("n_CRAM", "CRAM count"),
]

for col, label in metrics_to_test:
    m = df_simple[["latitude", col]].notna().all(axis=1)
    if m.sum() > 5:
        r, p = stats.pearsonr(df_simple.loc[m, "latitude"], df_simple.loc[m, col])
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        results.append(f"  {label:<25s} {r:>10.3f} {p:>12.4f} {sig:>20s}")

results.append("")
results.append("  If CIR_CRAM has the strongest |r|, it outperforms simpler metrics.")
results.append("  If a simpler metric matches or exceeds CIR_CRAM, CIR adds no value.")

# ═══════════════════════════════════════════════════════════════
# Outlier leverage check on CRAM-latitude
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("Outlier leverage check: CRAM-latitude correlation")
results.append("-" * 70)

m = df_s[["latitude", "CIR_CRAM"]].notna().all(axis=1)
if m.sum() > 5:
    lat = df_s.loc[m, "latitude"].values
    cir = df_s.loc[m, "CIR_CRAM"].values
    sites = df_s.loc[m, "site_id"].values

    r_full, p_full = stats.pearsonr(lat, cir)
    results.append(f"  Full dataset: r = {r_full:.3f}, p = {p_full:.6f}, n = {m.sum()}")

    # Leave-one-site-out
    unique_sites = df_s.loc[m, "site_id"].unique()
    loo_results = []
    for site in unique_sites:
        mask = df_s.loc[m, "site_id"] != site
        r_loo, _ = stats.pearsonr(lat[mask.values], cir[mask.values])
        loo_results.append((site, r_loo, r_full - r_loo))

    loo_results.sort(key=lambda x: abs(x[2]), reverse=True)
    results.append(f"  Leave-one-site-out analysis ({len(unique_sites)} sites):")
    results.append(f"    {'Site':<12s} {'r without':>12s} {'Change':>10s}")
    for site, r_loo, delta in loo_results[:10]:
        results.append(f"    {site:<12s} {r_loo:>12.3f} {delta:>+10.3f}")

    r_values = [x[1] for x in loo_results]
    results.append(f"  LOO r range: [{min(r_values):.3f}, {max(r_values):.3f}]")
    results.append(f"  Max influence: removing {loo_results[0][0]} changes r by {loo_results[0][2]:+.3f}")

    if max(r_values) < -0.5:
        results.append("  VERDICT: Correlation ROBUST. Even removing the most influential site,")
        results.append(f"           r remains strong (worst case: {max(r_values):.3f}).")
    else:
        results.append("  VERDICT: Correlation may be driven by specific sites. Investigate.")

# ═══════════════════════════════════════════════════════════════
# OUTPUT
# ═══════════════════════════════════════════════════════════════
output = NL.join(results)
print(output)

out_path = os.path.join(OUT, "robustness_checks.txt")
with open(out_path, "w") as f:
    f.write(output + NL)
print(f"\nSaved: {out_path}")
