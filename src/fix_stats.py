"""
Fix pseudoreplication and verify all numbers for manuscript.
Gemini review identified: correlations at n=80 triplicates are pseudoreplicated
because latitude is a site-level variable. Must use site means (n=27).
"""
import os, warnings
import pandas as pd
import numpy as np
from scipy import stats
warnings.filterwarnings("ignore")

BASE = os.path.expanduser(os.path.join("~", "OneDrive", "Desktop", "ISOMERs"))
OUT = os.path.join(BASE, "output")
NL = "\n"

df_s = pd.read_csv(os.path.join(OUT, "sample_level_CIR.csv"))
df_f = pd.read_csv(os.path.join(OUT, "formula_level_CIR.csv"))

results = []
results.append("=" * 70)
results.append("STATISTICAL CORRECTIONS — Addressing Gemini Referee Report")
results.append("=" * 70)

# ═══════════════════════════════════════════════════════════════
# 1. CHECK: Table 1 formula counts
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("1. TABLE 1 FORMULA COUNTS")
results.append("-" * 70)
cc_counts = df_f["compound_class"].value_counts()
total = cc_counts.sum()
results.append(f"Total assigned formulae: {total}")
for cc, n in cc_counts.sort_values(ascending=False).items():
    results.append(f"  {cc:25s} {n:6d} ({100*n/total:.1f}%)")
results.append(f"  {'SUM':25s} {total:6d} ({cc_counts.sum()/total*100:.1f}%)")
results.append("")
results.append(f"  'other' class present: {'YES' if 'other' in cc_counts else 'NO'}")
results.append(f"  'other' count: {cc_counts.get('other', 0)}")
results.append(f"  Sum WITHOUT other: {total - cc_counts.get('other', 0)}")
results.append("  FIX: Table 1 must include 'other' class to sum to 8,972.")

# ═══════════════════════════════════════════════════════════════
# 2. CHECK: Class ordering by median CIR
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("2. CLASS ORDERING BY MEDIAN CIR")
results.append("-" * 70)
medians = df_f.groupby("compound_class")["log_CIR"].median().sort_values()
for cc, med in medians.items():
    isomers = 10**med
    results.append(f"  {cc:25s} median={med:.3f}  (~{isomers:.1f} isomers)")
results.append("")
results.append(f"  Lowest class: {medians.index[0]} ({medians.iloc[0]:.3f})")
results.append(f"  NOTE: aromatic (0.881) < condensed_aromatic (0.896)")
results.append("  FIX: Text should say 'aromatic and condensed aromatic are the most")
results.append("  constrained' rather than singling out condensed_aromatic as lowest.")

# ═══════════════════════════════════════════════════════════════
# 3. FIX PSEUDOREPLICATION: Site-level means
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("3. PSEUDOREPLICATION FIX: Site-level means (n=27)")
results.append("-" * 70)

# Aggregate triplicates to site means
site_cols = ["mean_CIR", "sd_CIR", "prop_high_CIR", "richness", "latitude",
             "longitude", "temperature_C", "DOC_mgL"]
cir_class_cols = [c for c in df_s.columns if c.startswith("CIR_")]
all_numeric = site_cols + cir_class_cols

df_site = df_s.groupby("site_id")[all_numeric].mean().reset_index()
results.append(f"  Sites: {len(df_site)}")
results.append(f"  (Aggregated from {len(df_s)} triplicates)")

results.append("")
results.append("  ORIGINAL (n=80 pseudoreplicated) vs CORRECTED (n=27 site means):")
results.append(f"  {'Variable':<25s} {'r(n=80)':>10s} {'p(n=80)':>12s} {'r(n=27)':>10s} {'p(n=27)':>12s} {'Change':>8s}")
results.append(f"  {'-'*25} {'-'*10} {'-'*12} {'-'*10} {'-'*12} {'-'*8}")

env_vars = [
    ("latitude", "mean_CIR", "Lat vs mean_CIR"),
    ("temperature_C", "mean_CIR", "Temp vs mean_CIR"),
    ("DOC_mgL", "mean_CIR", "DOC vs mean_CIR"),
    ("richness", "mean_CIR", "Richness vs mean_CIR"),
]

# Add CIR_CRAM if present
if "CIR_CRAM" in df_s.columns:
    env_vars.extend([
        ("latitude", "CIR_CRAM", "Lat vs CIR_CRAM"),
        ("temperature_C", "CIR_CRAM", "Temp vs CIR_CRAM"),
        ("DOC_mgL", "CIR_CRAM", "DOC vs CIR_CRAM"),
    ])

for xvar, yvar, label in env_vars:
    # Original (n=80)
    m80 = df_s[[xvar, yvar]].notna().all(axis=1)
    if m80.sum() > 5:
        r80, p80 = stats.pearsonr(df_s.loc[m80, xvar], df_s.loc[m80, yvar])
    else:
        r80, p80 = np.nan, np.nan

    # Corrected (n=27)
    m27 = df_site[[xvar, yvar]].notna().all(axis=1)
    if m27.sum() > 5:
        r27, p27 = stats.pearsonr(df_site.loc[m27, xvar], df_site.loc[m27, yvar])
    else:
        r27, p27 = np.nan, np.nan

    delta = r27 - r80 if not (np.isnan(r80) or np.isnan(r27)) else np.nan
    p80_s = f"{p80:.4f}" if p80 < 0.01 else f"{p80:.3f}"
    p27_s = f"{p27:.4f}" if p27 < 0.01 else f"{p27:.3f}"
    results.append(f"  {label:<25s} {r80:>10.3f} {p80_s:>12s} {r27:>10.3f} {p27_s:>12s} {delta:>+8.3f}")

# ═══════════════════════════════════════════════════════════════
# 4. Leave-one-out at SITE level
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("4. LEAVE-ONE-OUT (SITE LEVEL, n=27)")
results.append("-" * 70)
if "CIR_CRAM" in df_site.columns:
    m = df_site[["latitude", "CIR_CRAM"]].notna().all(axis=1)
    lat = df_site.loc[m, "latitude"].values
    cir = df_site.loc[m, "CIR_CRAM"].values
    sites = df_site.loc[m, "site_id"].values

    r_full, p_full = stats.pearsonr(lat, cir)
    results.append(f"  Full: r = {r_full:.3f}, p = {p_full:.6f}, n = {m.sum()}")

    loo = []
    for i in range(len(sites)):
        mask = np.ones(len(sites), dtype=bool)
        mask[i] = False
        r_loo, _ = stats.pearsonr(lat[mask], cir[mask])
        loo.append((sites[i], r_loo, r_full - r_loo))

    loo.sort(key=lambda x: abs(x[2]), reverse=True)
    results.append(f"  Top 5 most influential sites:")
    for site, r_loo, delta in loo[:5]:
        results.append(f"    {site}: r_without = {r_loo:.3f} (change: {delta:+.3f})")
    r_vals = [x[1] for x in loo]
    results.append(f"  LOO range: [{min(r_vals):.3f}, {max(r_vals):.3f}]")

# ═══════════════════════════════════════════════════════════════
# 5. Comparison with simpler metrics AT SITE LEVEL
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("5. CIR_CRAM vs SIMPLER METRICS (SITE LEVEL, n=27)")
results.append("-" * 70)

# Recompute simpler CRAM metrics at site level from sample data
# We need to average the CRAM-specific metrics across triplicates
DATA_DIR = os.path.join(BASE, "WHONDRS_SWMB_DataPackage_2018-11_Rev1", "data")
FTICR_PATH = os.path.join(DATA_DIR, "WHONDRS_SW_Kits_FTICR", "WHONDRS_SW_Kits_FTICR_Aug2018.csv")
df_raw = pd.read_csv(FTICR_PATH)
meta_cols = ["Mass","C","H","O","N","C13","S","P","Na","El_comp","NeutralMass","Error_ppm","Candidates"]
sample_cols = [c for c in df_raw.columns if c not in meta_cols]
df_raw = df_raw[df_raw["C"] > 0].copy()

# Build formula strings
import re
def make_formula(row):
    parts = []
    for el in ["C","H","O","N","S","P"]:
        v = int(row.get(el, 0)) if pd.notna(row.get(el, 0)) else 0
        if v > 0:
            parts.append(f"{el}{v}" if v > 1 else el)
    return "".join(parts)
df_raw["MolecularFormula"] = df_raw.apply(make_formula, axis=1)
df_f_idx = df_f.set_index("MolecularFormula")

simple_rows = []
for _, srow in df_s.iterrows():
    sample = srow["sample"]
    if sample not in df_raw.columns:
        continue
    intensities = df_raw[sample]
    det_formulas = df_raw.loc[intensities > 0, "MolecularFormula"]
    sf = df_f_idx.loc[df_f_idx.index.isin(det_formulas)]
    cram = sf[sf["compound_class"] == "CRAM"]
    if len(cram) < 5:
        continue
    simple_rows.append({
        "sample": sample, "site_id": srow.get("site_id", ""),
        "mean_DBE_CRAM": cram["DBE"].mean(),
        "mean_HC_CRAM": cram["HC"].mean(),
        "mean_OC_CRAM": cram["OC"].mean(),
        "mean_AI_CRAM": cram["AI_mod"].mean(),
        "mean_Mass_CRAM": cram["Mass"].mean() if "Mass" in cram.columns else np.nan,
        "CIR_CRAM": srow.get("CIR_CRAM", np.nan),
        "latitude": srow.get("latitude", np.nan),
    })

df_simp = pd.DataFrame(simple_rows)
# Aggregate to site level
df_simp_site = df_simp.groupby("site_id").mean(numeric_only=True).reset_index()

results.append(f"  {'Metric':<25s} {'r(n=27)':>10s} {'p':>12s}")
results.append(f"  {'-'*25} {'-'*10} {'-'*12}")
for col, label in [("CIR_CRAM", "CIR_CRAM"), ("mean_HC_CRAM", "Mean H/C"),
                    ("mean_DBE_CRAM", "Mean DBE"), ("mean_Mass_CRAM", "Mean Mass"),
                    ("mean_AI_CRAM", "Mean AI_mod")]:
    m = df_simp_site[["latitude", col]].notna().all(axis=1)
    if m.sum() > 5:
        r, p = stats.pearsonr(df_simp_site.loc[m, "latitude"], df_simp_site.loc[m, col])
        results.append(f"  {label:<25s} {r:>10.3f} {p:>12.4f}")

# ═══════════════════════════════════════════════════════════════
# 6. CHECK: 63C temperature
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("6. TEMPERATURE CHECK")
results.append("-" * 70)
temps = df_s[["site_id", "temperature_C"]].drop_duplicates()
temps_sorted = temps.sort_values("temperature_C", ascending=False).head(5)
for _, row in temps_sorted.iterrows():
    results.append(f"  {row['site_id']}: {row['temperature_C']:.1f} C")
max_t = df_s["temperature_C"].max()
if max_t > 50:
    results.append(f"  WARNING: Max temperature {max_t:.1f}C is unusually high.")
    results.append("  This may be a hot spring site or a data entry error. Verify.")

# ═══════════════════════════════════════════════════════════════
# 7. Replicate-filtered anti-correlation at SITE level
# ═══════════════════════════════════════════════════════════════
results.append("")
results.append("7. CIR-RICHNESS ANTI-CORRELATION AT SITE LEVEL")
results.append("-" * 70)
m = df_site[["richness", "mean_CIR"]].notna().all(axis=1)
if m.sum() > 5:
    r, p = stats.pearsonr(df_site.loc[m, "richness"], df_site.loc[m, "mean_CIR"])
    results.append(f"  Site-level (n=27): r = {r:.3f}, p = {p:.4f}")

# ═══════════════════════════════════════════════════════════════
# OUTPUT
# ═══════════════════════════════════════════════════════════════
output = NL.join(results)
print(output)

out_path = os.path.join(OUT, "statistical_corrections.txt")
with open(out_path, "w") as f:
    f.write(output + NL)
print(f"\nSaved: {out_path}")
