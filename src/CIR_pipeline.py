"""CIR Pipeline - WHONDRS FT-ICR-MS"""
import os, sys, re, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from tqdm import tqdm
warnings.filterwarnings("ignore")

BASE = os.path.expanduser(os.path.join("~", "OneDrive", "Desktop", "ISOMERs"))
DATA_DIR = os.path.join(BASE, "WHONDRS_SWMB_DataPackage_2018-11_Rev1", "data")
FTICR_PATH = os.path.join(DATA_DIR, "WHONDRS_SW_Kits_FTICR", "WHONDRS_SW_Kits_FTICR_Aug2018.csv")
META_PATH = os.path.join(DATA_DIR, "WHONDRS_SW_Kits_Metadata_20Nov2018.csv")
GEOCHEM_PATH = os.path.join(DATA_DIR, "WHONDRS_SW_Kits_Geochem", "WHONDRS_SW_Kits_Geochem_20Nov2018.csv")
OUT_DIR = os.path.join(BASE, "output")
os.makedirs(OUT_DIR, exist_ok=True)
NL = chr(10)

def calculate_indices(row):
    C, H, O, N, S, P = row["C"], row["H"], row["O"], row["N"], row["S"], row["P"]
    if C == 0:
        return pd.Series({"DBE": np.nan, "HC": np.nan, "OC": np.nan, "NC": np.nan,
                          "AI_mod": np.nan, "NOSC": np.nan, "compound_class": "unassigned"})
    DBE = 1 + C - (H / 2) + (N / 2)
    HC, OC, NC = H / C, O / C, N / C
    ai_num = 1 + C - (O / 2) - S - (H / 2)
    ai_den = C - (O / 2) - S - N - P
    AI_mod = max(0, ai_num / ai_den) if ai_den > 0 else 0
    NOSC = -((4*C + H - 3*N - 2*O + 5*P - 2*S) / C) + 4
    if AI_mod >= 0.66: cc = "condensed_aromatic"
    elif AI_mod >= 0.50: cc = "aromatic"
    elif HC < 1.5 and OC < 0.5: cc = "CRAM"
    elif HC >= 1.5 and OC >= 0.5: cc = "carbohydrate_like"
    elif HC >= 1.5 and OC < 0.5: cc = "lipid_like"
    elif N > 0 and HC > 1.3: cc = "protein_like"
    elif N > 0: cc = "N_containing"
    else: cc = "other"
    return pd.Series({"DBE": DBE, "HC": HC, "OC": OC, "NC": NC,
                      "AI_mod": AI_mod, "NOSC": NOSC, "compound_class": cc})

def get_constraints(row):
    cc, DBE, AI_mod, OC = row["compound_class"], row["DBE"], row["AI_mod"], row["OC"]
    O, N = row["O"], row["N"]
    c = {"forbidden": ["peroxide", "triple_bond", "ring_3", "ring_4"],
         "required": [], "max_rings": None, "allow_aromatic": False}
    if cc == "condensed_aromatic":
        c["allow_aromatic"] = True
        c["required"].append("aromatic_ring")
        c["max_rings"] = int(DBE * 0.8)
    elif cc == "aromatic":
        c["allow_aromatic"] = True
        c["required"].append("aromatic_ring")
        if OC > 0.3: c["required"].append("hydroxyl_or_methoxy")
        c["max_rings"] = max(1, int(DBE * 0.6))
    elif cc == "CRAM":
        if O >= 2: c["required"].append("carboxyl_group")
        c["max_rings"] = min(3, int(DBE * 0.7))
    elif cc == "carbohydrate_like":
        c["required"].append("hydroxyl_group")
        if DBE >= 1: c["required"].append("ring_oxygen")
        c["max_rings"] = 2
        c["forbidden"].append("double_bond_CC")
    elif cc == "lipid_like":
        if OC < 0.2: c["required"].append("carboxyl_or_ester")
        c["max_rings"] = max(0, int(DBE * 0.5))
    elif cc == "protein_like":
        if N > 0: c["required"].append("amine_or_amide")
        c["max_rings"] = 1
    elif cc == "N_containing":
        c["required"].append("nitrogen_heterocycle_or_amine")
    return c

def estimate_CIR(row, constraints):
    C, O, N, DBE = row["C"], row["O"], row["N"], row["DBE"]
    AI_mod, cc = row["AI_mod"], row["compound_class"]
    if C == 0 or np.isnan(DBE): return np.nan
    carbon_complexity = np.log10(max(1, C - 3)) * 0.8
    dbe_contribution = max(0.1, -0.15 * (DBE - 4)**2 + 1.2)
    if cc == "carbohydrate_like":
        oxygen_contribution = np.log10(max(1, O)) * 0.4
    else:
        oxygen_contribution = np.log10(max(1, O)) * 0.7
    nitrogen_contribution = np.log10(max(1, N + 1)) * 0.5
    required_penalty = len(constraints["required"]) * 0.12
    forbidden_penalty = len(constraints["forbidden"]) * 0.05
    aromaticity_penalty = AI_mod * 0.8 if (constraints["allow_aromatic"] and AI_mod >= 0.5) else 0
    ring_penalty = max(0, (DBE - constraints["max_rings"])) * 0.1 if constraints["max_rings"] is not None else 0
    log_CIR = (carbon_complexity + dbe_contribution + oxygen_contribution
               + nitrogen_contribution - required_penalty - forbidden_penalty
               - aromaticity_penalty - ring_penalty)
    return np.clip(0.78 + (log_CIR / 3.0) * 0.82, 0.3, 2.5)

def aggregate_to_sample(df, sample_cols):
    results = []
    for sample in tqdm(sample_cols, desc="Aggregating"):
        intensities = df[sample]
        detected = intensities > 0
        sf = df.loc[detected]
        si = intensities[detected]
        if len(sf) == 0: continue
        valid = sf["log_CIR"].notna()
        if valid.sum() == 0: continue
        w = si[valid] / si[valid].sum()
        row = {"sample": sample, "n_formulas": len(sf),
               "mean_CIR": np.average(sf.loc[valid, "log_CIR"], weights=w),
               "sd_CIR": sf["log_CIR"].std(),
               "prop_high_CIR": (sf["log_CIR"] > 1.2).mean(),
               "richness": len(sf)}
        for k, v in sf.groupby("compound_class")["log_CIR"].mean().items():
            row[f"CIR_{k}"] = v
        for k, v in sf["compound_class"].value_counts(normalize=True).items():
            row[f"prop_{k}"] = v
        results.append(row)
    return pd.DataFrame(results)

def run_pipeline():
    print("=" * 70)
    print("CONSTRAINED ISOMERIC RICHNESS (CIR) PIPELINE")
    print("=" * 70)
    print()
    print("[1/6] Loading FT-ICR-MS data...")
    df = pd.read_csv(FTICR_PATH)
    print(f"  Raw: {len(df)} rows x {len(df.columns)} cols")
    meta_cols = ["Mass","C","H","O","N","C13","S","P","Na","El_comp","NeutralMass","Error_ppm","Candidates"]
    sample_cols = [c for c in df.columns if c not in meta_cols]
    print(f"  {len(sample_cols)} sample replicates")
    df = df[df["C"] > 0].copy()
    print(f"  {len(df)} assigned formulas")
    for el in ["C","H","O","N","S","P"]:
        df[el] = pd.to_numeric(df[el], errors="coerce").fillna(0).astype(int)
    def make_formula(row):
        parts = []
        for el in ["C","H","O","N","S","P"]:
            if row[el] > 0:
                parts.append(f"{el}{row[el]}" if row[el] > 1 else el)
        return "".join(parts)
    df["MolecularFormula"] = df.apply(make_formula, axis=1)
    print()
    print("[2/6] Calculating indices...")
    indices = df.apply(calculate_indices, axis=1)
    df = pd.concat([df, indices], axis=1)
    df = df[df["compound_class"] != "unassigned"].copy()
    print("  Compound classes:")
    for cc, n in df["compound_class"].value_counts().items():
        print(f"    {cc:25s} {n:5d} ({100*n/len(df):.1f}%)")
    print()
    print("[3/6] Estimating CIR...")
    cir_vals = []
    for _, r in tqdm(df.iterrows(), total=len(df), desc="CIR"):
        cir_vals.append(estimate_CIR(r, get_constraints(r)))
    df["log_CIR"] = cir_vals
    print()
    print("CIR by class:")
    print(df.groupby("compound_class")["log_CIR"].describe().round(3).to_string())
    cir_med = df.groupby("compound_class")["log_CIR"].median()
    if "condensed_aromatic" in cir_med and "CRAM" in cir_med:
        ok = cir_med["condensed_aromatic"] < cir_med["CRAM"]
        tag = "PASS" if ok else "WARN"
        print(f"  Validation: condensed_aromatic < CRAM: {tag}")
    rng = df["log_CIR"].dropna()
    print(f"  Range: {rng.min():.2f} - {rng.max():.2f}")
    nan_count = df["log_CIR"].isna().sum()
    print(f"  NaN CIR: {nan_count}")
    print()
    print("[4/6] Aggregating to samples...")
    df_s = aggregate_to_sample(df, sample_cols)
    print(f"  {len(df_s)} samples")
    print()
    print("[5/6] Merging metadata...")
    df_s["site_id"] = df_s["sample"].str.extract(r"(S\d+)")[0]
    df_s["replicate"] = df_s["sample"].str.extract(r"\.(\d+)$")[0].astype(int)
    meta = pd.read_csv(META_PATH, encoding="latin-1")
    meta = meta.rename(columns={
        "Sample_ID": "site_id",
        "Latitude.decimal.degrees.": "latitude",
        "Longitude.decimal.degrees": "longitude",
        "Water.temperature.at.50.depth.deg.C.": "temperature_C",
        "Name.of.river": "river",
        "Sampling.location.State.Province": "state",
        "Sampling.location.Country": "country"})
    meta["temperature_C"] = pd.to_numeric(meta["temperature_C"], errors="coerce")
    geo = pd.read_csv(GEOCHEM_PATH)
    geo = geo.rename(columns={"Sample_ID": "site_id", "X681_NPOC_mg_per_L_as_C": "DOC_mgL"})
    geo["DOC_mgL"] = pd.to_numeric(geo["DOC_mgL"], errors="coerce")
    df_s = df_s.merge(meta[["site_id","latitude","longitude","temperature_C","river","state","country"]],
                      on="site_id", how="left")
    df_s = df_s.merge(geo[["site_id","DOC_mgL"]], on="site_id", how="left")
    print("  Done")
    print()
    print("[6/6] Saving...")
    out_cols = ["MolecularFormula","Mass","C","H","O","N","S","P","El_comp",
                "DBE","HC","OC","NC","AI_mod","NOSC","compound_class","log_CIR"]
    df[out_cols].to_csv(os.path.join(OUT_DIR, "formula_level_CIR.csv"), index=False)
    df_s.to_csv(os.path.join(OUT_DIR, "sample_level_CIR.csv"), index=False)
    with open(os.path.join(OUT_DIR, "CIR_summary_stats.txt"), "w") as f:
        f.write("CIR SUMMARY STATISTICS" + NL + "=" * 60 + NL + NL)
        f.write("FORMULA-LEVEL" + NL + "-" * 60 + NL)
        f.write(df.groupby("compound_class")["log_CIR"].describe().round(3).to_string())
        f.write(NL + NL + "SAMPLE-LEVEL" + NL + "-" * 60 + NL)
        for c in ["mean_CIR","sd_CIR","prop_high_CIR","n_formulas"]:
            if c in df_s.columns:
                f.write(f"{c:20s}: mean={df_s[c].mean():.3f} sd={df_s[c].std():.3f} "
                        f"[{df_s[c].min():.3f}, {df_s[c].max():.3f}]" + NL)
    print("  Saved all data files")
    return df, df_s

def make_plots(df_f, df_s):
    print()
    print("Plotting...")
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    ax = axes[0, 0]
    order = df_f.groupby("compound_class")["log_CIR"].median().sort_values().index
    sns.boxplot(data=df_f, x="compound_class", y="log_CIR", order=order, ax=ax, palette="Set2")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Compound Class")
    ax.set_ylabel("log10(CIR)")
    ax.set_title("A) CIR by Class")
    ax = axes[0, 1]
    v = df_f.dropna(subset=["OC", "HC", "log_CIR"])
    sc = ax.scatter(v["OC"], v["HC"], c=v["log_CIR"], cmap="plasma", alpha=0.3, s=3, rasterized=True)
    plt.colorbar(sc, ax=ax, label="log10(CIR)", shrink=0.8)
    ax.set_xlabel("O/C")
    ax.set_ylabel("H/C")
    ax.set_xlim(0, 1.2)
    ax.set_ylim(0, 2.5)
    ax.set_title("B) Van Krevelen")
    ax = axes[0, 2]
    ax.hist(df_s["mean_CIR"].dropna(), bins=25, edgecolor="black", color="steelblue", alpha=0.8)
    med = df_s["mean_CIR"].median()
    ax.axvline(med, color="red", ls="--", label=f"median={med:.2f}")
    ax.set_xlabel("Mean log10(CIR)")
    ax.set_ylabel("Samples")
    ax.set_title("C) Sample CIR")
    ax.legend()
    panels = [
        (axes[1, 0], "temperature_C", "Temperature (C)", "D"),
        (axes[1, 1], "latitude", "Latitude", "E"),
        (axes[1, 2], "DOC_mgL", "DOC (mg/L)", "F")]
    for ax, var, lab, let in panels:
        m = df_s[[var, "mean_CIR"]].notna().all(axis=1)
        if m.sum() > 2:
            ax.scatter(df_s.loc[m, var], df_s.loc[m, "mean_CIR"], alpha=0.6, s=30)
            sl, ic, r, p, se = stats.linregress(df_s.loc[m, var], df_s.loc[m, "mean_CIR"])
            xl = np.linspace(df_s.loc[m, var].min(), df_s.loc[m, var].max(), 100)
            ax.plot(xl, sl * xl + ic, "r-", label=f"r={r:.3f}, p={p:.3f}")
            ax.legend()
        ax.set_xlabel(lab)
        ax.set_ylabel("Mean log10(CIR)")
        ax.set_title(f"{let}) CIR vs {lab}")
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "CIR_results.pdf"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, "CIR_results.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved CIR_results.pdf and CIR_results.png")

def eco_analyses(df_f, df_s):
    print()
    print("=" * 70)
    print("ECOLOGICAL ANALYSES")
    print("=" * 70)
    res = []
    tests = [("Latitude", "latitude"), ("Temperature", "temperature_C"),
             ("DOC", "DOC_mgL"), ("Richness", "richness")]
    for name, x_var in tests:
        m = df_s[[x_var, "mean_CIR"]].notna().all(axis=1)
        if m.sum() > 5:
            r, p = stats.pearsonr(df_s.loc[m, x_var], df_s.loc[m, "mean_CIR"])
            rho, ps = stats.spearmanr(df_s.loc[m, x_var], df_s.loc[m, "mean_CIR"])
            msg = f"{name} vs CIR: r={r:.3f} (p={p:.4f}), rho={rho:.3f} (p={ps:.4f}), n={m.sum()}"
            print(f"  {msg}")
            res.append(msg)
    if "CIR_CRAM" in df_s.columns:
        print()
        print("  CRAM-class CIR:")
        for var in ["temperature_C", "DOC_mgL", "latitude"]:
            m = df_s[[var, "CIR_CRAM"]].notna().all(axis=1)
            if m.sum() > 5:
                r, p = stats.pearsonr(df_s.loc[m, var], df_s.loc[m, "CIR_CRAM"])
                msg = f"  CIR_CRAM vs {var}: r={r:.3f}, p={p:.4f}, n={m.sum()}"
                print(msg)
                res.append(msg)
    site_st = df_s.groupby("site_id")["mean_CIR"].agg(["mean", "std", "count"])
    site_st = site_st[site_st["count"] > 1]
    if len(site_st) > 0:
        cv = (site_st["std"] / site_st["mean"]).mean()
        msg = f"Replicate CV: {cv:.4f} ({len(site_st)} sites)"
        print(f"  {msg}")
        res.append(msg)
    with open(os.path.join(OUT_DIR, "ecological_analysis_results.txt"), "w") as f:
        f.write("ECOLOGICAL ANALYSIS RESULTS" + NL + "=" * 60 + NL + NL)
        for line in res:
            f.write(line + NL)
    print("  Saved ecological_analysis_results.txt")

if __name__ == "__main__":
    df_f, df_s = run_pipeline()
    make_plots(df_f, df_s)
    eco_analyses(df_f, df_s)
    print()
    print("=" * 70)
    print("PIPELINE COMPLETE")
    print("Outputs: " + OUT_DIR)
    print("=" * 70)
