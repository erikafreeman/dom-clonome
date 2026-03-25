"""Cross-ecosystem CIR comparison: River vs Marine DOM"""
import os, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings("ignore")

BASE = os.path.expanduser(os.path.join("~", "OneDrive", "Desktop", "ISOMERs"))
OUT = os.path.join(BASE, "output")

# Load both datasets
r = pd.read_csv(os.path.join(OUT, "formula_level_CIR.csv"))
m = pd.read_csv(os.path.join(OUT, "marine_DOM", "Marine_DOM_PANGAEA_formula_CIR.csv"))

r["ecosystem"] = "River (WHONDRS 2018)"
m["ecosystem"] = "Marine (PANGAEA 2025)"
combined = pd.concat([r, m], ignore_index=True)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Panel A: CIR by class, split by ecosystem
ax = axes[0]
order = combined.groupby("compound_class")["log_CIR"].median().sort_values().index
# Filter to main classes only
main_classes = ["aromatic", "condensed_aromatic", "CRAM", "protein_like",
                "carbohydrate_like", "lipid_like"]
subset = combined[combined["compound_class"].isin(main_classes)]
order_main = [c for c in order if c in main_classes]
sns.boxplot(data=subset, x="compound_class", y="log_CIR", hue="ecosystem",
            order=order_main, ax=ax, palette=["#4C72B0", "#DD8452"], linewidth=0.7,
            fliersize=1)
ax.set_xticklabels([c.replace("_", "\n") for c in order_main], fontsize=8)
ax.set_xlabel("")
ax.set_ylabel("log$_{10}$(CIR)", fontsize=11)
ax.set_title("A) CIR by compound class", fontsize=12, fontweight="bold")
ax.legend(fontsize=8, loc="upper left")

# Panel B: Van Krevelen side by side
ax = axes[1]
rv = r.dropna(subset=["OC", "HC", "log_CIR"])
mv = m.dropna(subset=["OC", "HC", "log_CIR"])
ax.scatter(rv["OC"], rv["HC"], c="#4C72B0", alpha=0.15, s=2, label="River", rasterized=True)
ax.scatter(mv["OC"], mv["HC"], c="#DD8452", alpha=0.15, s=2, label="Marine", rasterized=True)
ax.set_xlabel("O/C", fontsize=11)
ax.set_ylabel("H/C", fontsize=11)
ax.set_xlim(0, 1.2)
ax.set_ylim(0, 2.5)
ax.set_title("B) Van Krevelen: river vs marine", fontsize=12, fontweight="bold")
ax.legend(fontsize=9, markerscale=5)

# Panel C: Density of CIR by ecosystem
ax = axes[2]
bins = np.linspace(0.3, 1.8, 50)
ax.hist(r["log_CIR"].dropna(), bins=bins, alpha=0.6, color="#4C72B0",
        label=f"River (n={len(r):,})", density=True, edgecolor="white", linewidth=0.3)
ax.hist(m["log_CIR"].dropna(), bins=bins, alpha=0.6, color="#DD8452",
        label=f"Marine (n={len(m):,})", density=True, edgecolor="white", linewidth=0.3)
ax.axvline(r["log_CIR"].median(), color="#4C72B0", ls="--", lw=1.5)
ax.axvline(m["log_CIR"].median(), color="#DD8452", ls="--", lw=1.5)
ax.set_xlabel("log$_{10}$(CIR)", fontsize=11)
ax.set_ylabel("Density", fontsize=11)
ax.set_title("C) CIR distribution: river vs marine", fontsize=12, fontweight="bold")
ax.legend(fontsize=9)

plt.tight_layout()
fig_path = os.path.join(OUT, "CIR_cross_ecosystem.png")
fig.savefig(fig_path, dpi=200, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "CIR_cross_ecosystem.pdf"), dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: {fig_path}")

# Summary stats
print("\nCROSS-ECOSYSTEM SUMMARY")
print("=" * 60)
for eco, df in [("River", r), ("Marine", m)]:
    print(f"\n{eco}: {len(df):,} formulae")
    meds = df.groupby("compound_class")["log_CIR"].median().sort_values()
    for cc, med in meds.items():
        print(f"  {cc:25s} {med:.3f}  (~{10**med:.0f} isomers)")
    print(f"  {'Overall':25s} {df['log_CIR'].median():.3f}")

# CRAM comparison
rc = r.loc[r["compound_class"]=="CRAM","log_CIR"]
mc = m.loc[m["compound_class"]=="CRAM","log_CIR"]
print(f"\nCRAM isomeric depth:")
print(f"  River:  median={rc.median():.3f} (~{10**rc.median():.0f} isomers), n={len(rc):,}")
print(f"  Marine: median={mc.median():.3f} (~{10**mc.median():.0f} isomers), n={len(mc):,}")
pct = (10**mc.median() / 10**rc.median() - 1) * 100
print(f"  Marine CRAM has {pct:.0f}% more estimated isomers per formula")

# Class proportion comparison
print(f"\nClass proportions:")
for eco, df in [("River", r), ("Marine", m)]:
    props = df["compound_class"].value_counts(normalize=True)
    top3 = props.head(3)
    print(f"  {eco}: " + ", ".join(f"{c} {v*100:.0f}%" for c, v in top3.items()))

# Save summary
with open(os.path.join(OUT, "cross_ecosystem_summary.txt"), "w") as f:
    f.write("CROSS-ECOSYSTEM CIR COMPARISON\n")
    f.write("River DOM (WHONDRS 2018) vs Marine DOM (PANGAEA Moye et al. 2025)\n")
    f.write("=" * 60 + "\n\n")
    for eco, df in [("River", r), ("Marine", m)]:
        f.write(f"{eco}: {len(df):,} formulae\n")
        meds = df.groupby("compound_class")["log_CIR"].median().sort_values()
        for cc, med in meds.items():
            f.write(f"  {cc:25s} {med:.3f}  (~{10**med:.0f} isomers)\n")
        f.write("\n")
    f.write(f"CRAM comparison:\n")
    f.write(f"  River CRAM median:  {rc.median():.3f} (~{10**rc.median():.0f} isomers)\n")
    f.write(f"  Marine CRAM median: {mc.median():.3f} (~{10**mc.median():.0f} isomers)\n")
    f.write(f"  Marine CRAM has {pct:.0f}% more isomers per formula\n")

print(f"\nSaved cross_ecosystem_summary.txt")
