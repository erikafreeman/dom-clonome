# CIR Analysis — Handoff Package for Independent Verification

## What this is

This package contains everything to reproduce the Constrained Isomeric Richness
(CIR) analysis of WHONDRS FT-ICR-MS dissolved organic matter data.

CIR estimates how many structurally distinct isomers hide behind each molecular
formula detected by mass spectrometry. It is a heuristic index, not a direct
measurement. Its value is that it captures ecological information (especially in
the CRAM compound class) that formula-level metrics miss entirely.

---

## Files

| File | Purpose |
|------|---------|
| `CIR_pipeline.py` | Complete analysis pipeline (306 lines). **Run this.** |
| `HANDOFF_README.md` | This file |
| `VERIFICATION_CHECKLIST.md` | Exact numbers to check against |
| `CIR_METHODS_DETAIL.md` | Every equation and constraint rule, line by line |
| `formula_level_CIR.csv` | Output: 8,972 formulae with CIR scores |
| `sample_level_CIR.csv` | Output: 80 samples with aggregated CIR + metadata |
| `CIR_results.png` | Output: 6-panel diagnostic figure |
| `CIR_summary_stats.txt` | Output: summary statistics |
| `ecological_analysis_results.txt` | Output: correlation results |
| `CIR_Paper_Final.html` | Manuscript (L&O Letters format) |

## Data source

WHONDRS 2018 Surface Water campaign, ESS-DIVE repository.
Three CSVs needed (see VERIFICATION_CHECKLIST.md for exact paths).

## How to run

```bash
pip install pandas numpy matplotlib seaborn scipy tqdm
python CIR_pipeline.py
```

Runtime: ~2 seconds. All outputs appear in `output/`.
Edit lines 13-17 of CIR_pipeline.py if your directory paths differ.
