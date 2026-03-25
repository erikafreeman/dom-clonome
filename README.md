# DOM Clonome: Constrained Isomeric Richness of Dissolved Organic Matter

Estimates the hidden structural diversity (isomeric richness) within molecular formulae detected by FT-ICR-MS in dissolved organic matter.

## Key Result

Applied to **691,229 molecular formulae** across **7 datasets** spanning river, sediment, soil, Arctic, and marine environments:

- **Class ordering conserved across all ecosystems**: aromatic always lowest, lipid-like always highest
- **CRAM latitudinal gradient**: r = -0.79 (p < 0.00001, n = 27 sites) -- strongest among all tested metrics
- **Structural depth and compositional breadth are partially decoupled**: r = -0.40
- **Processing intensity gradient**: Arctic < sediment < river < soil < marine

## Quick Start

```bash
pip install -r requirements.txt
python src/CIR_pipeline.py
```

Runtime: ~2 seconds for the primary WHONDRS 2018 dataset.

## Data

Download from ESS-DIVE (place in `data/raw/`):

| Dataset | DOI | Formulae |
|---|---|---|
| WHONDRS 2018 Surface Water | 10.15485/1484811 | 8,972 |
| WHONDRS S19S 2019 Surface Water | 10.15485/1603775 | 166,262 |
| WHONDRS S19S 2019 Sediment | 10.15485/1729719 | 127,136 |
| WHONDRS WROL 2019 | ESS-DIVE | 150,999 |
| WHONDRS YDE21 Arctic | ESS-DIVE | 86,338 |
| Soil DOM (Xia & Liu 2025) | 10.5061/dryad.pvmcvdnw4 | 121,324 |
| Marine DOM (Moye et al. 2025) | 10.1594/PANGAEA.974714 | 30,198 |

## Structure

```
src/
  CIR_pipeline.py          # Primary pipeline (WHONDRS 2018, 306 lines)
  CIR_universal.py          # Universal pipeline (any FT-ICR-MS data)
  formularity_python.py     # Formula assignment from raw XML peak lists
  robustness_checks.py      # Sensitivity analyses
  fix_stats.py              # Pseudoreplication correction
  cross_ecosystem.py        # 7-dataset comparison
results/
  formula_level_CIR.csv     # 8,972 formulae with CIR
  sample_level_CIR.csv      # 80 samples with metadata
  figures/                  # Publication figures
  tables/                   # Summary statistics
  robustness/               # Sensitivity outputs
  cross_ecosystem/          # Per-dataset results
manuscript/
  CIR_Paper_Final.html      # L&O Letters format
docs/
  CIR_METHODS_DETAIL.md     # Full technical spec
  VERIFICATION_CHECKLIST.md # Reproducibility checklist
```

## Scripts

- `src/CIR_pipeline.py` -- Primary analysis (WHONDRS 2018)
- `src/CIR_universal.py` -- Run CIR on any FT-ICR-MS CSV:
  `python src/CIR_universal.py data.csv --label MyData`
- `src/formularity_python.py` -- Assign formulae from Bruker XML peak lists (replaces Formularity R package)

## Citation

> [Authors]. Constrained isomeric richness reveals a latitudinal gradient in hidden DOM structural complexity. Limnology and Oceanography Letters (in preparation).

## License

MIT
