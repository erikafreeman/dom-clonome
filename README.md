# DOM Clonome: Constrained Isomeric Richness of Dissolved Organic Matter

Estimates the hidden structural diversity (isomeric richness) within molecular formulae detected by FT-ICR-MS in dissolved organic matter (DOM).

## Key Result

Applied to **691,229 molecular formulae** across **7 datasets** spanning river, sediment, soil, Arctic, and marine environments, CIR reveals:

- **Class ordering conserved across all ecosystems**: aromatic classes always lowest, lipid-like always highest
- **CRAM latitudinal gradient** (r = -0.79, p < 10^-5, n = 27 sites) — strongest among all tested metrics
- **Structural depth and compositional breadth are partially decoupled** (r = -0.40)
- **Processing intensity gradient**: Arctic < sediment < river < soil < marine

## Quick Start



Runtime: ~2 seconds for the primary WHONDRS 2018 dataset.

## Data

The pipeline expects WHONDRS data in . Download from ESS-DIVE:

| Dataset | DOI | What |
|---|---|---|
| WHONDRS 2018 Surface Water | [10.15485/1484811](https://data.ess-dive.lbl.gov/view/doi:10.15485/1484811) | Primary dataset (27 sites, pre-processed CSV) |
| WHONDRS S19S 2019 Surface Water | [10.15485/1603775](https://data.ess-dive.lbl.gov/view/doi:10.15485/1603775) | 97 global sites (raw XML, use ) |
| WHONDRS S19S 2019 Sediment | [10.15485/1729719](https://data.ess-dive.lbl.gov/view/doi:10.15485/1729719) | Sediment porewater (raw XML) |
| WHONDRS WROL 2019 | ESS-DIVE | Watershed resilience (raw XML) |
| WHONDRS YDE21 Arctic | ESS-DIVE | Yukon Delta Arctic sediment (raw XML) |
| Soil DOM | [10.5061/dryad.pvmcvdnw4](https://datadryad.org/dataset/doi:10.5061/dryad.pvmcvdnw4) | Paddy soil (Xia & Liu 2025, xlsx) |
| Marine DOM | [10.1594/PANGAEA.974714](https://doi.pangaea.de/10.1594/PANGAEA.974714) | Weddell/North Sea (tab-delimited) |

## Project Structure



## Pipeline Scripts

### 
The primary analysis pipeline for the WHONDRS 2018 dataset. Reads CSVs, computes chemical indices, estimates CIR, aggregates to samples, merges metadata, generates plots and statistics.

### 
Universal CIR pipeline for any FT-ICR-MS dataset. Auto-detects format (CSV with element columns, formula strings, or wide/long format).



### 
Python implementation of molecular formula assignment from Bruker DataAnalysis XML peak lists. Replaces the Formularity R package. Builds a 3.7M-formula database and assigns formulae at 1 ppm mass accuracy.



## Citation

If you use this code or the CIR metric, please cite:

> [Authors]. Constrained isomeric richness reveals a latitudinal gradient in hidden DOM structural complexity. *Limnology and Oceanography Letters* (in preparation).

## References

- Leyva et al. (2019) *Faraday Discussions* 218:431-440 (TIMS calibration)
- Koch & Dittmar (2006) *RCMS* 20:926-932 (AI_mod)
- Kim et al. (2003) *Anal Chem* 75:5336-5344 (van Krevelen classes)
- Hertkorn et al. (2006) *GCA* 70:2990-3010 (CRAM)
- Stegen et al. (2022) *mSystems* 7:e00151-22 (WHONDRS)

## License

MIT
