# CIR Methods — Complete Technical Specification
## Every equation, threshold, and constraint rule

---

## 1. Formula parsing

Element columns (C, H, O, N, S, P) are read directly from the WHONDRS CSV.
Rows with C = 0 are excluded (unassigned mass features).
All element counts cast to integer after coercing non-numeric values to 0.

## 2. Chemical index calculation

### DBE (Double Bond Equivalence)
```
DBE = 1 + C - H/2 + N/2
```

### AI_mod (Modified Aromaticity Index)
Koch & Dittmar (2006) Rapid Commun. Mass Spectrom. 20:926-932
```
numerator   = 1 + C - O/2 - S - H/2
denominator = C - O/2 - S - N - P
AI_mod = max(0, numerator / denominator)   if denominator > 0
AI_mod = 0                                  if denominator <= 0
```

### NOSC (Nominal Oxidation State of Carbon)
LaRowe & Van Cappellen (2011) Geochim. Cosmochim. Acta 75:2030-2042
```
NOSC = 4 - (4C + H - 3N - 2O + 5P - 2S) / C
```

### Elemental ratios
```
H/C = H / C
O/C = O / C
N/C = N / C
```

## 3. Compound class assignment

Following Kim et al. (2003) Anal. Chem. 75:5336-5344, modified by
Sleighter & Hatcher (2007) J. Mass Spectrom. 42:559-574.

Applied in this priority order (first match wins):
```
1. AI_mod >= 0.66                        → condensed_aromatic
2. AI_mod >= 0.50                        → aromatic
3. H/C < 1.5 AND O/C < 0.5              → CRAM
4. H/C >= 1.5 AND O/C >= 0.5            → carbohydrate_like
5. H/C >= 1.5 AND O/C < 0.5             → lipid_like
6. N > 0 AND H/C > 1.3                  → protein_like
7. N > 0                                → N_containing
8. everything else                       → other
```

## 4. Constraint rules per compound class

### Universal constraints (all classes)
```
Forbidden: peroxide (O-O bonds), triple bonds, 3-membered rings, 4-membered rings
Required:  (none)
Max rings: (none)
Allow aromatic: False
```

### condensed_aromatic
```
Allow aromatic: True
Required: aromatic_ring
Max rings: int(DBE * 0.8)
```

### aromatic
```
Allow aromatic: True
Required: aromatic_ring
Required (if O/C > 0.3): hydroxyl_or_methoxy
Max rings: max(1, int(DBE * 0.6))
```

### CRAM
```
Required (if O >= 2): carboxyl_group
Max rings: min(3, int(DBE * 0.7))
```

### carbohydrate_like
```
Required: hydroxyl_group
Required (if DBE >= 1): ring_oxygen
Max rings: 2
Additional forbidden: double_bond_CC
```

### lipid_like
```
Required (if O/C < 0.2): carboxyl_or_ester
Max rings: max(0, int(DBE * 0.5))
```

### protein_like
```
Required (if N > 0): amine_or_amide
Max rings: 1
```

### N_containing
```
Required: nitrogen_heterocycle_or_amine
```

## 5. CIR scoring function

### Components (all additive before penalties)
```
carbon_complexity     = log10(max(1, C - 3)) * 0.8
dbe_contribution      = max(0.1, -0.15 * (DBE - 4)^2 + 1.2)
oxygen_contribution   = log10(max(1, O)) * 0.7          [general]
                      = log10(max(1, O)) * 0.4          [carbohydrate_like only]
nitrogen_contribution = log10(max(1, N + 1)) * 0.5
```

### Penalties (subtracted)
```
required_penalty    = count(required substructures) * 0.12
forbidden_penalty   = count(forbidden substructures) * 0.05
aromaticity_penalty = AI_mod * 0.8    [only if allow_aromatic=True AND AI_mod >= 0.5]
                    = 0               [otherwise]
ring_penalty        = max(0, DBE - max_rings) * 0.1    [if max_rings is set]
                    = 0                                  [if max_rings is None]
```

### Assembly and calibration
```
log_CIR_raw = carbon_complexity + dbe_contribution + oxygen_contribution
              + nitrogen_contribution - required_penalty - forbidden_penalty
              - aromaticity_penalty - ring_penalty

log_CIR_calibrated = 0.78 + (log_CIR_raw / 3.0) * 0.82

log_CIR = clip(log_CIR_calibrated, 0.3, 2.5)
```

Calibration maps to TIMS-observed range:
- log10(6)  = 0.778  ≈ 0.78  (minimum observed isomers per formula)
- log10(40) = 1.602  ≈ 1.60  (maximum observed isomers per formula)
- Source: Leyva et al. (2019) Faraday Discuss. 218:431-440

### Worked example: C15H20O7 (a CRAM formula)
```
C=15, H=20, O=7, N=0, S=0, P=0
DBE = 1 + 15 - 10 + 0 = 6
H/C = 1.333, O/C = 0.467 → class = CRAM (H/C<1.5, O/C<0.5)
AI_mod = (1 + 15 - 3.5 - 0 - 10) / (15 - 3.5 - 0 - 0 - 0) = 2.5/11.5 = 0.217

Constraints for CRAM:
  forbidden = [peroxide, triple_bond, ring_3, ring_4] → 4 items
  required = [carboxyl_group] (O>=2) → 1 item
  max_rings = min(3, int(6*0.7)) = min(3, 4) = 3

Scoring:
  carbon_complexity     = log10(max(1, 12)) * 0.8 = 1.079 * 0.8 = 0.863
  dbe_contribution      = max(0.1, -0.15*(6-4)^2 + 1.2) = max(0.1, 0.6) = 0.600
  oxygen_contribution   = log10(max(1, 7)) * 0.7 = 0.845 * 0.7 = 0.592
  nitrogen_contribution = log10(max(1, 1)) * 0.5 = 0.0
  required_penalty      = 1 * 0.12 = 0.120
  forbidden_penalty     = 4 * 0.05 = 0.200
  aromaticity_penalty   = 0 (allow_aromatic is False)
  ring_penalty          = max(0, 6 - 3) * 0.1 = 0.300

  log_CIR_raw = 0.863 + 0.600 + 0.592 + 0.0 - 0.120 - 0.200 - 0.0 - 0.300 = 1.435
  log_CIR_calibrated = 0.78 + (1.435 / 3.0) * 0.82 = 0.78 + 0.392 = 1.172
  log_CIR = clip(1.172, 0.3, 2.5) = 1.172

  Estimated isomers = 10^1.172 ≈ 14.9
  (Within TIMS range of 6-40: PASS)
```

## 6. Sample-level aggregation

For each sample column (80 total):
1. Select formulae with intensity > 0 in that sample
2. Compute weights = intensity / sum(intensity) for formulae with valid CIR
3. mean_CIR = weighted average of log_CIR
4. sd_CIR = unweighted standard deviation of log_CIR
5. prop_high_CIR = fraction of formulae with log_CIR > 1.2
6. Compute mean CIR per compound class (unweighted)
7. Compute compound class proportions

## 7. Metadata merge

Sample columns like "S000003.1" are parsed:
- site_id = "S000003" (regex: S\d+)
- replicate = 1 (regex: \.(\d+)$)

Merged with metadata on site_id to get latitude, longitude, temperature, river, state.
Merged with geochemistry on site_id to get DOC (NPOC).

## 8. Statistical tests

- Pearson r and p-value for all continuous correlations
- Spearman rho and p-value for rank correlations
- Coefficient of variation across triplicates per site
