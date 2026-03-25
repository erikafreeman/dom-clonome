# Verification Checklist
## Exact numbers the pipeline should produce

### Step 1: Data loading
- Raw CSV: 34,021 rows x 93 columns
- After filtering C > 0: **8,972 assigned formulae**
- Sample columns: **80** replicates across **27** sites
- Mass range: 116.07 — 891.72 Da
- Carbon range: 4 — 57

### Step 2: Compound class counts (must match exactly)
```
CRAM                       2735  (30.5%)
other                      1802  (20.1%)
lipid_like                 1225  (13.7%)
aromatic                   1136  (12.7%)
carbohydrate_like           951  (10.6%)
N_containing                472   (5.3%)
condensed_aromatic          446   (5.0%)
protein_like                205   (2.3%)
```

### Step 3: CIR by compound class (median log10 CIR)
```
aromatic            0.881
condensed_aromatic  0.896
CRAM                0.988
protein_like        1.086
carbohydrate_like   1.199
N_containing        1.237
other               1.252
lipid_like          1.267
```
**Validation: condensed_aromatic (0.896) < CRAM (0.988) = PASS**
**Range: 0.30 — 1.58**
**NaN count: 0**

### Step 4: Sample-level summary
```
mean_CIR:       mean=1.151  sd=0.011  range=[1.135, 1.183]
sd_CIR:         mean=0.152  sd=0.007  range=[0.135, 0.163]
prop_high_CIR:  mean=0.457  sd=0.019  range=[0.412, 0.497]
n_formulas:     mean=2262   sd=925    range=[434, 3871]
```

### Step 5: Ecological correlations (KEY RESULTS)
```
Latitude    vs mean_CIR:   r = -0.383   p = 0.0004
Temperature vs mean_CIR:   r =  0.274   p = 0.0140
DOC         vs mean_CIR:   r =  0.126   p = 0.2658  (NOT significant)
Richness    vs mean_CIR:   r = -0.397   p = 0.0003

CIR_CRAM vs latitude:      r = -0.743   p < 0.0001  *** STRONGEST RESULT ***
CIR_CRAM vs temperature:   r =  0.593   p < 0.0001
CIR_CRAM vs DOC:           r =  0.397   p = 0.0003

Replicate CV: 0.0032 (0.32%) across 27 sites
```

### Things that would indicate a bug
- [ ] Any NaN CIR values for formulae with C > 0
- [ ] Condensed aromatics having HIGHER CIR than CRAM
- [ ] CIR range outside [0.3, 2.5] (hard bounds in code)
- [ ] Replicate CV > 5% (would indicate instability)
- [ ] Fewer than 8,900 or more than 9,000 assigned formulae
