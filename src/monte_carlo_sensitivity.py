"""
Proper Monte Carlo sensitivity analysis.
Recomputes CIR from scratch for every formula at each iteration.
"""
import os, warnings
import pandas as pd
import numpy as np
from scipy import stats
warnings.filterwarnings("ignore")

BASE = os.path.join(os.path.expanduser('~'), 'OneDrive', 'Desktop', 'ISOMERs', 'output')
df_f = pd.read_csv(os.path.join(BASE, 'formula_level_CIR.csv'))
df_s = pd.read_csv(os.path.join(BASE, 'sample_level_CIR.csv'))

# Precompute arrays
C = df_f['C'].values.astype(float)
H = df_f['H'].values.astype(float)
O = df_f['O'].values.astype(float)
N = df_f['N'].values.astype(float)
S = df_f['S'].values.astype(float)
P = df_f['P'].values.astype(float)
DBE = df_f['DBE'].values
AI = df_f['AI_mod'].values
cc = df_f['compound_class'].values
HC = df_f['HC'].values
OC = df_f['OC'].values

# Precompute constraint counts per formula
n_forb = np.full(len(df_f), 4.0)
n_req = np.zeros(len(df_f))
max_rings = np.full(len(df_f), np.nan)
allow_arom = np.zeros(len(df_f), dtype=bool)

m = cc == 'condensed_aromatic'
allow_arom[m] = True; n_req[m] = 1; max_rings[m] = np.floor(DBE[m] * 0.8)
m = cc == 'aromatic'
allow_arom[m] = True; n_req[m] = 1; n_req[m & (OC > 0.3)] += 1; max_rings[m] = np.maximum(1, np.floor(DBE[m] * 0.6))
m = cc == 'CRAM'
n_req[m & (O >= 2)] = 1; max_rings[m] = np.minimum(3, np.floor(DBE[m] * 0.7))
m = cc == 'carbohydrate_like'
n_req[m] = 1; n_req[m & (DBE >= 1)] += 1; n_forb[m] += 1; max_rings[m] = 2
m = cc == 'lipid_like'
n_req[m & (OC < 0.2)] = 1; max_rings[m] = np.maximum(0, np.floor(DBE[m] * 0.5))
m = cc == 'protein_like'
n_req[m & (N > 0)] = 1; max_rings[m] = 1
m = cc == 'N_containing'
n_req[m] = 1

# Load raw FTICR for sample-level recomputation
DATA_DIR = os.path.join(os.path.expanduser('~'), 'OneDrive', 'Desktop', 'ISOMERs',
                        'WHONDRS_SWMB_DataPackage_2018-11_Rev1', 'data')
df_raw = pd.read_csv(os.path.join(DATA_DIR, 'WHONDRS_SW_Kits_FTICR', 'WHONDRS_SW_Kits_FTICR_Aug2018.csv'))
raw_meta = ['Mass','C','H','O','N','C13','S','P','Na','El_comp','NeutralMass','Error_ppm','Candidates']
sample_cols = [c for c in df_raw.columns if c not in raw_meta]
df_raw = df_raw[df_raw['C'] > 0].copy().reset_index(drop=True)

# Build formula lookup
def make_f(row):
    parts = []
    for el in ['C','H','O','N','S','P']:
        v = int(row.get(el, 0)) if pd.notna(row.get(el, 0)) else 0
        if v > 0: parts.append(f"{el}{v}" if v > 1 else el)
    return ''.join(parts)
df_raw['MF'] = df_raw.apply(make_f, axis=1)

# Map formula to df_f index
mf_to_idx = {}
for i, mf in enumerate(df_f['MolecularFormula']):
    mf_to_idx[mf] = i

# Precompute sample structure: for each sample, which df_f indices are detected
sample_structure = {}
import re
site_map = {}
for col in sample_cols:
    m_site = re.match(r"(S\d+)\.\d+", col)
    if not m_site: continue
    site = m_site.group(1)
    site_map[col] = site
    detected = df_raw[col] > 0
    det_mfs = df_raw.loc[detected, 'MF']
    det_ints = df_raw.loc[detected, col].values
    # Map to df_f indices
    f_idx = []
    f_int = []
    for mf, intensity in zip(det_mfs, det_ints):
        if mf in mf_to_idx:
            f_idx.append(mf_to_idx[mf])
            f_int.append(intensity)
    if f_idx:
        sample_structure[col] = (np.array(f_idx), np.array(f_int, dtype=float))

# Get site latitudes
site_lat = df_s.groupby('site_id')['latitude'].first().to_dict()

def compute_cir_full(w_c, w_dbe, w_o, w_o_carb, w_n, w_req, w_forb, w_arom, w_ring):
    """Fully vectorised CIR for all formulae."""
    carbon = np.log10(np.maximum(1, C - 3)) * w_c
    dbe_c = np.maximum(0.1, w_dbe * (DBE - 4)**2 + 1.2)
    oxy = np.where(cc == 'carbohydrate_like',
                   np.log10(np.maximum(1, O)) * w_o_carb,
                   np.log10(np.maximum(1, O)) * w_o)
    nit = np.log10(np.maximum(1, N + 1)) * w_n
    rp = n_req * w_req
    fp = n_forb * w_forb
    ap = np.where(allow_arom & (AI >= 0.5), AI * w_arom, 0)
    rgp = np.where(~np.isnan(max_rings), np.maximum(0, (DBE - max_rings)) * w_ring, 0)
    raw = carbon + dbe_c + oxy + nit - rp - fp - ap - rgp
    return np.clip(0.78 + (raw / 3.0) * 0.82, 0.3, 1.70)

# Monte Carlo
np.random.seed(42)
N_ITER = 500
defaults = [0.8, -0.15, 0.7, 0.4, 0.5, 0.12, 0.05, 0.8, 0.1]
names = ['w_carbon', 'w_dbe', 'w_oxy', 'w_oxy_carb', 'w_nitrogen', 'w_req', 'w_forb', 'w_arom', 'w_ring']

mc_results = []
for i in range(N_ITER):
    # Perturb each coefficient +-10%
    perturbed = [v * np.random.uniform(0.9, 1.1) for v in defaults]
    cir = compute_cir_full(*perturbed)

    # Aggregate to site-level CRAM means
    site_cram_cir = {}
    site_counts = {}
    for col, (fidx, fint) in sample_structure.items():
        site = site_map[col]
        # Get CRAM formulae
        cram_mask = cc[fidx] == 'CRAM'
        if cram_mask.sum() == 0: continue
        cram_cir = cir[fidx[cram_mask]]
        cram_int = fint[cram_mask]
        w = cram_int / cram_int.sum()
        wcir = np.average(cram_cir, weights=w)
        if site not in site_cram_cir:
            site_cram_cir[site] = []
        site_cram_cir[site].append(wcir)

    # Average triplicates
    lats, cirs = [], []
    for site, vals in site_cram_cir.items():
        if site in site_lat and not np.isnan(site_lat[site]):
            lats.append(site_lat[site])
            cirs.append(np.mean(vals))

    if len(lats) > 5:
        r, p = stats.pearsonr(lats, cirs)
        mc_results.append({'iter': i, 'r': r, 'p': p})

mc_df = pd.DataFrame(mc_results)

print('MONTE CARLO SENSITIVITY ANALYSIS (FULL RECOMPUTATION)')
print('='*60)
print(f'Iterations: {len(mc_df)}')
print(f'Each iteration: perturb all 9 coefficients by U[-10%, +10%]')
print(f'Recompute CIR for all 8,972 formulae, reaggregate to sites')
print()
print(f'CRAM-latitude Pearson r:')
print(f'  Original:      r = -0.785')
print(f'  MC mean:       r = {mc_df["r"].mean():.3f}')
print(f'  MC sd:         {mc_df["r"].std():.4f}')
print(f'  MC 2.5%:       r = {mc_df["r"].quantile(0.025):.3f}')
print(f'  MC 97.5%:      r = {mc_df["r"].quantile(0.975):.3f}')
print(f'  MC range:      [{mc_df["r"].min():.3f}, {mc_df["r"].max():.3f}]')
print(f'  % |r| > 0.7:   {(mc_df["r"].abs() > 0.7).mean()*100:.1f}%')
print(f'  % |r| > 0.6:   {(mc_df["r"].abs() > 0.6).mean()*100:.1f}%')
print(f'  % p < 0.001:   {(mc_df["p"] < 0.001).mean()*100:.1f}%')
print(f'  % p < 0.01:    {(mc_df["p"] < 0.01).mean()*100:.1f}%')

# Also check richness correlation
print()
print('CIR-richness Pearson r under perturbation:')
rich_results = []
for i in range(min(100, N_ITER)):
    perturbed = [v * np.random.uniform(0.9, 1.1) for v in defaults]
    cir = compute_cir_full(*perturbed)

    site_mean_cir = {}
    site_richness = {}
    for col, (fidx, fint) in sample_structure.items():
        site = site_map[col]
        valid = ~np.isnan(cir[fidx])
        if valid.sum() == 0: continue
        w = fint[valid] / fint[valid].sum()
        mcir = np.average(cir[fidx[valid]], weights=w)
        if site not in site_mean_cir:
            site_mean_cir[site] = []
            site_richness[site] = []
        site_mean_cir[site].append(mcir)
        site_richness[site].append(len(fidx))

    cirs_s, richs_s = [], []
    for site in site_mean_cir:
        cirs_s.append(np.mean(site_mean_cir[site]))
        richs_s.append(np.mean(site_richness[site]))
    if len(cirs_s) > 5:
        r, p = stats.pearsonr(richs_s, cirs_s)
        rich_results.append(r)

rich_arr = np.array(rich_results)
print(f'  MC mean r: {rich_arr.mean():.3f}, range: [{rich_arr.min():.3f}, {rich_arr.max():.3f}]')

# Save
out = os.path.join(BASE, '..', 'dom-clonome', 'results', 'robustness', 'monte_carlo_sensitivity.txt')
with open(out, 'w') as f:
    f.write('Monte Carlo Sensitivity Analysis\n')
    f.write(f'Iterations: {len(mc_df)}, Perturbation: all coefficients +/-10%\n')
    f.write(f'CRAM-lat r: mean={mc_df["r"].mean():.3f}, sd={mc_df["r"].std():.4f}, ')
    f.write(f'95% CI=[{mc_df["r"].quantile(0.025):.3f}, {mc_df["r"].quantile(0.975):.3f}]\n')
    f.write(f'% |r|>0.7: {(mc_df["r"].abs() > 0.7).mean()*100:.1f}%\n')
    f.write(f'% p<0.001: {(mc_df["p"] < 0.001).mean()*100:.1f}%\n')
print(f'\nSaved: {out}')
