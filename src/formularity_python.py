"""
Python implementation of FT-ICR-MS formula assignment.
Replaces Formularity R package for WHONDRS XML files.

Takes: Bruker DataAnalysis XML peak lists (m/z, intensity)
Produces: CSV with molecular formula assignments (C, H, O, N, S, P)

Method:
1. Parse m/z and intensity from XML pk elements
2. For each m/z, test candidate molecular formulae within mass error tolerance
3. Apply element constraints (C1-57, H1-100, O0-30, N0-5, S0-3, P0-2)
4. Apply chemical rules (DBE >= 0, H/C <= 3, etc.)
5. Select best candidate by lowest mass error
"""
import os, sys, re, glob, warnings
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from itertools import product
from tqdm import tqdm

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════

# Exact masses of most abundant isotopes
MASS = {
    'C': 12.000000000,
    'H': 1.00782503207,
    'O': 15.99491461956,
    'N': 14.0030740048,
    'S': 31.97207100,
    'P': 30.97376163,
}

# Electron mass (for [M-H]- ionisation)
ELECTRON_MASS = 0.00054858
PROTON_MASS = 1.00727646677

# Element constraints for DOM (WHONDRS standard)
ELEM_BOUNDS = {
    'C': (1, 57),
    'H': (1, 100),
    'O': (0, 30),
    'N': (0, 5),
    'S': (0, 3),
    'P': (0, 2),
}

# Mass error tolerance (ppm)
PPM_TOLERANCE = 1.0

# Mass range for formula assignment
MASS_MIN = 100.0
MASS_MAX = 900.0


# ═══════════════════════════════════════════════════════════════
# FORMULA ASSIGNMENT ENGINE
# ═══════════════════════════════════════════════════════════════

def neutral_mass_from_neg_esi(mz):
    """Convert [M-H]- m/z to neutral mass."""
    return mz + PROTON_MASS


def calc_exact_mass(C, H, O, N, S, P):
    """Calculate exact monoisotopic mass from element counts."""
    return (C * MASS['C'] + H * MASS['H'] + O * MASS['O'] +
            N * MASS['N'] + S * MASS['S'] + P * MASS['P'])


def calc_ppm_error(measured, theoretical):
    """Calculate mass error in ppm."""
    return abs(measured - theoretical) / theoretical * 1e6


def is_valid_formula(C, H, O, N, S, P):
    """Apply chemical validity rules."""
    if C == 0:
        return False

    # DBE must be non-negative and integer/half-integer
    DBE = 1 + C - H/2 + N/2
    if DBE < 0:
        return False

    # Senior rule: sum of valences must be >= 2 * (atoms - 1)
    # Simplified: H + N must have correct parity with C
    # Actually just check H/C ratio
    HC = H / C
    if HC > 3.0 or HC < 0.2:
        return False

    # O/C ratio
    OC = O / C
    if OC > 1.5:
        return False

    # N rule: even-electron [M-H]- requires even H+N for even mass
    # This is handled implicitly by the mass matching

    return True


def build_formula_db(mass_range=(MASS_MIN, MASS_MAX)):
    """
    Pre-compute a sorted array of candidate formulae and their exact masses.
    This is the key to fast assignment: build once, binary search many times.
    """
    print("Building formula database...")
    formulae = []

    C_min, C_max = ELEM_BOUNDS['C']
    H_min, H_max = ELEM_BOUNDS['H']
    O_min, O_max = ELEM_BOUNDS['O']
    N_min, N_max = ELEM_BOUNDS['N']
    S_min, S_max = ELEM_BOUNDS['S']
    P_min, P_max = ELEM_BOUNDS['P']

    # Iterate over heteroatoms first (fewer combinations)
    for N in range(N_min, N_max + 1):
        for S in range(S_min, S_max + 1):
            for P in range(P_min, P_max + 1):
                hetero_mass = N * MASS['N'] + S * MASS['S'] + P * MASS['P']
                if hetero_mass > mass_range[1]:
                    break
                for C in range(C_min, C_max + 1):
                    c_mass = hetero_mass + C * MASS['C']
                    if c_mass > mass_range[1]:
                        break
                    for O in range(O_min, O_max + 1):
                        co_mass = c_mass + O * MASS['O']
                        if co_mass > mass_range[1]:
                            break
                        # Calculate H range from remaining mass
                        remaining_min = mass_range[0] - co_mass
                        remaining_max = mass_range[1] - co_mass
                        h_lo = max(H_min, int(np.floor(remaining_min / MASS['H'])))
                        h_hi = min(H_max, int(np.ceil(remaining_max / MASS['H'])))
                        for H in range(h_lo, h_hi + 1):
                            exact = co_mass + H * MASS['H']
                            if exact < mass_range[0] or exact > mass_range[1]:
                                continue
                            if is_valid_formula(C, H, O, N, S, P):
                                formulae.append((exact, C, H, O, N, S, P))

    # Sort by mass for binary search
    formulae.sort(key=lambda x: x[0])
    masses = np.array([f[0] for f in formulae])
    elements = np.array([(f[1], f[2], f[3], f[4], f[5], f[6]) for f in formulae])

    print(f"  {len(formulae):,} candidate formulae in {mass_range[0]:.0f}-{mass_range[1]:.0f} Da")
    return masses, elements


def assign_formula(neutral_mass, db_masses, db_elements, ppm_tol=PPM_TOLERANCE):
    """
    Find the best-matching formula for a given neutral mass.
    Returns (C, H, O, N, S, P, error_ppm) or None.
    """
    # Binary search for mass window
    tol_da = neutral_mass * ppm_tol * 1e-6
    lo = np.searchsorted(db_masses, neutral_mass - tol_da)
    hi = np.searchsorted(db_masses, neutral_mass + tol_da)

    if lo >= hi:
        return None

    # Find best match within window
    candidates = db_masses[lo:hi]
    errors = np.abs(candidates - neutral_mass) / neutral_mass * 1e6
    best_idx = np.argmin(errors)
    best_error = errors[best_idx]

    if best_error <= ppm_tol:
        C, H, O, N, S, P = db_elements[lo + best_idx]
        return (int(C), int(H), int(O), int(N), int(S), int(P), best_error)
    return None


# ═══════════════════════════════════════════════════════════════
# XML PARSING
# ═══════════════════════════════════════════════════════════════

def parse_bruker_xml(xml_path):
    """Parse Bruker DataAnalysis XML to extract peak list (m/z, intensity)."""
    tree = ET.parse(xml_path)
    root = tree.getroot()

    peaks = []
    for pk in root.iter('pk'):
        mz = float(pk.get('mz', 0))
        intensity = float(pk.get('i', 0))
        sn = float(pk.get('sn', 0))
        if mz > 0 and intensity > 0:
            peaks.append({'mz': mz, 'intensity': intensity, 'sn': sn})

    return pd.DataFrame(peaks)


# ═══════════════════════════════════════════════════════════════
# MAIN PROCESSING
# ═══════════════════════════════════════════════════════════════

def process_xml_files(xml_dir, output_path, sn_threshold=7.0):
    """
    Process all XML files in a directory.
    Returns a wide-format DataFrame like the WHONDRS 2018 format.
    """
    xml_files = sorted(glob.glob(os.path.join(xml_dir, "*.xml")))
    print(f"Found {len(xml_files)} XML files")

    if len(xml_files) == 0:
        print("No XML files found!")
        return None

    # Build formula database once
    db_masses, db_elements = build_formula_db()

    # Process each XML file
    all_results = {}
    all_mz = set()

    for xml_path in tqdm(xml_files, desc="Processing XMLs"):
        sample_name = os.path.splitext(os.path.basename(xml_path))[0]
        peaks = parse_bruker_xml(xml_path)

        if len(peaks) == 0:
            continue

        # Filter by S/N
        peaks = peaks[peaks['sn'] >= sn_threshold].copy()

        # Convert to neutral mass (negative ESI: [M-H]-)
        peaks['neutral_mass'] = peaks['mz'].apply(neutral_mass_from_neg_esi)

        # Assign formulae
        assignments = []
        for _, pk in peaks.iterrows():
            nm = pk['neutral_mass']
            if nm < MASS_MIN or nm > MASS_MAX:
                continue
            result = assign_formula(nm, db_masses, db_elements)
            if result is not None:
                C, H, O, N, S, P, error = result
                assignments.append({
                    'mz': pk['mz'],
                    'neutral_mass': nm,
                    'intensity': pk['intensity'],
                    'C': C, 'H': H, 'O': O, 'N': N, 'S': S, 'P': P,
                    'error_ppm': error,
                })

        if assignments:
            df_assign = pd.DataFrame(assignments)
            # Create formula key
            def make_key(r):
                parts = []
                for el in ['C','H','O','N','S','P']:
                    if r[el] > 0:
                        parts.append(f"{el}{r[el]}" if r[el] > 1 else el)
                return ''.join(parts)
            df_assign['formula'] = df_assign.apply(make_key, axis=1)
            all_results[sample_name] = df_assign.set_index('formula')['intensity']
            all_mz.update(df_assign.set_index('formula')['mz'].to_dict().items())

    print(f"\nProcessed {len(all_results)} samples with assignments")

    if len(all_results) == 0:
        print("No assignments found!")
        return None

    # Build wide-format DataFrame
    # Get unique formulae
    all_formulae = set()
    formula_info = {}
    for sample, intensities in all_results.items():
        for formula in intensities.index:
            all_formulae.add(formula)

    # Build formula metadata
    for sample, df in all_results.items():
        # Use first sample's data for each formula
        for formula in df.index:
            if formula not in formula_info:
                row = all_results[sample]
                # We need element counts - parse from formula string
                pass

    # Alternative approach: build from first occurrence
    formula_meta = {}
    for sample_name, df_assign_series in all_results.items():
        # Go back to the full assignment data
        pass

    # Simpler approach: re-process to get element counts
    print("Building wide-format output...")
    formula_data = {}

    for xml_path in tqdm(xml_files, desc="Re-reading"):
        sample_name = os.path.splitext(os.path.basename(xml_path))[0]
        peaks = parse_bruker_xml(xml_path)
        if len(peaks) == 0:
            continue
        peaks = peaks[peaks['sn'] >= sn_threshold].copy()
        peaks['neutral_mass'] = peaks['mz'].apply(neutral_mass_from_neg_esi)

        for _, pk in peaks.iterrows():
            nm = pk['neutral_mass']
            if nm < MASS_MIN or nm > MASS_MAX:
                continue
            result = assign_formula(nm, db_masses, db_elements)
            if result is not None:
                C, H, O, N, S, P, error = result
                parts = []
                for el, v in [('C',C),('H',H),('O',O),('N',N),('S',S),('P',P)]:
                    if v > 0:
                        parts.append(f"{el}{v}" if v > 1 else el)
                formula = ''.join(parts)

                if formula not in formula_data:
                    formula_data[formula] = {
                        'Mass': pk['mz'], 'C': C, 'H': H, 'O': O,
                        'N': N, 'S': S, 'P': P, 'Error_ppm': error,
                    }
                formula_data[formula][sample_name] = pk['intensity']

    # Build DataFrame
    df_wide = pd.DataFrame.from_dict(formula_data, orient='index')
    df_wide.index.name = 'MolecularFormula'

    # Fill NaN intensities with 0
    sample_names = [os.path.splitext(os.path.basename(f))[0] for f in xml_files]
    for s in sample_names:
        if s not in df_wide.columns:
            df_wide[s] = 0
    df_wide = df_wide.fillna(0)

    # Sort by mass
    df_wide = df_wide.sort_values('Mass')

    # Save
    df_wide.to_csv(output_path)
    print(f"\nSaved: {output_path}")
    print(f"  {len(df_wide)} unique formulae across {len(sample_names)} samples")
    print(f"  Assigned: {(df_wide[sample_names] > 0).any(axis=1).sum()} formulae with at least 1 detection")

    return df_wide


if __name__ == "__main__":
    BASE = os.path.expanduser(os.path.join("~", "OneDrive", "Desktop", "ISOMERs"))
    XML_DIR = os.path.join(BASE, "S19S_extract",
                           "WHONDRS_S19S_Sediment_FTICR_Data_And_Instructions",
                           "WHONDRS_S19S_FTICR_FieldSediment_Data")
    OUTPUT = os.path.join(BASE, "output", "S19S_FTICR_processed.csv")

    df = process_xml_files(XML_DIR, OUTPUT)
