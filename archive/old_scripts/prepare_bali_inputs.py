"""
prepare_bali_inputs.py

Reads an Excel file of rock inputs, normalizes headers and computes a few derived
columns used by the Bali model (absolute log10(fO2), NaCl molality and its log10).

Dependencies: pandas, numpy, math

Usage (CLI):
    python prepare_bali_inputs.py --in_xlsx /path/to/input.xlsx --sheet 0

If no args provided the script defaults to
    /Users/polsuka/Desktop/Bali et al 2012 model/input.xlsx (first sheet)

Outputs:
    - prints a preview (first 20 rows) of the normalized table
    - writes inputs_normalized.csv in the current working directory

"""

import argparse
import math
from typing import Tuple

import numpy as np
import pandas as pd


# -------- constants --------
DM = {"Mo": 0.025, "Ce": 0.772, "W": 0.0024, "U": 0.0047, "Th": 0.0137, "Nb": 0.21, "La": 0.234}  # ppm
M_MO = 95.95
M_W = 183.84
M_NACL = 58.44
CE_D_CPX, CE_D_GRT, CE_D_RUT = 2.0, 0.4, 2.0
TH_D_CPX, TH_D_GRT, TH_D_RUT = 1.0, 1.19, 0.10
NB_D_CPX, NB_D_GRT, NB_D_RUT = 0.172, 0.204, 200.0
LA_D_CPX, LA_D_GRT, LA_D_RUT = 1.429, 0.204, 1.250
W_D_RUT = 1.250


# Helper: safe log10
def log10_safe(x: float) -> float:
    """Return math.log10(x) for x>0 else NaN."""
    try:
        if x is None:
            return float("nan")
        x = float(x)
        if x > 0.0:
            return math.log10(x)
    except Exception:
        pass
    return float("nan")


def load_excel(path: str, sheet: int | str = 0) -> pd.DataFrame:
    """Read an Excel file and return a raw DataFrame."""
    return pd.read_excel(path, sheet_name=sheet)


def _normalize_colname(col: str) -> str:
    return str(col).strip().lower().replace('\u00b0', '').replace('°', '').replace(' ', '').replace('_', '')


def normalize_headers(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize headers to the internal naming scheme and coerce numeric types.

    Returns a new DataFrame (copy) with renamed columns.
    """
    df = df.copy()
    mapping = {}

    # canonical name keys -> internal name
    key_map = {
        'input': 'rock_id',
        'rock': 'rock_id',
        'p(bar)': 'P_bar', 'pbar': 'P_bar', 'pressure(bar)': 'P_bar', 'pressure': 'P_bar', 'p': 'P_bar',
        'temperature(k)': 'T_K', 'temperature(k)': 'T_K', 'temperaturek': 'T_K', 'tk': 'T_K', 't(k)': 'T_K',
        'temperature(c)': 'T_C', 'temperaturec': 'T_C', 'tc': 'T_C', 'temperature': 'T_K',
        'dfmq(logunits)': 'dFMQ', 'dfmq': 'dFMQ', 'd_fmq': 'dFMQ', 'Δfmq': 'dFMQ',
        'grt': 'mode_grt', 'garnet': 'mode_grt',
        'cpx': 'mode_cpx', 'clinopyroxene': 'mode_cpx',
        'rut': 'mode_rut', 'rutile': 'mode_rut',
        'mo': 'C0_Mo_ppm', 'molibdenum': 'C0_Mo_ppm',
        'ce': 'C0_Ce_ppm', 'cerium': 'C0_Ce_ppm',
        'w': 'C0_W_ppm', 'tungsten': 'C0_W_ppm',
        'u': 'C0_U_ppm', 'uranium': 'C0_U_ppm',
        'th': 'C0_Th_ppm', 'thorium': 'C0_Th_ppm',
        'nb': 'C0_Nb_ppm', 'niobium': 'C0_Nb_ppm',
        'la': 'C0_La_ppm', 'lanthanum': 'C0_La_ppm',
        'nacl_m': 'NaCl_m', 'nacl(mol/kg)': 'NaCl_m', 'naclmolality': 'NaCl_m',
        'nacl_wt_pct': 'NaCl_wt_pct', 'naclwt%': 'NaCl_wt_pct', 'naclwtpct': 'NaCl_wt_pct', 'naclwt': 'NaCl_wt_pct'
    }

    for col in list(df.columns):
        n = _normalize_colname(col)
        if n in key_map:
            mapping[col] = key_map[n]
        else:
            # if already matches internal name (case-insensitive)
            for internal in ['rock_id','P_bar','T_K','dFMQ','mode_grt','mode_cpx','mode_rut',
                             'C0_Mo_ppm','C0_Ce_ppm','C0_W_ppm','C0_U_ppm','C0_Th_ppm','C0_Nb_ppm','C0_La_ppm',
                             'NaCl_m','NaCl_wt_pct','T_C']:
                if n == internal.lower():
                    mapping[col] = internal
                    break

    df = df.rename(columns=mapping)

    # Coerce numeric columns that we will use
    numeric_cols = [
        'P_bar', 'T_K', 'T_C', 'dFMQ',
        'mode_grt','mode_cpx','mode_rut',
        'C0_Mo_ppm','C0_Ce_ppm','C0_W_ppm','C0_U_ppm','C0_Th_ppm','C0_Nb_ppm','C0_La_ppm',
        'NaCl_m','NaCl_wt_pct'
    ]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')

    # If T_K not present but T_C is present, compute T_K
    if 'T_K' not in df.columns and 'T_C' in df.columns:
        df['T_K'] = df['T_C'] + 273.15

    # Ensure rock_id exists
    if 'rock_id' not in df.columns:
        df.insert(0, 'rock_id', [f'rock_{i+1}' for i in range(len(df))])

    return df


def compute_derived(df: pd.DataFrame) -> pd.DataFrame:
    """Compute P_kbar, FMQ base, absolute log10(fO2), NaCl_m and log10_NaCl.

    The function returns a new DataFrame with added columns.
    """
    df = df.copy()

    # Required columns check
    required = ['rock_id', 'P_bar', 'T_K', 'dFMQ']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns after normalization: {', '.join(missing)}")

    # Compute P_kbar
    df['P_kbar'] = df['P_bar'] / 1000.0

    # FMQ base term: base = -25096.3 / T_K + 8.735 + 0.11 * (P_kbar - 1.0) / T_K
    df['log10_fO2_FMQ'] = df['T_K'].apply(lambda t: -25096.3 / t + 8.735 if not pd.isna(t) and t != 0 else float('nan'))
    # add the small pressure term
    df['log10_fO2_FMQ'] = df.apply(lambda r: r['log10_fO2_FMQ'] + (0.11 * (r['P_kbar'] - 1.0) / r['T_K'] if (not pd.isna(r['P_kbar']) and not pd.isna(r['T_K']) and r['T_K'] != 0) else 0.0), axis=1)

    # absolute log10(fO2)
    df['log10_fO2_abs'] = df.apply(lambda r: (r['log10_fO2_FMQ'] + r['dFMQ']) if (not pd.isna(r['log10_fO2_FMQ']) and not pd.isna(r['dFMQ'])) else float('nan'), axis=1)

    # NaCl handling: if NaCl_m present and >0 take log10, else if NaCl_wt_pct present convert
    if 'NaCl_m' in df.columns:
        df['log10_NaCl'] = df['NaCl_m'].apply(lambda x: log10_safe(x) if not pd.isna(x) and x > 0 else float('nan'))
    else:
        df['NaCl_m'] = np.nan
        df['log10_NaCl'] = np.nan

    if 'NaCl_wt_pct' in df.columns and df['NaCl_wt_pct'].notna().any():
        # compute NaCl_m from wt% where wt% > 0
        def wt_to_molality(wt):
            try:
                wt = float(wt)
            except Exception:
                return float('nan')
            if wt <= 0.0 or np.isnan(wt):
                return float('nan')
            return 1000.0 * wt / (M_NACL * (100.0 - wt))

        # fill NaCl_m only where missing
        df['NaCl_m'] = df.apply(lambda r: r['NaCl_m'] if (not pd.isna(r.get('NaCl_m', np.nan)) and r.get('NaCl_m', np.nan) > 0) else wt_to_molality(r['NaCl_wt_pct']) if (not pd.isna(r['NaCl_wt_pct'])) else np.nan, axis=1)
        df['log10_NaCl'] = df['NaCl_m'].apply(lambda x: log10_safe(x) if not pd.isna(x) and x > 0 else float('nan'))

    # Coerce numeric outputs
    out_numeric = ['P_kbar','log10_fO2_FMQ','log10_fO2_abs','NaCl_m','log10_NaCl']
    for c in out_numeric:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    return df


def preview(df: pd.DataFrame, max_rows: int = 20) -> None:
    """Print a clean preview of the DataFrame (first max_rows rows) without index."""
    cols = [
        'rock_id','P_bar','P_kbar','T_K','dFMQ','log10_fO2_FMQ','log10_fO2_abs',
        'mode_grt','mode_cpx','mode_rut',
        'C0_Mo_ppm','C0_Ce_ppm','C0_W_ppm','C0_U_ppm','C0_Th_ppm','C0_Nb_ppm','C0_La_ppm',
        'NaCl_m','log10_NaCl'
    ]
    # Ensure columns exist in df
    cols_present = [c for c in cols if c in df.columns]
    print(df[cols_present].head(max_rows).to_string(index=False))


def main(argv: list | None = None) -> None:
    ap = argparse.ArgumentParser(description='Prepare and normalize Bali inputs')
    ap.add_argument('--in_xlsx', '-i', default='/Users/polsuka/Desktop/Bali et al 2012 model/input.xlsx')
    ap.add_argument('--sheet', '-s', default=0)
    args = ap.parse_args(argv)

    df_raw = load_excel(args.in_xlsx, sheet=args.sheet)
    df_norm = normalize_headers(df_raw)

    try:
        df_calc = compute_derived(df_norm)
    except ValueError as e:
        # clear error for missing required columns
        raise

    # write CSV with the required output columns in the requested order
    out_cols = [
        'rock_id','P_bar','P_kbar','T_K','dFMQ','log10_fO2_FMQ','log10_fO2_abs',
        'mode_grt','mode_cpx','mode_rut',
        'C0_Mo_ppm','C0_Ce_ppm','C0_W_ppm','C0_U_ppm','C0_Th_ppm','C0_Nb_ppm','C0_La_ppm',
        'NaCl_m','log10_NaCl'
    ]
    out_cols_present = [c for c in out_cols if c in df_calc.columns]
    df_calc[out_cols_present].to_csv('inputs_normalized.csv', index=False)

    print(f"Wrote inputs_normalized.csv with columns: {', '.join(out_cols_present)}\n")
    print("Preview (first 20 rows):")
    preview(df_calc, max_rows=20)


if __name__ == '__main__':
    main()
