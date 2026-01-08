#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compute_logfO2_abs.py

- Reads an Excel sheet with columns like:
    'input', 'P (bar)', 'Temperature (k)', 'ΔFMQ (log units)' (headers are matched flexibly)
- Shows the table with an Excel-style row number (header=1, first data row=2)
- Asks which Excel line (>=2) to use
- Computes ONLY:
    log10(fO2)_FMQ = (-25096.3 / T_K) + 8.735 + 0.11 * ((P_bar - 1) / T_K)   <-- P in BAR
    log10(fO2)_abs = log10(fO2)_FMQ + ΔFMQ  (if ΔFMQ exists)
"""

import argparse
import math
import pandas as pd
import numpy as np

# -------- constants --------
DM = {"MO": 0.025, "Ce": 0.772, "W": 0.0024, "U": 0.0047, "Th": 0.0137, "Nb": 0.21, "La": 0.234}  # ppm
M_MO = 95.95
M_W = 183.84
M_NACL = 58.44
CE_D_CPX, CE_D_GRT, CE_D_RUT = 2.0, 0.4, 2.0
TH_D_CPX, TH_D_GRT, TH_D_RUT = 1.19, 0.613, 0.10
NB_D_CPX, NB_D_GRT, NB_D_RUT = 0.172, 0.204, 200.0
LA_D_CPX, LA_D_GRT, LA_D_RUT = 1.429, 0.204, 1.250
W_D_RUT = 1.250

# Precompute scale log10 constants once so we don't recompute them inside loops
try:
    SCALE_LOG10_M_MO = math.log10(M_MO * 1000.0)
except Exception:
    SCALE_LOG10_M_MO = float('nan')

try:
    SCALE_LOG10_M_W = math.log10(M_W * 1000.0)
except Exception:
    SCALE_LOG10_M_W = float('nan')

# -------- header normalization helpers --------
def _norm(s: str) -> str:
    if not isinstance(s, str): return ""
    x = s.strip().lower()
    for ch in ("°", " ", "\t", "(", ")", "_", "-"):
        x = x.replace(ch, "")
    return x

HEADER_MAP = {
    # identifiers
    "input": "rock_id", "id": "rock_id", "sample": "rock_id",
    # pressure (BAR)
    "pbar": "P_bar", "p(bar)": "P_bar", "pressurebar": "P_bar", "p": "P_bar",
    # temperature
    "temperaturek": "T_K", "temperature(k)": "T_K", "tk": "T_K",
    "temperaturec": "T_C", "temperature(c)": "T_C", "tc": "T_C",
    # ΔFMQ
    "Δfmq(logunits)": "dFMQ", "dfmq(logunits)": "dFMQ", "Δfmq": "dFMQ", "dfmq": "dFMQ",
}

# -------- I/O --------
def load_excel(path: str, sheet=0) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet)

def normalize_headers(df: pd.DataFrame) -> pd.DataFrame:
    # build normalized key map for original columns
    norm_map = {col: _norm(col) for col in df.columns}

    ren = {}
    for col, key in norm_map.items():
        if key in HEADER_MAP:
            ren[col] = HEADER_MAP[key]
    out = df.rename(columns=ren).copy()

    # Fallback: if dFMQ wasn't matched, look for any column whose normalized name contains 'dfmq' or 'fmq'
    if 'dFMQ' not in out.columns:
        for col, key in norm_map.items():
            if 'dfmq' in key or ('fmq' in key and 'd' in key):
                out['dFMQ'] = df[col]
                break
        else:
            # also accept columns that include 'fmq' anywhere
            for col, key in norm_map.items():
                if 'fmq' in key:
                    out['dFMQ'] = df[col]
                    break

    # Ensure rock_id
    if "rock_id" not in out.columns:
        if "input" in out.columns:
            out["rock_id"] = out["input"].astype(str)
        else:
            out.insert(0, "rock_id", [f"rock_{i+1}" for i in range(len(out))])

    # Celsius -> Kelvin if needed
    if "T_K" not in out.columns and "T_C" in out.columns:
        out["T_K"] = pd.to_numeric(out["T_C"], errors="coerce") + 273.15

    # Coerce numerics (handle unicode minus and stray text) - robust for negative dFMQ
    for c in ("P_bar", "T_K", "dFMQ"):
        if c in out.columns:
            # convert to string safely
            s = out[c].astype(str)
            # replace unicode minus (U+2212) with ASCII hyphen-minus and strip whitespace
            s = s.str.replace('\u2212', '-', regex=False).str.strip()
            # allow decimal commas -> convert to dot
            s = s.str.replace(',', '.', regex=False)
            # common textual NA values -> empty string
            s = s.replace({'none': '', 'na': '', 'n/a': '', 'nan': '', '': ''})
            # finally coerce to numeric (invalid -> NaN)
            out[c] = pd.to_numeric(s, errors='coerce')

    # Required columns
    for req in ("rock_id", "P_bar", "T_K"):
        if req not in out.columns:
            raise ValueError(f"Missing required column: {req}")

    return out

    # end normalize_headers

# -------- calculations (ONLY up to log10(fO2)_abs) --------
def compute_base_fmq_log10(T_K: float, P_bar: float) -> float:
    """
    FMQ base (log10) with P input in BAR (do NOT convert to kbar):
      base = (-25096.3 / T_K) + 8.735 + 0.11 * ((P_bar - 1.0) / T_K)
    """
    # use P_bar directly (pressure must be provided in BAR)
    try:
        p_bar = float(P_bar)
    except Exception:
        p_bar = float('nan')
    return (-25096.3 / T_K) + 8.735 + 0.11 * ((p_bar - 1.0) / T_K)

# -------- main --------
def main():
    ap = argparse.ArgumentParser(description="Pick an Excel line (>=2) and compute ONLY log10(fO2)_FMQ and log10(fO2)_abs.")
    ap.add_argument("--in_xlsx", default="/Users/polsuka/Desktop/Bali et al 2012 model/input.xlsx")
    ap.add_argument("--sheet", default=0)
    # Non-interactive options
    ap.add_argument("--row", type=int, help="Excel row number to use (start at 2). Overrides interactive prompt.")
    ap.add_argument("--salinity-index", type=int, choices=range(1, 7), help="Select salinity option by index (1-5) or 6 for ALL.")
    ap.add_argument("--salinity-wt", type=float, help="Directly specify a single salinity wt.% NaCl to use (overrides index).")
    ap.add_argument("--all-salinity", action="store_true", help="Use all predefined salinity options (same as choosing ALL interactively).")
    # Saving of results to disk has been disabled per user request; do not add a --save flag
    args = ap.parse_args()

    df_raw = load_excel(args.in_xlsx, args.sheet)
    df = normalize_headers(df_raw)

    # Show table with Excel-style row numbers (header=1, first data row=2)
    disp = df.copy()
    disp.insert(0, "Excel_row", range(2, 2 + len(disp)))
    # show_cols = [c for c in ("Excel_row","rock_id","P_bar","T_K","dFMQ") if c in disp.columns]
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print("\nTable (Excel_row starts at 2) — showing all detected columns:")
        print(disp.to_string(index=False))

    # Choose Excel row (allow non-interactive override via --row)
    if args.row is not None:
        sel_i = args.row
    else:
        sel = input("\nWhich Excel line do you want to use for absolute log10(fO2)? (start at 2): ")
        try:
            sel_i = int(sel)
        except Exception:
            print(f"Invalid selection: {sel!r} — please enter an integer Excel row number starting at 2.")
            return

    idx = sel_i - 2
    if idx < 0 or idx >= len(df):
        print(f"Selected Excel line {sel_i} is out of range (data rows are 2..{len(df)+1}).")
        return

    row = df.iloc[idx]

    # Extract inputs
    rock_id = row.get("rock_id")
    T_K = row.get("T_K")
    P_bar = row.get("P_bar")

    # helper: find a value in the row by normalized header variants
    def find_val(row, want_names):
        for col in row.index:
            if _norm(col) in want_names:
                return row.get(col)
        return None

    # extract modal proportions and trace elements (accept many header variants)
    modal_grt = find_val(row, {"modegrt", "grt", "garnet"})
    modal_cpx = find_val(row, {"modecpx", "cpx", "clinopyroxene"})
    modal_rut = find_val(row, {"moderut", "rut", "rutile"})

    # trace elements
    MO_v = find_val(row, {"mo", "molibdenum", "c0moppm", "c0_mo_ppm", "c0_mo"})
    Ce_v = find_val(row, {"ce", "cerium", "c0ceppm", "c0_ce_ppm", "c0_ce"})
    W_v = find_val(row, {"w", "tungsten", "c0wppm", "c0_w_ppm", "c0_w"})
    Th_v = find_val(row, {"th", "thorium", "c0thppm", "c0_th_ppm", "c0_th"})
    Nb_v = find_val(row, {"nb", "niobium", "c0nbppm", "c0_nb_ppm", "c0_nb"})
    La_v = find_val(row, {"la", "lanthanum", "c0lappm", "c0_la_ppm", "c0_la"})
    U_v = find_val(row, {"u", "uranium", "c0uppm", "c0_u_ppm", "c0_u"})

    # coerce trace/modals to floats where possible
    def to_float_or_nan(x):
        try:
            if x is None:
                return float('nan')
            if isinstance(x, (float, int, np.number)):
                return float(x)
            s = str(x).strip()
            s = s.replace('\u2212', '-')
            s = s.replace(',', '.')
            if s == '' or s.lower() in ('nan', 'na', 'n/a', 'none'):
                return float('nan')
            return float(s)
        except Exception:
            return float('nan')

    modal_grt = to_float_or_nan(modal_grt)
    modal_cpx = to_float_or_nan(modal_cpx)
    modal_rut = to_float_or_nan(modal_rut)
    MO_v = to_float_or_nan(MO_v)
    Ce_v = to_float_or_nan(Ce_v)
    W_v = to_float_or_nan(W_v)
    U_v = to_float_or_nan(U_v)
    Th_v = to_float_or_nan(Th_v)
    Nb_v = to_float_or_nan(Nb_v)
    La_v = to_float_or_nan(La_v)

    # Ensure dFMQ is always a float (NaN if missing/invalid). Accept negative values.
    if "dFMQ" in row.index:
        raw_dfmq = row.get("dFMQ")
        try:
            dFMQ = float(raw_dfmq)
        except Exception:
            dFMQ = float('nan')
    else:
        dFMQ = float('nan')

    # Basic checks
    if pd.isna(T_K) or pd.isna(P_bar):
        print("Missing numeric T_K or P_bar for selected row — cannot compute FMQ.")
        print(row.to_dict())
        return

    # compute FMQ base using T_K and P_bar provided (P in BAR)
    base = compute_base_fmq_log10(float(T_K), float(P_bar))
    base = float(base)  # ensure Python float (allow negatives)

    # compute absolute log10 fO2 as float: log10_fO2_abs = log10_fO2_FMQ + dFMQ
    if not pd.isna(dFMQ):
        log10_fO2_abs = float(base) + float(dFMQ)
    else:
        log10_fO2_abs = float('nan')

    # (single-row preview removed; final results will be shown after salinity adjustments)

    # Diagnostic output when absolute value is missing/NaN
    if pd.isna(log10_fO2_abs):
        try:
            p_bar_dbg = float(P_bar)
        except Exception:
            p_bar_dbg = float('nan')
        print("\nDEBUG: Missing/invalid absolute log10(fO2). Diagnostic values:")
        print(f"  rock_id: {rock_id!r}")
        print(f"  T_K: {T_K!r} (type: {type(T_K)})")
        print(f"  P_bar: {P_bar!r} (type: {type(P_bar)}) -> P_bar (float): {p_bar_dbg!r}")
        print(f"  FMQ base: {base!r}")
        print(f"  dFMQ: {dFMQ!r} (type: {type(dFMQ)})")

    # Ask user for salinity choice and compute NaCl molality(s)
    sal_options = [0.01, 5.0, 10.0, 15.0, 20.0]
    print("\nSalinity options (wt.% NaCl):")
    for i, v in enumerate(sal_options, start=1):
        print(f"  {i}) {v}%")
    print(f"  {len(sal_options)+1}) ALL of the above")
    # Determine salinity selection (allow non-interactive overrides)
    if args.salinity_wt is not None:
        sel_wts = [args.salinity_wt]
    elif args.salinity_index is not None:
        sel_sal_i = args.salinity_index
        if 1 <= sel_sal_i <= len(sal_options):
            sel_wts = [sal_options[sel_sal_i - 1]]
        else:
            # index '6' or other allowed value -> ALL
            sel_wts = sal_options[:]
    elif args.all_salinity:
        sel_wts = sal_options[:]
    else:
        sel_sal = input(f"Choose salinity option (1-{len(sal_options)+1}): ")
        try:
            sel_sal_i = int(sel_sal)
        except Exception:
            print(f"Invalid salinity selection {sel_sal!r}, defaulting to ALL.")
            sel_sal_i = len(sal_options) + 1

        if sel_sal_i >= 1 and sel_sal_i <= len(sal_options):
            sel_wts = [sal_options[sel_sal_i - 1]]
        else:
            # ALL
            sel_wts = sal_options[:]

    def wt_to_molality(wt):
        try:
            wt = float(wt)
        except Exception:
            return float('nan')
        if wt <= 0.0 or wt >= 100.0 or math.isnan(wt):
            return float('nan')
        return 1000.0 * wt / (M_NACL * (100.0 - wt))

    results = []
    # compute base FMQ once (depends only on T_K and P_bar)
    base = compute_base_fmq_log10(float(T_K), float(P_bar))
    base = float(base)

    for wt in sel_wts:
        NaCl_m = wt_to_molality(wt)
        try:
            log10_NaCl = math.log10(NaCl_m) if (not math.isnan(NaCl_m) and NaCl_m > 0.0) else float('nan')
        except Exception:
            log10_NaCl = float('nan')

        if not pd.isna(dFMQ):
            log10_fO2_abs = float(base) + float(dFMQ)
        else:
            log10_fO2_abs = float('nan')

        # compute Mo-related quantities using chosen log10_fO2_abs, log10_NaCl and T_K
        if not (pd.isna(log10_NaCl) or pd.isna(log10_fO2_abs) or pd.isna(T_K)):
            try:
                # LogMO: linear combination (log10 scale) using the precomputed variables
                LogMO = (0.435 * float(log10_fO2_abs)) + (0.42 * float(log10_NaCl)) - (1.8 * (1000.0 / float(T_K))) + 4.8

                # Compute E_MO using the exact requested formula:
                #   E_MO = 10**LogMO * M_MO * 1000
                # To avoid exponent overflow, check combined log10(E_MO) = LogMO + log10(M_MO*1000)
                # use precomputed SCALE_LOG10_M_MO
                scale_log10 = SCALE_LOG10_M_MO
                combined_log = LogMO + scale_log10
                if (not math.isfinite(LogMO)) or (not math.isfinite(scale_log10)) or (not math.isfinite(combined_log)) or combined_log > 308 or combined_log < -308:
                    E_MO = float('nan')
                    MO_D_CPX = float('nan')
                    MO_D_GRT = float('nan')
                    MO_D_RUT = float('nan')
                else:
                    E_MO = (10.0 ** LogMO) * M_MO * 1000.0
                    if not math.isfinite(E_MO) or E_MO == 0.0:
                        MO_D_CPX = float('nan')
                        MO_D_GRT = float('nan')
                        MO_D_RUT = float('nan')
                    else:
                        MO_D_CPX = 40.0 / E_MO
                        MO_D_GRT = 12.0 / E_MO
                        MO_D_RUT = 87670.0 / E_MO
            except Exception:
                LogMO = float('nan')
                E_MO = float('nan')
                MO_D_CPX = float('nan')
                MO_D_GRT = float('nan')
                MO_D_RUT = float('nan')
        else:
            LogMO = float('nan')
            E_MO = float('nan')
            MO_D_CPX = float('nan')
            MO_D_GRT = float('nan')
            MO_D_RUT = float('nan')

        # --- W (tungsten) calculations ---
        # Moved here so bulk calculations below see W_D_* values
        if not (pd.isna(log10_fO2_abs) or pd.isna(T_K)):
            try:
                # LogW as specified
                LogW = (0.07 * float(log10_fO2_abs)) - (4.7236 * (1000.0 / float(T_K))) + 4.4271

                # Compute E_W = 10**LogW * M_W * 1000 with overflow guarding
                # use precomputed SCALE_LOG10_M_W
                scale_log10_w = SCALE_LOG10_M_W
                combined_log_w = LogW + scale_log10_w
                if (not math.isfinite(LogW)) or (not math.isfinite(scale_log10_w)) or (not math.isfinite(combined_log_w)) or combined_log_w > 308 or combined_log_w < -308:
                    E_W = float('nan')
                    W_D_CPX = float('nan')
                    W_D_GRT = float('nan')
                else:
                    E_W = (10.0 ** LogW) * M_W * 1000.0
                    if not math.isfinite(E_W) or E_W == 0.0:
                        W_D_CPX = float('nan')
                        W_D_GRT = float('nan')
                    else:
                        W_D_CPX = 60.0 / E_W
                        W_D_GRT = 12.0 / E_W
            except Exception:
                LogW = float('nan')
                E_W = float('nan')
                W_D_CPX = float('nan')
                W_D_GRT = float('nan')
        else:
            LogW = float('nan')
            E_W = float('nan')
            W_D_CPX = float('nan')
            W_D_GRT = float('nan')

        # --- U (uranium) calculations ---
        # Moved here so bulk calculations below see U_D_* values
        if not (pd.isna(log10_fO2_abs) or pd.isna(NaCl_m)):
            try:
                # Two LogU variants per user specification:
                #  - LogU_cg for CPX and GRT
                #  - LogU_rut for RUT
                LogU_cg = 2.681 + (0.1433 * float(log10_fO2_abs)) + (0.594 * float(NaCl_m))
                LogU_rut = 1.7954 + (0.1433 * float(log10_fO2_abs)) + (0.594 * float(NaCl_m))

                # Compute E_U for both variants (no molecular weight scaling)
                if not math.isfinite(LogU_cg) or LogU_cg > 308 or LogU_cg < -308:
                    E_U_cg = float('nan')
                else:
                    E_U_cg = 10.0 ** LogU_cg

                if not math.isfinite(LogU_rut) or LogU_rut > 308 or LogU_rut < -308:
                    E_U_rut = float('nan')
                else:
                    E_U_rut = 10.0 ** LogU_rut

                # Derive partition coefficients:
                # CPX and GRT use the CPX/GRT variant; RUT uses the RUT variant
                try:
                    E_for_U_cg = float(E_U_cg)
                    if not math.isfinite(E_for_U_cg) or E_for_U_cg == 0.0:
                        U_D_CPX = float('nan')
                        U_D_GRT = float('nan')
                    else:
                        U_D_CPX = 11.0 / E_for_U_cg
                        U_D_GRT = 40.0 / E_for_U_cg
                except Exception:
                    U_D_CPX = float('nan')
                    U_D_GRT = float('nan')

                try:
                    E_for_U_rut = float(E_U_rut)
                    if not math.isfinite(E_for_U_rut) or E_for_U_rut == 0.0:
                        U_D_RUT = float('nan')
                    else:
                        U_D_RUT = 94.0 / E_for_U_rut
                except Exception:
                    U_D_RUT = float('nan')

                # For backward compatibility keep LogU/E_U variables (set to RUT variant)
                LogU = float(LogU_rut) if math.isfinite(LogU_rut) else float('nan')
                E_U = float(E_U_rut) if (isinstance(E_U_rut, float) and math.isfinite(E_U_rut)) else float('nan')
            except Exception:
                LogU_cg = float('nan')
                LogU_rut = float('nan')
                E_U_cg = float('nan')
                E_U_rut = float('nan')
                LogU = float('nan')
                E_U = float('nan')
                U_D_CPX = float('nan')
                U_D_GRT = float('nan')
                U_D_RUT = float('nan')
        else:
            LogU_cg = float('nan')
            LogU_rut = float('nan')
            E_U_cg = float('nan')
            E_U_rut = float('nan')
            LogU = float('nan')
            E_U = float('nan')
            U_D_CPX = float('nan')
            U_D_GRT = float('nan')
            U_D_RUT = float('nan')

        # Compute bulk partitioning (weighted by modal proportions) for MO and other elements
        # Treat missing modal values as zero so partial modal lists still work
        try:
            def _safe_modal(x):
                return float(x) if (not pd.isna(x) and math.isfinite(float(x))) else 0.0
            m_cpx = _safe_modal(modal_cpx)
            m_grt = _safe_modal(modal_grt)
            m_rut = _safe_modal(modal_rut)

            # MO bulk
            if not (math.isfinite(MO_D_CPX) or math.isfinite(MO_D_GRT) or math.isfinite(MO_D_RUT)):
                MO_Dbulk = float('nan')
            else:
                MO_Dbulk = 0.0
                if math.isfinite(MO_D_CPX):
                    MO_Dbulk += MO_D_CPX * m_cpx
                if math.isfinite(MO_D_GRT):
                    MO_Dbulk += MO_D_GRT * m_grt
                if math.isfinite(MO_D_RUT):
                    MO_Dbulk += MO_D_RUT * m_rut
                if not math.isfinite(MO_Dbulk):
                    MO_Dbulk = float('nan')

            # Ce bulk (constants CE_D_CPX, CE_D_GRT, CE_D_RUT)
            if not (math.isfinite(CE_D_CPX) or math.isfinite(CE_D_GRT) or math.isfinite(CE_D_RUT)):
                CE_Dbulk = float('nan')
            else:
                CE_Dbulk = 0.0
                if math.isfinite(CE_D_CPX):
                    CE_Dbulk += CE_D_CPX * m_cpx
                if math.isfinite(CE_D_GRT):
                    CE_Dbulk += CE_D_GRT * m_grt
                if math.isfinite(CE_D_RUT):
                    CE_Dbulk += CE_D_RUT * m_rut
                if not math.isfinite(CE_Dbulk):
                    CE_Dbulk = float('nan')

            # W bulk (uses computed W_D_CPX/W_D_GRT and constant W_D_RUT)
            if not (math.isfinite(W_D_CPX) or math.isfinite(W_D_GRT) or math.isfinite(W_D_RUT)):
                W_Dbulk = float('nan')
            else:
                W_Dbulk = 0.0
                if math.isfinite(W_D_CPX):
                    W_Dbulk += W_D_CPX * m_cpx
                if math.isfinite(W_D_GRT):
                    W_Dbulk += W_D_GRT * m_grt
                if math.isfinite(W_D_RUT):
                    W_Dbulk += W_D_RUT * m_rut
                if not math.isfinite(W_Dbulk):
                    W_Dbulk = float('nan')

            # U bulk (uses computed U_D_CPX/U_D_GRT/U_D_RUT)
            if not (math.isfinite(U_D_CPX) or math.isfinite(U_D_GRT) or math.isfinite(U_D_RUT)):
                U_Dbulk = float('nan')
            else:
                U_Dbulk = 0.0
                if math.isfinite(U_D_CPX):
                    U_Dbulk += U_D_CPX * m_cpx
                if math.isfinite(U_D_GRT):
                    U_Dbulk += U_D_GRT * m_grt
                if math.isfinite(U_D_RUT):
                    U_Dbulk += U_D_RUT * m_rut
                if not math.isfinite(U_Dbulk):
                    U_Dbulk = float('nan')

            # Th bulk
            if not (math.isfinite(TH_D_CPX) or math.isfinite(TH_D_GRT) or math.isfinite(TH_D_RUT)):
                TH_Dbulk = float('nan')
            else:
                TH_Dbulk = 0.0
                if math.isfinite(TH_D_CPX):
                    TH_Dbulk += TH_D_CPX * m_cpx
                if math.isfinite(TH_D_GRT):
                    TH_Dbulk += TH_D_GRT * m_grt
                if math.isfinite(TH_D_RUT):
                    TH_Dbulk += TH_D_RUT * m_rut
                if not math.isfinite(TH_Dbulk):
                    TH_Dbulk = float('nan')

            # Nb bulk
            if not (math.isfinite(NB_D_CPX) or math.isfinite(NB_D_GRT) or math.isfinite(NB_D_RUT)):
                NB_Dbulk = float('nan')
            else:
                NB_Dbulk = 0.0
                if math.isfinite(NB_D_CPX):
                    NB_Dbulk += NB_D_CPX * m_cpx
                if math.isfinite(NB_D_GRT):
                    NB_Dbulk += NB_D_GRT * m_grt
                if math.isfinite(NB_D_RUT):
                    NB_Dbulk += NB_D_RUT * m_rut
                if not math.isfinite(NB_Dbulk):
                    NB_Dbulk = float('nan')

            # La bulk
            if not (math.isfinite(LA_D_CPX) or math.isfinite(LA_D_GRT) or math.isfinite(LA_D_RUT)):
                LA_Dbulk = float('nan')
            else:
                LA_Dbulk = 0.0
                if math.isfinite(LA_D_CPX):
                    LA_Dbulk += LA_D_CPX * m_cpx
                if math.isfinite(LA_D_GRT):
                    LA_Dbulk += LA_D_GRT * m_grt
                if math.isfinite(LA_D_RUT):
                    LA_Dbulk += LA_D_RUT * m_rut
                if not math.isfinite(LA_Dbulk):
                    LA_Dbulk = float('nan')
        except Exception:
            # Ensure all bulk variables exist even if the above computation failed
            MO_Dbulk = float('nan')
            CE_Dbulk = float('nan')
            W_Dbulk = float('nan')
            U_Dbulk = float('nan')
            TH_Dbulk = float('nan')
            NB_Dbulk = float('nan')
            LA_Dbulk = float('nan')

        # Compute MO fluid enmember: MO_F_EM = C0_MO / MO_Dbulk (C0_MO is MO_v from selected rock)
        # Defensive initialization of all fluid enmember variables
        MO_F_EM = CE_F_EM = W_F_EM = U_F_EM = TH_F_EM = NB_F_EM = LA_F_EM = float('nan')
        try:
            if not math.isfinite(MO_Dbulk) or MO_Dbulk == 0.0 or not math.isfinite(MO_v):
                MO_F_EM = float('nan')
            else:
                MO_F_EM = float(MO_v) / float(MO_Dbulk)
        except Exception:
            MO_F_EM = float('nan')

        # Compute fluid enmembers for other elements (CE, W, U, TH, NB, LA)
        try:
            if not math.isfinite(CE_Dbulk) or CE_Dbulk == 0.0 or not math.isfinite(Ce_v):
                CE_F_EM = float('nan')
            else:
                CE_F_EM = float(Ce_v) / float(CE_Dbulk)
        except Exception:
            CE_F_EM = float('nan')

        try:
            if not math.isfinite(W_Dbulk) or W_Dbulk == 0.0 or not math.isfinite(W_v):
                W_F_EM = float('nan')
            else:
                W_F_EM = float(W_v) / float(W_Dbulk)
        except Exception:
            W_F_EM = float('nan')

        try:
            if not math.isfinite(U_Dbulk) or U_Dbulk == 0.0 or not math.isfinite(U_v):
                U_F_EM = float('nan')
            else:
                U_F_EM = float(U_v) / float(U_Dbulk)
        except Exception:
            U_F_EM = float('nan')

        try:
            if not math.isfinite(TH_Dbulk) or TH_Dbulk == 0.0 or not math.isfinite(Th_v):
                TH_F_EM = float('nan')
            else:
                TH_F_EM = float(Th_v) / float(TH_Dbulk)
        except Exception:
            TH_F_EM = float('nan')

        try:
            if not math.isfinite(NB_Dbulk) or NB_Dbulk == 0.0 or not math.isfinite(Nb_v):
                NB_F_EM = float('nan')
            else:
                NB_F_EM = float(Nb_v) / float(NB_Dbulk)
        except Exception:
            NB_F_EM = float('nan')

        try:
            if not math.isfinite(LA_Dbulk) or LA_Dbulk == 0.0 or not math.isfinite(La_v):
                LA_F_EM = float('nan')
            else:
                LA_F_EM = float(La_v) / float(LA_Dbulk)
        except Exception:
            LA_F_EM = float('nan')

        # After computing bulks and fluid enmembers, append a result entry so the DataFrame has the requested columns
        results.append({
             "rock_id": rock_id,
             "P_bar": P_bar,
             "T_K": T_K,
             "dFMQ": dFMQ,
             "NaCl_wt_pct": wt,
             "NaCl_m": NaCl_m,
             "log10_NaCl": log10_NaCl,
             "log10_fO2_FMQ": base,
             "log10_fO2_abs": log10_fO2_abs,
             "C0_MO": MO_v,
             "C0_Ce": Ce_v,
             "C0_W": W_v,
             "C0_U": U_v,
             "C0_Th": Th_v,
             "C0_Nb": Nb_v,
             "C0_La": La_v,
             "MO_D_CPX": MO_D_CPX,
             "MO_D_GRT": MO_D_GRT,
             "MO_D_RUT": MO_D_RUT,
             "W_D_CPX": W_D_CPX,
             "W_D_GRT": W_D_GRT,
             "U_D_CPX": U_D_CPX,
             "U_D_GRT": U_D_GRT,
             "U_D_RUT": U_D_RUT,
             "MO_Dbulk": MO_Dbulk,
             "MO_F_EM": MO_F_EM,
             "W_Dbulk": W_Dbulk,
             "U_Dbulk": U_Dbulk,
             "CE_Dbulk": CE_Dbulk,
             "TH_Dbulk": TH_Dbulk,
             "NB_Dbulk": NB_Dbulk,
             "LA_Dbulk": LA_Dbulk,
             "CE_F_EM": CE_F_EM,
             "W_F_EM": W_F_EM,
             "U_F_EM": U_F_EM,
             "TH_F_EM": TH_F_EM,
             "NB_F_EM": NB_F_EM,
             "LA_F_EM": LA_F_EM,
         })

     # end for wt in sel_wts

    # Quick debug: show how many results were appended and sample keys
    try:
        print(f"DEBUG: len(results)={len(results)}")
        if len(results) > 0:
            sample_keys = list(results[0].keys())
            print(f"DEBUG: sample result keys: {sample_keys}")
    except Exception:
        pass

    res_df = pd.DataFrame(results)

    # Compute DM constants and endmember ratios from DM dict (DM defined at top of file)
    def _get_dm(key):
        try:
            return float(DM.get(key, float('nan')))
        except Exception:
            return float('nan')

    DM_MO = _get_dm("MO")
    DM_CE = _get_dm("Ce")
    DM_W = _get_dm("W")
    DM_U = _get_dm("U")
    DM_TH = _get_dm("Th")
    DM_NB = _get_dm("Nb")
    DM_LA = _get_dm("La")

    def _safe_div(a, b):
        try:
            if (not math.isfinite(a)) or (not math.isfinite(b)) or b == 0.0:
                return float('nan')
            return a / b
        except Exception:
            return float('nan')

    # DM ratios (constants)
    DM_MO_CE = _safe_div(DM_MO, DM_CE)
    DM_U_TH = _safe_div(DM_U, DM_TH)
    DM_W_TH = _safe_div(DM_W, DM_TH)
    DM_MO_W = _safe_div(DM_MO, DM_W)
    DM_NB_LA = _safe_div(DM_NB, DM_LA)

    # Build Table 4 including both DM ratios (constants) and fluid endmember ratios (F_*)
    # The fluid ratios are computed per-result row (one row per salinity). If no results exist,
    # create a single-row table with DM ratios and NaN fluid ratios.
    table4_rows = []
    if res_df is None or res_df.empty:
        table4_rows.append({
            "NaCl_wt_pct": float('nan'),
            "DM_MO_CE": DM_MO_CE,
            "DM_U_TH": DM_U_TH,
            "DM_W_TH": DM_W_TH,
            "DM_MO_W": DM_MO_W,
            "DM_NB_LA": DM_NB_LA,
            "F_MO_CE": float('nan'),
            "F_U_TH": float('nan'),
            "F_W_TH": float('nan'),
            "F_MO_W": float('nan'),
            "F_NB_LA": float('nan'),
        })
    else:
        for _, r in res_df.iterrows():
            mo_em = r.get("MO_F_EM", float('nan'))
            ce_em = r.get("CE_F_EM", float('nan'))
            u_em = r.get("U_F_EM", float('nan'))
            th_em = r.get("TH_F_EM", float('nan'))
            w_em = r.get("W_F_EM", float('nan'))
            nb_em = r.get("NB_F_EM", float('nan'))
            la_em = r.get("LA_F_EM", float('nan'))

            F_MO_CE = _safe_div(mo_em, ce_em)
            F_U_TH = _safe_div(u_em, th_em)
            F_W_TH = _safe_div(w_em, th_em)
            F_MO_W = _safe_div(mo_em, w_em)
            F_NB_LA = _safe_div(nb_em, la_em)

            table4_rows.append({
                "NaCl_wt_pct": r.get("NaCl_wt_pct"),
                "DM_MO_CE": DM_MO_CE,
                "DM_U_TH": DM_U_TH,
                "DM_W_TH": DM_W_TH,
                "DM_MO_W": DM_MO_W,
                "DM_NB_LA": DM_NB_LA,
                "F_MO_CE": F_MO_CE,
                "F_U_TH": F_U_TH,
                "F_W_TH": F_W_TH,
                "F_MO_W": F_MO_W,
                "F_NB_LA": F_NB_LA,
            })

    table4_df = pd.DataFrame(table4_rows)

    # Exact column lists per user request
    table1_cols = [
        "rock_id", "P_bar", "T_K", "dFMQ", "NaCl_wt_pct", "NaCl_m", "log10_NaCl",
        "log10_fO2_FMQ", "log10_fO2_abs",
        "C0_MO", "C0_Ce", "C0_W", "C0_U", "C0_Th", "C0_Nb", "C0_La",
    ]

    table2_cols = [
        "MO_D_CPX", "MO_D_GRT", "MO_D_RUT",
        "W_D_CPX", "W_D_GRT",
        "U_D_CPX", "U_D_GRT", "U_D_RUT",
    ]

    table3_cols = [
        "MO_Dbulk", "CE_Dbulk", "W_Dbulk", "U_Dbulk",
        "TH_Dbulk", "NB_Dbulk", "LA_Dbulk",
        "MO_F_EM", "CE_F_EM", "W_F_EM", "U_F_EM", "TH_F_EM", "NB_F_EM", "LA_F_EM",
    ]

    def _present(df, cols):
        return [c for c in cols if c in df.columns]

    with pd.option_context("display.max_columns", None, "display.width", 200):
        # Table 1
        print("\nTable 1 — Inputs, logs and input element concentrations:")
        cols1 = _present(res_df, table1_cols)
        if cols1:
            print(res_df[cols1].to_string(index=False))
        else:
            print("  (no columns available for Table 1)")

        # Table 2
        print("\nTable 2 — Per-mineral partition coefficients:")
        cols2 = _present(res_df, table2_cols)
        if cols2:
            print(res_df[cols2].to_string(index=False))
        else:
            print("  (no columns available for Table 2)")

        # Table 3
        print("\nTable 3 — Bulk (modal-weighted) partition coefficients:")
        cols3 = _present(res_df, table3_cols)
        if cols3:
            print(res_df[cols3].to_string(index=False))
        else:
            print("  (no columns available for Table 3)")

        # Table 4
        print("\nTable 4 — Endmember ratios (from DM constants):")
        try:
            print(table4_df.to_string(index=False))
        except Exception:
            print("  (no data available for Table 4)")


if __name__ == '__main__':
    main()
