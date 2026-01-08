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
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings

# Set publication-quality plot defaults
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams.update({
    'figure.figsize': (10, 8),
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

# -------- Enhanced Automation Functions --------

def create_publication_plot(all_mixes, output_dir=None, show_plot=True, save_plot=True):
    """
    Create publication-quality DM-F mixing plots with professional styling.
    
    Parameters:
    -----------
    all_mixes : list of tuples
        List of (salinity, mix_df) tuples from mixing calculations
    output_dir : str or Path, optional
        Directory to save plots. If None, saves to current directory
    show_plot : bool
        Whether to display the plot
    save_plot : bool
        Whether to save the plot to file
    
    Returns:
    --------
    fig : matplotlib.figure.Figure
        The created figure object
    """
    
    # Create figure with publication styling
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Set logarithmic scales
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Professional axis limits - optimized for geochemical data
    ax.set_ylim(0.01, 1.0)
    ax.set_xlim(0.1, 10.0)
    
    # Custom formatter for clean axis labels
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    
    def log_formatter(val, pos):
        """Custom formatter for log scales"""
        if val >= 1:
            return f'{val:.0f}'
        elif val >= 0.1:
            return f'{val:.1f}'
        else:
            return f'{val:.2f}'
    
    formatter = FuncFormatter(log_formatter)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    
    # Set up professional grid
    ax.grid(True, which="major", linestyle='-', linewidth=0.8, alpha=0.3)
    ax.grid(True, which="minor", linestyle='--', linewidth=0.5, alpha=0.2)
    
    # Color scheme for different salinities
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    salinity_colors = {}
    
    # f values to mark and their labels
    label_fs = [0.01, 0.05, 0.1]
    label_percent = {0.01: '1%', 0.05: '5%', 0.1: '10%'}
    
    # Get axis limits for bounds checking
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    
    # Plot mixing curves for each salinity
    for i, (sal, mix_df) in enumerate(all_mixes):
        if mix_df is None or mix_df.empty:
            continue
            
        color = colors[i % len(colors)]
        salinity_colors[sal] = color
        
        # Filter for valid data points
        valid_mask = (np.isfinite(mix_df['U_TH_MIX']) & 
                     np.isfinite(mix_df['MO_CE_MIX']) & 
                     (mix_df['U_TH_MIX'] > 0) & 
                     (mix_df['MO_CE_MIX'] > 0))
        
        if not np.any(valid_mask):
            continue
            
        plot_df = mix_df.loc[valid_mask].copy()
        plot_df = plot_df.sort_values('U_TH_MIX')
        
        x = plot_df['U_TH_MIX'].values
        y = plot_df['MO_CE_MIX'].values
        
        # Plot mixing curve
        label = f'{sal:.1f}% NaCl' if not pd.isna(sal) else 'Unknown salinity'
        ax.plot(x, y, color=color, linestyle='-', linewidth=2.5, 
               label=label, alpha=0.8)
        
        # Add fraction markers
        for f in label_fs:
            sel = plot_df[np.isclose(plot_df['f'].values, f, atol=1e-6)]
            if len(sel) == 0:
                continue
                
            try:
                sx = float(sel['U_TH_MIX'].iloc[0])
                sy = float(sel['MO_CE_MIX'].iloc[0])
            except (IndexError, ValueError):
                continue
                
            if not (math.isfinite(sx) and math.isfinite(sy)):
                continue
                
            # Only plot if within bounds
            if x_min <= sx <= x_max and y_min <= sy <= y_max:
                ax.plot(sx, sy, marker='s', color=color, markersize=6, 
                       markeredgewidth=1, markeredgecolor='white', alpha=0.9)
        
        # Plot DM point (f = 0)
        try:
            sel_dm = plot_df[np.isclose(plot_df['f'].values, 0.0, atol=1e-6)]
            if len(sel_dm) > 0:
                sxd = float(sel_dm['U_TH_MIX'].iloc[0])
                syd = float(sel_dm['MO_CE_MIX'].iloc[0])
                
                if (math.isfinite(sxd) and math.isfinite(syd) and 
                    x_min <= sxd <= x_max and y_min <= syd <= y_max):
                    ax.plot(sxd, syd, marker='o', markerfacecolor='#8B0000', 
                           markeredgecolor='white', markersize=10, 
                           markeredgewidth=2, alpha=0.9, zorder=10)
        except (IndexError, ValueError):
            pass
    
    # Professional labeling
    ax.set_xlabel('U/Th', fontweight='bold')
    ax.set_ylabel('Mo/Ce', fontweight='bold')
    ax.set_title('Depleted Mantle - Fluid Mixing Model\n'
                'Trace Element Signatures vs Salinity', 
                fontweight='bold', pad=20)
    
    # Professional legend
    legend = ax.legend(title='Salinity', loc='upper right', 
                      frameon=True, fancybox=True, shadow=True)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.9)
    
    # Add DM point to legend
    ax.plot([], [], marker='o', markerfacecolor='#8B0000', 
           markeredgecolor='white', markersize=8, markeredgewidth=1.5,
           linestyle='None', label='Depleted Mantle', alpha=0.9)
    
    # Update legend
    legend = ax.legend(title='Salinity', loc='upper right', 
                      frameon=True, fancybox=True, shadow=True)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.9)
    
    plt.tight_layout()
    
    # Save plot if requested
    if save_plot:
        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(exist_ok=True)
        else:
            output_path = Path.cwd()
            
        filename = output_path / 'DM_fluid_mixing_model.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"Plot saved to: {filename}")
        
        # Also save as PDF for publications
        pdf_filename = output_path / 'DM_fluid_mixing_model.pdf'
        plt.savefig(pdf_filename, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"PDF version saved to: {pdf_filename}")
    
    if show_plot:
        plt.show()
    
    return fig

def automated_batch_analysis(df, output_dir=None, create_plots=True):
    """
    Fully automated batch analysis for all salinity conditions.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Normalized input dataframe
    output_dir : str or Path, optional
        Directory for outputs
    create_plots : bool
        Whether to generate plots
        
    Returns:
    --------
    dict : Complete results with all tables and plots
    """
    
    print("="*60)
    print("AUTOMATED DM-FLUID MIXING ANALYSIS")
    print("="*60)
    
    # Standard salinity range for your workflow
    sal_options = [0.01, 5.0, 10.0, 15.0, 20.0]
    
    print(f"Processing {len(df)} samples across {len(sal_options)} salinity conditions...")
    print(f"Salinity range: {sal_options} wt% NaCl")
    
    # Run vectorized processing
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Suppress minor numerical warnings
        results_df = process_all_vectorized(df, sal_options)
    
    if results_df.empty:
        print("ERROR: No valid results generated!")
        return None
    
    print(f"✓ Generated {len(results_df)} result rows")
    
    # Organize results by tables
    results = {
        'raw_data': results_df,
        'summary_stats': {},
        'tables': {},
        'plots': {}
    }
    
    # Generate summary statistics
    numeric_cols = results_df.select_dtypes(include=[np.number]).columns
    results['summary_stats'] = {
        'mean': results_df[numeric_cols].mean(),
        'std': results_df[numeric_cols].std(),
        'min': results_df[numeric_cols].min(),
        'max': results_df[numeric_cols].max()
    }
    
    # Create publication tables
    table_configs = {
        'inputs': ['rock_id', 'P_bar', 'T_K', 'dFMQ', 'NaCl_wt_pct', 'NaCl_m', 
                  'log10_NaCl', 'log10_fO2_FMQ', 'log10_fO2_abs'],
        'concentrations': ['C0_MO', 'C0_Ce', 'C0_W', 'C0_U', 'C0_Th', 'C0_Nb', 'C0_La'],
        'partitioning': ['MO_D_CPX', 'MO_D_GRT', 'MO_D_RUT', 'W_D_CPX', 'W_D_GRT',
                        'U_D_CPX', 'U_D_GRT', 'U_D_RUT'],
        'bulk_coefficients': ['MO_Dbulk', 'CE_Dbulk', 'W_Dbulk', 'U_Dbulk',
                             'TH_Dbulk', 'NB_Dbulk', 'LA_Dbulk'],
        'fluid_endmembers': ['MO_F_EM', 'CE_F_EM', 'W_F_EM', 'U_F_EM', 
                           'TH_F_EM', 'NB_F_EM', 'LA_F_EM']
    }
    
    for table_name, cols in table_configs.items():
        available_cols = [c for c in cols if c in results_df.columns]
        if available_cols:
            results['tables'][table_name] = results_df[available_cols].copy()
    
    # Generate mixing analysis
    print("\n" + "="*40)
    print("GENERATING MIXING ANALYSIS")
    print("="*40)
    
    all_mixes = []
    for sal in sal_options:
        sal_data = results_df[results_df['NaCl_wt_pct'] == sal].iloc[0] if len(results_df[results_df['NaCl_wt_pct'] == sal]) > 0 else None
        if sal_data is None:
            continue
            
        mix_df = generate_mixing_table(sal_data)
        if mix_df is not None:
            all_mixes.append((sal, mix_df))
            print(f"✓ Generated mixing table for {sal}% NaCl")
    
    results['mixing_tables'] = all_mixes
    
    # Create plots if requested
    if create_plots and all_mixes:
        print("\n" + "="*40)
        print("CREATING PUBLICATION PLOTS")
        print("="*40)
        
        try:
            fig = create_publication_plot(all_mixes, output_dir=output_dir)
            results['plots']['main_mixing'] = fig
            print("✓ Publication-quality mixing plot created")
        except Exception as e:
            print(f"⚠ Plot creation failed: {e}")
    
    # Save results if output directory specified
    if output_dir:
        save_results_to_files(results, output_dir)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    
    return results

def generate_mixing_table(sal_data):
    """Generate mixing table for a single salinity condition."""
    
    # DM constants from your model
    DM = {"MO": 0.025, "Ce": 0.772, "W": 0.0024, "U": 0.0047, 
          "Th": 0.0137, "Nb": 0.21, "La": 0.234}
    
    f_values = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
    
    # Extract fluid endmembers
    endmembers = {
        'MO': sal_data.get('MO_F_EM', np.nan),
        'Ce': sal_data.get('CE_F_EM', np.nan),
        'W': sal_data.get('W_F_EM', np.nan),
        'U': sal_data.get('U_F_EM', np.nan),
        'Th': sal_data.get('TH_F_EM', np.nan),
        'Nb': sal_data.get('NB_F_EM', np.nan),
        'La': sal_data.get('LA_F_EM', np.nan)
    }
    
    def safe_mix(f_em, dm_val, frac):
        """Safely calculate mixing value."""
        try:
            if not np.isfinite(f_em):
                return np.nan
            return (float(f_em) * frac) + (float(dm_val) * (1.0 - frac))
        except:
            return np.nan
    
    def safe_ratio(a, b):
        """Safely calculate ratio."""
        try:
            if not (np.isfinite(a) and np.isfinite(b)) or b == 0:
                return np.nan
            return float(a) / float(b)
        except:
            return np.nan
    
    mix_rows = []
    for frac in f_values:
        # Calculate mixed compositions
        mo_mix = safe_mix(endmembers['MO'], DM['MO'], frac)
        ce_mix = safe_mix(endmembers['Ce'], DM['Ce'], frac)
        w_mix = safe_mix(endmembers['W'], DM['W'], frac)
        u_mix = safe_mix(endmembers['U'], DM['U'], frac)
        th_mix = safe_mix(endmembers['Th'], DM['Th'], frac)
        nb_mix = safe_mix(endmembers['Nb'], DM['Nb'], frac)
        la_mix = safe_mix(endmembers['La'], DM['La'], frac)
        
        # Calculate elemental ratios
        mo_ce = safe_ratio(mo_mix, ce_mix)
        u_th = safe_ratio(u_mix, th_mix)
        w_th = safe_ratio(w_mix, th_mix)
        mo_w = safe_ratio(mo_mix, w_mix)
        nb_la = safe_ratio(nb_mix, la_mix)
        
        mix_rows.append({
            'f': frac,
            'MO_DM_F_MIX': mo_mix,
            'CE_DM_F_MIX': ce_mix,
            'W_DM_F_MIX': w_mix,
            'U_DM_F_MIX': u_mix,
            'TH_DM_F_MIX': th_mix,
            'NB_DM_F_MIX': nb_mix,
            'LA_DM_F_MIX': la_mix,
            'MO_CE_MIX': mo_ce,
            'U_TH_MIX': u_th,
            'W_TH_MIX': w_th,
            'MO_W_MIX': mo_w,
            'NB_LA_MIX': nb_la,
        })
    
    return pd.DataFrame(mix_rows)

def save_results_to_files(results, output_dir):
    """Save all results to organized files."""
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Save main results
    if 'raw_data' in results:
        csv_path = output_path / 'complete_results.csv'
        results['raw_data'].to_csv(csv_path, index=False, float_format='%.6g')
        print(f"✓ Complete results saved to: {csv_path}")
    
    # Save individual tables
    if 'tables' in results:
        for table_name, table_df in results['tables'].items():
            table_path = output_path / f'table_{table_name}.csv'
            table_df.to_csv(table_path, index=False, float_format='%.6g')
            print(f"✓ {table_name} table saved to: {table_path}")
    
    # Save mixing tables
    if 'mixing_tables' in results:
        mixing_dir = output_path / 'mixing_tables'
        mixing_dir.mkdir(exist_ok=True)
        
        for sal, mix_df in results['mixing_tables']:
            mix_path = mixing_dir / f'mixing_table_{sal}pct_NaCl.csv'
            mix_df.to_csv(mix_path, index=False, float_format='%.6g')
        
        print(f"✓ Mixing tables saved to: {mixing_dir}")
    
    # Save summary statistics
    if 'summary_stats' in results:
        stats_path = output_path / 'summary_statistics.csv'
        stats_df = pd.DataFrame(results['summary_stats'])
        stats_df.to_csv(stats_path, float_format='%.6g')
        print(f"✓ Summary statistics saved to: {stats_path}")

# -------- constants --------
DM = {"MO": 0.025, "Ce": 0.772, "W": 0.0024, "U": 0.0047, "Th": 0.0137, "Nb": 0.21, "La": 0.234}  # ppm
M_MO = 95.95
M_W = 183.84
M_NACL = 58.44
CE_D_CPX, CE_D_GRT, CE_D_RUT = 2.0, 0.4, 2.0
TH_D_CPX, TH_D_GRT, TH_D_RUT = 1.19, 0.610, 0.10
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

# helper: convert NaCl weight percent to molality (mol/kg solvent)
def wt_to_molality(wt_pct):
    try:
        w = float(wt_pct)
        if not math.isfinite(w):
            return float('nan')
        # 0 wt% has molality 0 -> treat explicitly
        if w <= 0.0 or w >= 100.0:
            return float('nan')
        return 1000.0 * w / (M_NACL * (100.0 - w))
    except Exception:
        return float('nan')

# -------- header normalization helpers --------
def _norm(s: str) -> str:
    if not isinstance(s, str): return ""
    x = s.strip().lower()
    # Normalize common delta characters to ascii 'd' so headers like 'ΔFMQ' map to 'dfmq'
    x = x.replace('Δ', 'd').replace('δ', 'd')
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
    # modal proportion fallbacks (normalized keys)
    "modegrt": "mode_grt", "grt": "mode_grt", "garnet": "mode_grt",
    "modecpx": "mode_cpx", "cpx": "mode_cpx", "clinopyroxene": "mode_cpx",
    "moderut": "mode_rut", "rut": "mode_rut", "rutile": "mode_rut",
    # simple trace element column mappings (if present as exact names)
    "mo": "C0_MO", "molibdenum": "C0_MO",
    "ce": "C0_Ce", "cerium": "C0_Ce",
    "w": "C0_W", "tungsten": "C0_W",
    "u": "C0_U", "uranium": "C0_U",
    "th": "C0_Th", "thorium": "C0_Th",
    "nb": "C0_Nb", "niobium": "C0_Nb",
    "la": "C0_La", "lanthanum": "C0_La",
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

    # Try to populate modal and trace columns from flexible header variants using norm_map
    # This keeps backward compatibility with the row-wise find_val approach but creates proper
    # dataframe columns for vectorized processing.
    def _assign_first_variant(variants, colname):
        if colname in out.columns:
            return
        for col, key in norm_map.items():
            if key in variants or any(v in key for v in variants):
                out[colname] = df[col]
                return

    _assign_first_variant({'modegrt', 'grt', 'garnet'}, 'mode_grt')
    _assign_first_variant({'modecpx', 'cpx', 'clinopyroxene'}, 'mode_cpx')
    _assign_first_variant({'moderut', 'rut', 'rutile'}, 'mode_rut')

    _assign_first_variant({'mo', 'molibdenum', 'c0moppm', 'c0_mo_ppm', 'c0_mo'}, 'C0_MO')
    _assign_first_variant({'ce', 'cerium', 'c0ceppm', 'c0_ce_ppm', 'c0_ce'}, 'C0_Ce')
    _assign_first_variant({'w', 'tungsten', 'c0wppm', 'c0_w_ppm', 'c0_w'}, 'C0_W')
    _assign_first_variant({'u', 'uranium', 'c0uppm', 'c0_u_ppm', 'c0_u'}, 'C0_U')
    _assign_first_variant({'th', 'thorium', 'c0thppm', 'c0_th_ppm', 'c0_th'}, 'C0_Th')
    _assign_first_variant({'nb', 'niobium', 'c0nbppm', 'c0_nb_ppm', 'c0_nb'}, 'C0_Nb')
    _assign_first_variant({'la', 'lanthanum', 'c0lappm', 'c0_la_ppm', 'c0_la'}, 'C0_La')

    # Coerce numerics (handle unicode minus and stray text) - robust for negative dFMQ
    for c in ("P_bar", "T_K", "dFMQ", 'mode_grt', 'mode_cpx', 'mode_rut', 'C0_MO', 'C0_Ce', 'C0_W', 'C0_U', 'C0_Th', 'C0_Nb', 'C0_La'):
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

    This implementation validates T_K and returns NaN for invalid/non-positive T_K values
    to avoid unexpected exceptions or infinite results.
    """
    try:
        t = float(T_K)
        if not math.isfinite(t) or t <= 0.0:
            return float('nan')
    except Exception:
        return float('nan')
    try:
        p_bar = float(P_bar)
    except Exception:
        p_bar = float('nan')
    return (-25096.3 / t) + 8.735 + 0.11 * ((p_bar - 1.0) / t)


# Vectorized helper for FMQ base values (returns numpy array-like)
def compute_base_fmq_log10_vectorized(T_K_arr, P_bar_arr):
    """
    Vectorized FMQ base calculation. Returns an array of float values with NaN
    where inputs are invalid (non-finite or non-positive T_K).
    """
    TK = np.array(T_K_arr, dtype=float)
    PB = np.array(P_bar_arr, dtype=float)
    out = np.full_like(TK, np.nan, dtype=float)
    mask = np.isfinite(TK) & (TK > 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        out[mask] = (-25096.3 / TK[mask]) + 8.735 + 0.11 * ((PB[mask] - 1.0) / TK[mask])
    return out

# -------- process_all_vectorized --------
def process_all_vectorized(df: pd.DataFrame, sal_wts):
    """
    Vectorized (but dtype-safe) processing for all rows and a list of salinity wt.% values.
    Returns a pandas DataFrame with one row per input row per salinity.
    This version builds a per-salinity concatenation to preserve pandas dtypes
    and then coerces numeric columns explicitly to avoid object/string issues.
    """
    if df is None or df.empty:
        return pd.DataFrame()

    # Work on a copy
    base_df = df.copy()

    # Ensure expected columns exist
    expected = ['mode_grt', 'mode_cpx', 'mode_rut', 'C0_MO', 'C0_Ce', 'C0_W', 'C0_U', 'C0_Th', 'C0_Nb', 'C0_La', 'dFMQ']
    for c in expected:
        if c not in base_df.columns:
            base_df[c] = np.nan

    # Compute FMQ base and provisional absolute fO2 (dFMQ may be NaN)
    # Use the shared, validated vectorized helper to ensure consistent behavior
    base_df['log10_fO2_FMQ'] = compute_base_fmq_log10_vectorized(base_df['T_K'], base_df['P_bar'])
    base_df['log10_fO2_abs'] = base_df['log10_fO2_FMQ'] + base_df['dFMQ'].fillna(0.0)
    base_df.loc[base_df['dFMQ'].isna(), 'log10_fO2_abs'] = np.nan

    # Build per-salinity frames and concatenate to preserve dtypes
    frames = []
    for wt in sal_wts:
        tmp = base_df.copy()
        tmp['NaCl_wt_pct'] = float(wt)
        frames.append(tmp)
    rep_df = pd.concat(frames, ignore_index=True)

    # Coerce numeric columns to float explicitly to avoid object dtypes
    num_cols = ['P_bar', 'T_K', 'dFMQ', 'mode_grt', 'mode_cpx', 'mode_rut',
                'C0_MO', 'C0_Ce', 'C0_W', 'C0_U', 'C0_Th', 'C0_Nb', 'C0_La',
                'log10_fO2_FMQ', 'log10_fO2_abs']
    for c in num_cols:
        if c in rep_df.columns:
            rep_df[c] = pd.to_numeric(rep_df[c], errors='coerce')

    # molality and log10_NaCl
    wt_series = rep_df['NaCl_wt_pct'].astype(float)
    rep_df['NaCl_m'] = 1000.0 * wt_series / (float(M_NACL) * (100.0 - wt_series))
    rep_df['NaCl_m'].replace([np.inf, -np.inf], np.nan, inplace=True)
    rep_df['log10_NaCl'] = np.where((rep_df['NaCl_m'] > 0) & np.isfinite(rep_df['NaCl_m']), np.log10(rep_df['NaCl_m']), np.nan)

    # LogMO
    rep_df['LogMO'] = (0.435 * rep_df['log10_fO2_abs']) + (0.42 * rep_df['log10_NaCl']) - (1.8 * (1000.0 / rep_df['T_K'])) + 4.8

    # E_MO with overflow guards
    scale_log10 = float(SCALE_LOG10_M_MO) if np.isfinite(SCALE_LOG10_M_MO) else np.nan
    combined_log = rep_df['LogMO'] + scale_log10
    ok = np.isfinite(combined_log) & (combined_log <= 308) & (combined_log >= -308)
    rep_df['E_MO'] = np.nan
    rep_df.loc[ok, 'E_MO'] = (10.0 ** rep_df.loc[ok, 'LogMO']) * M_MO * 1000.0
    rep_df['MO_D_CPX'] = np.where(rep_df['E_MO'] > 0, 40.0 / rep_df['E_MO'], np.nan)
    rep_df['MO_D_GRT'] = np.where(rep_df['E_MO'] > 0, 12.0 / rep_df['E_MO'], np.nan)
    rep_df['MO_D_RUT'] = np.where(rep_df['E_MO'] > 0, 87670.0 / rep_df['E_MO'], np.nan)

    # LogW and W partitioning
    rep_df['LogW'] = (0.07 * rep_df['log10_fO2_abs']) - (4.7236 * (1000.0 / rep_df['T_K'])) + 4.4271
    scale_log10_w = float(SCALE_LOG10_M_W) if np.isfinite(SCALE_LOG10_M_W) else np.nan
    combined_log_w = rep_df['LogW'] + scale_log10_w
    okw = np.isfinite(combined_log_w) & (combined_log_w <= 308) & (combined_log_w >= -308)
    rep_df['E_W'] = np.nan
    rep_df.loc[okw, 'E_W'] = (10.0 ** rep_df.loc[okw, 'LogW']) * M_W * 1000.0
    rep_df['W_D_CPX'] = np.where(rep_df['E_W'] > 0, 60.0 / rep_df['E_W'], np.nan)
    rep_df['W_D_GRT'] = np.where(rep_df['E_W'] > 0, 12.0 / rep_df['E_W'], np.nan)

    # LogU & U partitioning
    rep_df['LogU_cg'] = 2.681 + (0.1433 * rep_df['log10_fO2_abs']) + (0.594 * rep_df['NaCl_m'])
    rep_df['LogU_rut'] = 1.7954 + (0.1433 * rep_df['log10_fO2_abs']) + (0.594 * rep_df['NaCl_m'])
    rep_df['E_U_cg'] = np.where(np.isfinite(rep_df['LogU_cg']) & (rep_df['LogU_cg'] <= 308) & (rep_df['LogU_cg'] >= -308), 10.0 ** rep_df['LogU_cg'], np.nan)
    rep_df['E_U_rut'] = np.where(np.isfinite(rep_df['LogU_rut']) & (rep_df['LogU_rut'] <= 308) & (rep_df['LogU_rut'] >= -308), 10.0 ** rep_df['LogU_rut'], np.nan)
    rep_df['U_D_CPX'] = np.where(rep_df['E_U_cg'] > 0, 11.0 / rep_df['E_U_cg'], np.nan)
    rep_df['U_D_GRT'] = np.where(rep_df['E_U_cg'] > 0, 40.0 / rep_df['E_U_cg'], np.nan)
    rep_df['U_D_RUT'] = np.where(rep_df['E_U_rut'] > 0, 94.0 / rep_df['E_U_rut'], np.nan)

    # Compute bulk coefficient (modal-weighted)
    for m in ('mode_cpx', 'mode_grt', 'mode_rut'):
        if m in rep_df.columns:
            rep_df[m] = rep_df[m].astype(float).fillna(0.0)
    m_cpx = rep_df.get('mode_cpx', pd.Series(0.0, index=rep_df.index))
    m_grt = rep_df.get('mode_grt', pd.Series(0.0, index=rep_df.index))
    m_rut = rep_df.get('mode_rut', pd.Series(0.0, index=rep_df.index))

    rep_df['MO_Dbulk'] = (rep_df['MO_D_CPX'].fillna(0.0) * m_cpx) + (rep_df['MO_D_GRT'].fillna(0.0) * m_grt) + (rep_df['MO_D_RUT'].fillna(0.0) * m_rut)
    rep_df['CE_Dbulk'] = (CE_D_CPX * m_cpx) + (CE_D_GRT * m_grt) + (CE_D_RUT * m_rut)
    rep_df['W_Dbulk'] = (rep_df.get('W_D_CPX', 0.0).fillna(0.0) * m_cpx) + (rep_df.get('W_D_GRT', 0.0).fillna(0.0) * m_grt) + (W_D_RUT * m_rut)
    rep_df['U_Dbulk'] = (rep_df.get('U_D_CPX', 0.0).fillna(0.0) * m_cpx) + (rep_df.get('U_D_GRT', 0.0).fillna(0.0) * m_grt) + (rep_df.get('U_D_RUT', 0.0).fillna(0.0) * m_rut)
    rep_df['TH_Dbulk'] = (TH_D_CPX * m_cpx) + (TH_D_GRT * m_grt) + (TH_D_RUT * m_rut)
    rep_df['NB_Dbulk'] = (NB_D_CPX * m_cpx) + (NB_D_GRT * m_grt) + (NB_D_RUT * m_rut)
    rep_df['LA_Dbulk'] = (LA_D_CPX * m_cpx) + (LA_D_GRT * m_grt) + (LA_D_RUT * m_rut)

    # safe division helper
    def safe_div(a, b):
        a = np.array(a, dtype=float)
        b = np.array(b, dtype=float)
        out = np.full_like(a, np.nan, dtype=float)
        mask = np.isfinite(a) & np.isfinite(b) & (b != 0)
        out[mask] = a[mask] / b[mask]
        return out

    rep_df['MO_F_EM'] = safe_div(rep_df.get('C0_MO', np.nan), rep_df['MO_Dbulk'])
    rep_df['CE_F_EM'] = safe_div(rep_df.get('C0_Ce', np.nan), rep_df['CE_Dbulk'])
    rep_df['W_F_EM'] = safe_div(rep_df.get('C0_W', np.nan), rep_df['W_Dbulk'])
    rep_df['U_F_EM'] = safe_div(rep_df.get('C0_U', np.nan), rep_df['U_Dbulk'])
    rep_df['TH_F_EM'] = safe_div(rep_df.get('C0_Th', np.nan), rep_df['TH_Dbulk'])
    rep_df['NB_F_EM'] = safe_div(rep_df.get('C0_Nb', np.nan), rep_df['NB_Dbulk'])
    rep_df['LA_F_EM'] = safe_div(rep_df.get('C0_La', np.nan), rep_df['LA_Dbulk'])

    # select and return relevant columns
    keep_cols = [
        'rock_id', 'P_bar', 'T_K', 'dFMQ', 'NaCl_wt_pct', 'NaCl_m', 'log10_NaCl',
        'log10_fO2_FMQ', 'log10_fO2_abs',
        'C0_MO', 'C0_Ce', 'C0_W', 'C0_U', 'C0_Th', 'C0_Nb', 'C0_La',
        'MO_D_CPX', 'MO_D_GRT', 'MO_D_RUT', 'W_D_CPX', 'W_D_GRT', 'U_D_CPX', 'U_D_GRT', 'U_D_RUT',
        'MO_Dbulk', 'CE_Dbulk', 'W_Dbulk', 'U_Dbulk', 'TH_Dbulk', 'NB_Dbulk', 'LA_Dbulk',
        'MO_F_EM', 'CE_F_EM', 'W_F_EM', 'U_F_EM', 'TH_F_EM', 'NB_F_EM', 'LA_F_EM'
    ]
    keep = [c for c in keep_cols if c in rep_df.columns]
    out = rep_df[keep].copy()
    out['rock_id'] = out['rock_id'].astype(str)
    return out

# -------- main --------
def main():
    ap = argparse.ArgumentParser(description="Automated DM-Fluid mixing analysis with publication-quality outputs.")
    ap.add_argument("--in_xlsx", default="/Users/polsuka/Desktop/Bali et al 2012 model/input.xlsx")
    ap.add_argument("--sheet", default=0)
    # Enhanced automation options
    ap.add_argument("--auto", action="store_true", help="Run fully automated analysis for all samples")
    ap.add_argument("--output_dir", default=None, help="Directory for organized output files and plots")
    ap.add_argument("--no_plots", action="store_true", help="Skip plot generation")
    # Original single-sample options (preserved for backward compatibility)
    ap.add_argument("--row", type=int, help="Excel row number to use (start at 2). Overrides interactive prompt.")
    ap.add_argument("--salinity-index", type=int, choices=range(1, 6), help="Select salinity option by index (1-5).")
    ap.add_argument("--salinity-wt", type=float, help="Directly specify a single salinity wt.% NaCl to use (overrides index).")
    ap.add_argument("--batch", action="store_true", help="Process all rows in the input file in batch (non-interactive).")
    ap.add_argument("--out_csv", default=None, help="Optional path to write results CSV when running in batch mode.")
    args = ap.parse_args()

    # Load and normalize data
    try:
        df_raw = load_excel(args.in_xlsx, args.sheet)
        df = normalize_headers(df_raw)
        print(f"✓ Loaded {len(df)} samples from {args.in_xlsx}")
    except Exception as e:
        print(f"ERROR: Failed to load input file: {e}")
        return

    # NEW: Fully automated workflow
    if args.auto:
        print("\n🚀 RUNNING AUTOMATED ANALYSIS")
        print("This will process all samples with all standard salinity conditions...")
        
        results = automated_batch_analysis(
            df, 
            output_dir=args.output_dir,
            create_plots=not args.no_plots
        )
        
        if results:
            print("\n✅ AUTOMATED ANALYSIS COMPLETED SUCCESSFULLY!")
            print("\nResults generated:")
            if 'tables' in results:
                for table_name in results['tables'].keys():
                    print(f"  • {table_name.replace('_', ' ').title()} table")
            if 'mixing_tables' in results:
                print(f"  • Mixing tables for {len(results['mixing_tables'])} salinity conditions")
            if 'plots' in results and results['plots']:
                print(f"  • {len(results['plots'])} publication-quality plot(s)")
            
            if args.output_dir:
                print(f"\n📁 All outputs saved to: {args.output_dir}")
        else:
            print("❌ Automated analysis failed!")
        return

    # Original workflow for single-sample analysis (preserved)
    sal_options = [0.01, 5.0, 10.0, 15.0, 20.0]
    if args.batch:
        # Determine salinity list for batch: explicit wt -> index -> default to all five options
        if args.salinity_wt is not None:
            sal_list = [args.salinity_wt]
        elif args.salinity_index is not None:
            if 1 <= args.salinity_index <= len(sal_options):
                sal_list = [sal_options[args.salinity_index - 1]]
            else:
                print(f"Invalid salinity-index {args.salinity_index}; using default salinity options.")
                sal_list = sal_options[:]
        else:
            sal_list = sal_options[:]

        res_df = process_all_vectorized(df, sal_list)
        float_fmt = lambda x: f"{x:.8g}" if (isinstance(x, (int, float, np.floating, np.integer)) and not pd.isna(x)) else str(x)
        with pd.option_context("display.max_columns", None, "display.width", 200):
            print(res_df.to_string(index=False, float_format=float_fmt))
        if args.out_csv:
            try:
                res_df.to_csv(args.out_csv, index=False, float_format='%.8g')
                print(f"Wrote results to {args.out_csv}")
            except Exception as e:
                print(f"Failed to write CSV: {e}")
        return

    # Show table with Excel-style row numbers (header=1, first data row=2)
    disp = df.copy()
    disp.insert(0, "Excel_row", range(2, 2 + len(disp)))
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
    # Prefer normalized/canonical column names produced by normalize_headers; fall back to legacy find_val when needed
    modal_grt = row.get('mode_grt') if 'mode_grt' in row.index else find_val(row, {"modegrt", "grt", "garnet"})
    modal_cpx = row.get('mode_cpx') if 'mode_cpx' in row.index else find_val(row, {"modecpx", "cpx", "clinopyroxene"})
    modal_rut = row.get('mode_rut') if 'mode_rut' in row.index else find_val(row, {"moderut", "rut", "rutile"})

    # trace elements: prefer C0_ canonical columns
    MO_v = row.get('C0_MO') if 'C0_MO' in row.index else find_val(row, {"mo", "molibdenum", "c0moppm", "c0_mo_ppm", "c0_mo"})
    Ce_v = row.get('C0_Ce') if 'C0_Ce' in row.index else find_val(row, {"ce", "cerium", "c0ceppm", "c0_ce_ppm", "c0_ce"})
    W_v = row.get('C0_W') if 'C0_W' in row.index else find_val(row, {"w", "tungsten", "c0wppm", "c0_w_ppm", "c0_w"})
    Th_v = row.get('C0_Th') if 'C0_Th' in row.index else find_val(row, {"th", "thorium", "c0thppm", "c0_th_ppm", "c0_th"})
    Nb_v = row.get('C0_Nb') if 'C0_Nb' in row.index else find_val(row, {"nb", "niobium", "c0nbppm", "c0_nb_ppm", "c0_nb"})
    La_v = row.get('C0_La') if 'C0_La' in row.index else find_val(row, {"la", "lanthanum", "c0lappm", "c0_la_ppm", "c0_la"})
    U_v = row.get('C0_U') if 'C0_U' in row.index else find_val(row, {"u", "uranium", "c0uppm", "c0_u_ppm", "c0_u"})

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
    # (console salinity menu removed to reduce clutter)

    # Determine salinity selection (non-interactive): compute all five unless overrides provided
    if args.salinity_wt is not None:
        sel_wts = [args.salinity_wt]
    elif args.salinity_index is not None:
        sel_sal_i = args.salinity_index
        if 1 <= sel_sal_i <= len(sal_options):
            sel_wts = [sal_options[sel_sal_i - 1]]
        else:
            print(f"Invalid salinity-index {sel_sal_i}, defaulting to all options.")
            sel_wts = sal_options[:]
    else:
        # Default behavior: compute all five standard salinities
        sel_wts = sal_options[:]

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
        except Exception as e:
            # Ensure all bulk variables exist even if the above computation failed
            import traceback
            traceback.print_exc()
            print(f"DEBUG: bulk partitioning computation failed: {e!r}")
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
            print(res_df[cols1].to_string(index=False, float_format=lambda x: f"{x:.8g}"))
        else:
            print("  (no columns available for Table 1)")

        # Table 2
        print("\nTable 2 — Per-mineral partition coefficients:")
        cols2 = _present(res_df, table2_cols)
        if cols2:
            print(res_df[cols2].to_string(index=False, float_format=lambda x: f"{x:.8g}"))
        else:
            print("  (no columns available for Table 2)")

        # Table 3
        print("\nTable 3 — Bulk (modal-weighted) partition coefficients:")
        cols3 = _present(res_df, table3_cols)
        if cols3:
            print(res_df[cols3].to_string(index=False, float_format=lambda x: f"{x:.8g}"))
        else:
            print("  (no columns available for Table 3)")

        # Table 4
        print("\nTable 4 — Endmember ratios (from DM constants):")
        try:
            print(table4_df.to_string(index=False, float_format=lambda x: f"{x:.8g}"))
        except Exception:
            print("  (no data available for Table 4)")

        # Table 5 — build DM-F mixing table per-salinity row using per-element fluid endmembers
        f_values = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1]

        def _safe_div_val(a, b):
            try:
                if not (math.isfinite(a) and math.isfinite(b)):
                    return float('nan')
                if float(b) == 0.0:
                    return float('nan')
                return float(a) / float(b)
            except Exception:
                return float('nan')

        # For each salinity/result row produce a mix table using MIX = (F_EM * f) + (DM * (1-f))
        all_mixes = []
        for _, r in res_df.iterrows():
            sal = r.get("NaCl_wt_pct")
            mo_em = r.get("MO_F_EM", float('nan'))
            ce_em = r.get("CE_F_EM", float('nan'))
            w_em = r.get("W_F_EM", float('nan'))
            u_em = r.get("U_F_EM", float('nan'))
            th_em = r.get("TH_F_EM", float('nan'))
            nb_em = r.get("NB_F_EM", float('nan'))
            la_em = r.get("LA_F_EM", float('nan'))

            mix_rows = []
            for frac in f_values:
                # compute per-element mix only when fluid endmember is finite; otherwise result is NaN
                def mix_val(f_em, dm_const):
                    try:
                        if not math.isfinite(f_em):
                            return float('nan')
                        return (float(f_em) * frac) + (float(dm_const) * (1.0 - frac))
                    except Exception:
                        return float('nan')

                mo_mix = mix_val(mo_em, DM_MO)
                ce_mix = mix_val(ce_em, DM_CE)
                w_mix = mix_val(w_em, DM_W)
                u_mix = mix_val(u_em, DM_U)
                th_mix = mix_val(th_em, DM_TH)
                nb_mix = mix_val(nb_em, DM_NB)
                la_mix = mix_val(la_em, DM_LA)

                # elemental ratios using safe division
                mo_ce = _safe_div_val(mo_mix, ce_mix)
                u_th = _safe_div_val(u_mix, th_mix)
                w_th = _safe_div_val(w_mix, th_mix)
                mo_w = _safe_div_val(mo_mix, w_mix)
                nb_la = _safe_div_val(nb_mix, la_mix)

                mix_rows.append({
                    "f": frac,
                    "MO_DM_F_MIX": mo_mix,
                    "CE_DM_F_MIX": ce_mix,
                    "W_DM_F_MIX": w_mix,
                    "U_DM_F_MIX": u_mix,
                    "TH_DM_F_MIX": th_mix,
                    "NB_DM_F_MIX": nb_mix,
                    "LA_DM_F_MIX": la_mix,
                    "MO_CE_MIX": mo_ce,
                    "U_TH_MIX": u_th,
                    "W_TH_MIX": w_th,
                    "MO_W_MIX": mo_w,
                    "NB_LA_MIX": nb_la,
                })

            mix_df = pd.DataFrame(mix_rows)
            all_mixes.append((sal, mix_df))
            sal_label = f" (NaCl wt% = {sal})" if (sal is not None and not pd.isna(sal)) else ""
            print(f"\nTable 5 — DM-F Mixing Table{sal_label} — columns: f, MO_DM_F_MIX, CE_DM_F_MIX, W_DM_F_MIX, U_DM_F_MIX, TH_DM_F_MIX, NB_DM_F_MIX, LA_DM_F_MIX, MO_CE_MIX, U_TH_MIX, W_TH_MIX, MO_W_MIX, NB_LA_MIX")
            print(mix_df.to_string(index=False, float_format=lambda x: f"{x:.8g}"))

        # After building mix tables for all salinities, create the requested plot
        # (no try/except wrapper here — keep errors visible during development)
        try:
            # Build figure
            fig, ax = plt.subplots(figsize=(7, 6))

            # Set logarithmic scales
            ax.set_xscale('log')
            ax.set_yscale('log')

            # Axis limits
            ax.set_ylim(0.01, 1.0)
            ax.set_xlim(0.1, 10.0)

            # Import formatter locally
            from matplotlib.ticker import FuncFormatter

            # Formatter that prints numbers in plain (non-scientific) form
            fmt = FuncFormatter(lambda val, pos: ('%g' % val))
            ax.xaxis.set_major_formatter(fmt)
            ax.yaxis.set_major_formatter(fmt)

            # f values to mark (fractions) and their percentage labels
            label_fs = [0.01, 0.05, 0.1]
            label_percent = {0.01: '1%', 0.05: '5%', 0.1: '10%'}

            # capture axis limits so markers outside the plotting area are not drawn
            x_min, x_max = ax.get_xlim()
            y_min, y_max = ax.get_ylim()

            for sal, mix_df in all_mixes:
                if mix_df is None or mix_df.empty:
                    continue

                # filter finite positive values for log plotting
                valid_mask = np.isfinite(mix_df['U_TH_MIX']) & np.isfinite(mix_df['MO_CE_MIX']) & (mix_df['U_TH_MIX'] > 0) & (mix_df['MO_CE_MIX'] > 0)
                if not np.any(valid_mask):
                    continue

                plot_df = mix_df.loc[valid_mask].copy()

                # sort by x (U/Th) to produce sensible connecting lines
                plot_df.sort_values('U_TH_MIX', inplace=True)
                x = plot_df['U_TH_MIX'].values
                y = plot_df['MO_CE_MIX'].values

                # plot only the connecting line (no markers for every point)
                ax.plot(x, y, color='k', linestyle='-', linewidth=1)

                # plot small 'X' markers for the requested fractions only if they lie within axis limits
                for f in label_fs:
                    sel = plot_df.loc[np.isclose(plot_df['f'].values, f)]
                    if len(sel) == 0:
                        continue
                    try:
                        sx = float(sel['U_TH_MIX'].iat[0])
                        sy = float(sel['MO_CE_MIX'].iat[0])
                    except Exception:
                        continue
                    if not (math.isfinite(sx) and math.isfinite(sy)):
                        continue

                    # skip marker if it falls outside the explicitly set axis limits
                    if not (x_min <= sx <= x_max and y_min <= sy <= y_max):
                        continue

                    # draw a smaller marker (cross)
                    ax.plot(sx, sy, marker='X', color='k', markersize=4, markeredgewidth=0.8)

                # Plot the DM point (f == 0) as a filled dark-red circle with black edge
                try:
                    sel_dm = plot_df.loc[np.isclose(plot_df['f'].values, 0.0)]
                    if len(sel_dm) > 0:
                        sxd = float(sel_dm['U_TH_MIX'].iat[0])
                        syd = float(sel_dm['MO_CE_MIX'].iat[0])
                        if math.isfinite(sxd) and math.isfinite(syd):
                            # only plot DM point if within axis limits
                            if x_min <= sxd <= x_max and y_min <= syd <= y_max:
                                ax.plot(sxd, syd, marker='o', markerfacecolor='#8B0000', markeredgecolor='k', markersize=8, markeredgewidth=1.5)
                except Exception:
                    pass

                # end of per-salinity plotting

            # Final plot adjustments
            ax.set_xlabel("U/Th")
            ax.set_ylabel("Mo/Ce")
            ax.set_title("DM-F Mixing Model Results")
            plt.grid(True, which="both", linestyle='--', linewidth=0.5)
            # legend removed (no labeled series) to avoid the 'Salinity (wt.% NaCl)' title appearing
            # plt.legend(title="Salinity (wt.% NaCl)", loc="best", fontsize=8)
            plt.tight_layout()

            # Show plot
            plt.show()
        except Exception as e:
            print(f"Failed to create plot: {e}")

    # end main

if __name__ == "__main__":
    main()
