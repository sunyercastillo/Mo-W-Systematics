#!/usr/bin/env python3
"""
Bali et al. 2012 Dashboard Interactive Interface

Interactive Streamlit dashboard for exploring the Bali et al. 2012 trace element 
partitioning model. Features:

- Real-time parameter adjustment (Temperature, ΔFMQ, Modal Mineralogy)
- Multi-salinity visualization (all 5 salinities plotted simultaneously)
- Mixing model visualization with custom fluid fractions
- Enhanced dark mode interface with interactive plots
- Fixed pressure at 26,100 bar following Bali et al. 2012

Usage: streamlit run dashboard_interactive_interface.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import math

# Import our modular calculation functions
from CORE_scientific_calculations_engine import (
    BaliConstants, 
    calculate_nacl_properties,
    calculate_molybdenum_partitioning,
    calculate_tungsten_partitioning, 
    calculate_uranium_partitioning,
    calculate_bulk_partitioning,
    calculate_fluid_endmember,
    compute_base_fmq_log10,
    load_excel,
    normalize_headers,
    process_all_vectorized
)

# ========== GLOBAL HELPER FUNCTIONS ==========
def safe_div(num, denom):
    """Safe division that handles zero denominator"""
    return num / denom if denom != 0 else float('nan')

def safe_ratio(num, den):
    """Safe ratio calculation that handles pandas NaN values and zero denominator"""
    if pd.notna(num) and pd.notna(den) and den != 0:
        return num / den
    return np.nan

# ========== PAGE CONFIGURATION ==========
st.set_page_config(
    page_title="Bali et al. 2012 Enhanced Model v5",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Set dark theme with readable text
st.markdown("""
<style>
    .stApp {
        background-color: #0E1117;
        color: #FFFFFF;
    }
    .stSidebar {
        background-color: #262730;
        color: #FFFFFF;
    }
    .stDataFrame {
        background-color: #1E1E1E;
        color: #FFFFFF;
    }
    .stSelectbox > div > div {
        background-color: #262730;
        color: #FFFFFF;
    }
    .stSlider > div > div > div {
        color: #FFFFFF;
    }
    .stCheckbox > label {
        color: #FFFFFF;
    }
    .stMetric {
        background-color: #262730;
        border: 1px solid #404040;
        border-radius: 4px;
        padding: 8px;
        color: #FFFFFF;
    }
    .stMarkdown {
        color: #FFFFFF;
    }
</style>
""", unsafe_allow_html=True)

# ========== CONSTANTS ==========
CONSTANTS = BaliConstants()
FIXED_PRESSURE_BAR = 26100  # Fixed pressure as requested

# Define composition styles with Okabe-Ito colorblind-safe palette (global variable)
composition_styles = {
    "Depleted Mantle (DM)": {'color': '#E69F00', 'symbol': 'star', 'name': 'DM'},  # Orange (Okabe-Ito)
    "Primitive Mantle (PM)": {'color': '#009E73', 'symbol': 'star', 'name': 'PM'}  # Bluish green (Okabe-Ito)
}

# Add Excel composition styles dynamically
def add_excel_composition_style(comp_name, rock_id):
    global composition_styles
    composition_styles[comp_name] = {'color': '#56B4E9', 'symbol': 'star', 'name': f'Excel_{rock_id}'}  # Sky blue (Okabe-Ito)

import os

# ========== EXCEL DATA LOADING ==========
@st.cache_data
def load_input_data():
    """Load and process Excel input data"""
    try:
        # Use relative path to Excel file (works on any computer)
        current_dir = os.path.dirname(os.path.abspath(__file__))
        excel_path = os.path.join(current_dir, "sample_input_data.xlsx")
        
        df_raw = load_excel(excel_path, sheet=0)
        df = normalize_headers(df_raw)
        
        # Add Excel row numbers for display
        df_display = df.copy()
        df_display.insert(0, "Excel_row", range(2, 2 + len(df_display)))
        
        return df, df_display
    except Exception as e:
        st.error(f"Error loading Excel data: {e}")
        return None, None

# Load Excel data
input_df, input_df_display = load_input_data()

# Initialize global variables
excel_mixing_data = None
excel_results_df = None

# ========== DASHBOARD TITLE ==========
st.title("⚛️ Mo-W-Th-U Systematics")
st.markdown("---")

# ========== EXCEL DATA INPUT SECTION ==========
if input_df is not None and input_df_display is not None:
    st.subheader("📊 Source parameters")
    
    # Single-select checkbox per Excel row (data editor with enforced single select)
    row_options = input_df_display["Excel_row"].tolist()
    if 'selected_excel_row' not in st.session_state and row_options:
        st.session_state.selected_excel_row = row_options[0]
    for row_id in row_options:
        key = f"excel_row_checkbox_{row_id}"
        if key not in st.session_state:
            st.session_state[key] = st.session_state.selected_excel_row == row_id

    # Build table with checkbox column first
    selection_flags = [st.session_state[f"excel_row_checkbox_{rid}"] for rid in row_options]
    table_display = input_df_display.copy()
    table_display.insert(0, "Select source", selection_flags)

    prev_selected = st.session_state.get("selected_excel_row")

    edited_table = st.data_editor(
        table_display,
        hide_index=True,
        column_config={"Select source": st.column_config.CheckboxColumn("Select source")},
        disabled=[col for col in table_display.columns if col != "Select source"],
        height=260,
        key="source_table_editor"
    )

    # Enforce single selection without recursive rerun
    selected_rows = edited_table[edited_table["Select source"]]["Excel_row"].tolist()
    
    # Identify the intended selection
    if len(selected_rows) > 0:
        # If multiple selected, pick the one that differs from previous (the new click)
        candidates = [r for r in selected_rows if r != prev_selected]
        new_selection = candidates[0] if candidates else prev_selected
    else:
        # If nothing selected, revert to previous (enforce at least one)
        new_selection = prev_selected

    # Determine if we need to update state or force a visual refresh
    selection_changed = (new_selection != prev_selected)
    visual_mismatch = (len(selected_rows) != 1) or (len(selected_rows)==1 and selected_rows[0] != new_selection)
    
    if selection_changed or visual_mismatch:
        st.session_state.selected_excel_row = new_selection
        for rid in row_options:
            st.session_state[f"excel_row_checkbox_{rid}"] = (rid == new_selection)
        st.rerun()
    
    st.markdown("---")
    
    # Use selected row data to override parameters
    selected_excel_row = st.session_state.get("selected_excel_row")
    if selected_excel_row:
        selected_data = input_df.iloc[selected_excel_row - 2]  # Convert back to 0-based index
        use_excel_data = True
    else:
        use_excel_data = False
        selected_data = None
else:
    st.warning("⚠️ No Excel input data found. Using manual parameter controls.")
    use_excel_data = False
    selected_data = None

# ========== DATA LOADING AND PARAMETER EXTRACTION ==========
# All parameters are extracted from Excel data (no sidebar)

if use_excel_data and selected_data is not None:
    # Use Excel data
    temperature_c = float(selected_data['T_K']) - 273.15
    pressure_bar = float(selected_data['P_bar'])
    dfmq = float(selected_data.get('dFMQ', 0.0))
    
    # Modal mineralogy from Excel
    mode_cpx = float(selected_data.get('mode_cpx', 0.3))
    mode_grt = float(selected_data.get('mode_grt', 0.7))
    mode_rut = float(selected_data.get('mode_rut', 0.0))
    
    # Create custom composition from Excel data
    custom_composition_name = f"Excel_{selected_data['rock_id']}"
    
    # Extract concentrations from Excel data
    excel_composition = {
        'MO': float(selected_data.get('C0_MO', 0.025)),
        'Ce': float(selected_data.get('C0_Ce', 0.772)),
        'W': float(selected_data.get('C0_W', 0.0024)),
        'U': float(selected_data.get('C0_U', 0.0047)),
        'Th': float(selected_data.get('C0_Th', 0.0137)),
        'Nb': float(selected_data.get('C0_Nb', 0.21)),
        'La': float(selected_data.get('C0_La', 0.234))
    }
    
    # Add to constants for calculations
    CONSTANTS.mantle_compositions[custom_composition_name] = excel_composition
    
    # Add to composition styles for plotting
    add_excel_composition_style(custom_composition_name, selected_data['rock_id'])
    
    # Use Excel composition
    selected_compositions = [custom_composition_name]
        
else:
    # Default values when no Excel data
    temperature_c = 1000
    pressure_bar = FIXED_PRESSURE_BAR
    dfmq = 0.0
    
    # Default modal mineralogy
    mode_cpx = 0.3
    mode_grt = 0.7
    mode_rut = 0.0
    
    # Default to Depleted Mantle
    selected_compositions = ["Depleted Mantle (DM)"]

# Convert temperature
temperature_k = temperature_c + 273.15

# Normalize modal proportions
total_modes = mode_cpx + mode_grt + mode_rut
if total_modes > 0:
    mode_cpx_norm = mode_cpx / total_modes
    mode_grt_norm = mode_grt / total_modes  
    mode_rut_norm = mode_rut / total_modes
else:
    mode_cpx_norm = mode_grt_norm = mode_rut_norm = 0

# Fixed fluid fraction steps as specified
fluid_fractions = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
highlight_fraction = 0.01  # 1% to be highlighted

# ========== SALINITY SELECTION (define early for calculations, display later) ==========
# Initialize defaults only once (these will be managed by checkbox keys)
if 'sal_0001' not in st.session_state:
    st.session_state.sal_0001 = True
if 'sal_5' not in st.session_state:
    st.session_state.sal_5 = True
if 'sal_10' not in st.session_state:
    st.session_state.sal_10 = True
if 'sal_15' not in st.session_state:
    st.session_state.sal_15 = True
if 'sal_20' not in st.session_state:
    st.session_state.sal_20 = True

# Configure selected salinities based on current checkbox states
selected_salinities = []
if st.session_state.sal_0001:
    selected_salinities.append(0.001)
if st.session_state.sal_5:
    selected_salinities.append(5.0)
if st.session_state.sal_10:
    selected_salinities.append(10.0)
if st.session_state.sal_15:
    selected_salinities.append(15.0)
if st.session_state.sal_20:
    selected_salinities.append(20.0)

# If no salinity selected, default to 0.001% and 5%
if not selected_salinities:
    selected_salinities = [0.001, 5.0]

# ========== CALCULATE SELECTED SALINITIES FOR ALL COMPOSITIONS ==========
@st.cache_data
def calculate_selected_compositions_salinities(temp_k, press_bar, dfmq_val, cpx_norm, grt_norm, rut_norm, 
                                             composition_names, salinity_list):
    """Calculate results for selected salinity values and all selected mantle compositions - cached for performance."""
    
    all_results = {}
    
    # Base calculations (same for all salinities and compositions)
    log10_fo2_fmq_base = compute_base_fmq_log10(temp_k, press_bar)
    log10_fo2_abs_base = log10_fo2_fmq_base + dfmq_val
    
    modal_props_base = {'cpx': cpx_norm, 'grt': grt_norm, 'rut': rut_norm}
    
    for comp_name in composition_names:
        comp_data = CONSTANTS.mantle_compositions[comp_name]
        composition_results = []
        
        # Extract initial concentrations for this composition
        c0_mo_val = comp_data['MO']
        c0_ce_val = comp_data['Ce']
        c0_w_val = comp_data['W']
        c0_u_val = comp_data['U']
        c0_th_val = comp_data['Th']
        c0_nb_val = comp_data['Nb']
        c0_la_val = comp_data['La']
        
        for sal_wt in salinity_list:
            # NaCl properties
            nacl_props = calculate_nacl_properties(sal_wt)
            nacl_mol = nacl_props['molality']
            log10_nacl = nacl_props['log10_nacl']
            
            # Element partitioning
            mo_res = calculate_molybdenum_partitioning(log10_fo2_abs_base, log10_nacl, temp_k)
            w_res = calculate_tungsten_partitioning(log10_fo2_abs_base, temp_k)
            u_res = calculate_uranium_partitioning(log10_fo2_abs_base, nacl_mol)
            
            # Bulk partition coefficients
            mo_bulk = calculate_bulk_partitioning(
                {'CPX': mo_res['D_CPX'], 'GRT': mo_res['D_GRT'], 'RUT': mo_res['D_RUT']},
                modal_props_base
            )
            w_bulk = calculate_bulk_partitioning(
                {'CPX': w_res['D_CPX'], 'GRT': w_res['D_GRT'], 'RUT': CONSTANTS.partition_constants['W']['RUT']},
                modal_props_base
            )
            u_bulk = calculate_bulk_partitioning(
                {'CPX': u_res['D_CPX'], 'GRT': u_res['D_GRT'], 'RUT': u_res['D_RUT']},
                modal_props_base
            )
            ce_bulk = calculate_bulk_partitioning(CONSTANTS.partition_constants['CE'], modal_props_base)
            th_bulk = calculate_bulk_partitioning(CONSTANTS.partition_constants['TH'], modal_props_base)
            nb_bulk = calculate_bulk_partitioning(CONSTANTS.partition_constants['NB'], modal_props_base)
            la_bulk = calculate_bulk_partitioning(CONSTANTS.partition_constants['LA'], modal_props_base)
            
            # Fluid endmembers
            mo_fluid = calculate_fluid_endmember(c0_mo_val, mo_bulk)
            ce_fluid = calculate_fluid_endmember(c0_ce_val, ce_bulk)
            w_fluid = calculate_fluid_endmember(c0_w_val, w_bulk)
            u_fluid = calculate_fluid_endmember(c0_u_val, u_bulk)
            th_fluid = calculate_fluid_endmember(c0_th_val, th_bulk)
            nb_fluid = calculate_fluid_endmember(c0_nb_val, nb_bulk)
            la_fluid = calculate_fluid_endmember(c0_la_val, la_bulk)
            
            # Calculate ratios
            mo_ce_ratio = mo_fluid / ce_fluid if (np.isfinite(mo_fluid) and np.isfinite(ce_fluid) and ce_fluid != 0) else np.nan
            u_th_ratio = u_fluid / th_fluid if (np.isfinite(u_fluid) and np.isfinite(th_fluid) and th_fluid != 0) else np.nan
            w_th_ratio = w_fluid / th_fluid if (np.isfinite(w_fluid) and np.isfinite(th_fluid) and th_fluid != 0) else np.nan
            nb_la_ratio = nb_fluid / la_fluid if (np.isfinite(nb_fluid) and np.isfinite(la_fluid) and la_fluid != 0) else np.nan
            mo_w_ratio = mo_fluid / w_fluid if (np.isfinite(mo_fluid) and np.isfinite(w_fluid) and w_fluid != 0) else np.nan
            
            composition_results.append({
                'salinity_wt': sal_wt,
                'nacl_molality': nacl_mol,
                'log10_nacl': log10_nacl,
                'mo_ce_ratio': mo_ce_ratio,
                'u_th_ratio': u_th_ratio,
                'w_th_ratio': w_th_ratio,
                'nb_la_ratio': nb_la_ratio,
                'mo_w_ratio': mo_w_ratio,
                'mo_fluid': mo_fluid,
                'ce_fluid': ce_fluid,
                'w_fluid': w_fluid,
                'u_fluid': u_fluid,
                'th_fluid': th_fluid,
                'nb_fluid': nb_fluid,
                'la_fluid': la_fluid,
                'mo_results': mo_res,
                'w_results': w_res,
                'u_results': u_res,
                'composition': comp_name,
                'initial_concentrations': comp_data
            })
        
        all_results[comp_name] = composition_results
    
    return all_results

# Calculate results for all selected compositions and salinities
if not (use_excel_data and excel_mixing_data is not None):
    # Use original calculation method when not using Excel data
    all_composition_results = calculate_selected_compositions_salinities(
        temperature_k, pressure_bar, dfmq, mode_cpx_norm, mode_grt_norm, mode_rut_norm,
        selected_compositions, selected_salinities
    )
else:
    # When using Excel data, create results structure compatible with existing plotting code
    # Use the excel_results_df calculated values instead of mixing data
    all_composition_results = {}
    
    # Create results structure compatible with existing plotting code using calculated fluid endmembers
    for comp_name in selected_compositions:
        comp_results = []
        
        for salinity in selected_salinities:
            # Get calculated fluid endmember values from excel_results_df for this salinity
            if excel_results_df is not None:
                salinity_row = excel_results_df[excel_results_df['NaCl_wt_pct'] == salinity]
                
                if not salinity_row.empty:
                    row = salinity_row.iloc[0]
                    
                    # Calculate ratios from fluid endmembers (these are the calculated values from tables)
                    def safe_ratio(num, den):
                        if pd.notna(num) and pd.notna(den) and den != 0:
                            return num / den
                        return np.nan
                    
                    # Use calculated fluid endmember values to compute ratios
                    mo_ce_ratio = safe_ratio(row['MO_F_EM'], row['CE_F_EM'])
                    u_th_ratio = safe_ratio(row['U_F_EM'], row['TH_F_EM'])
                    w_th_ratio = safe_ratio(row['W_F_EM'], row['TH_F_EM'])
                    nb_la_ratio = safe_ratio(row['NB_F_EM'], row['LA_F_EM'])
                    mo_w_ratio = safe_ratio(row['MO_F_EM'], row['W_F_EM'])
                    
                    # Create result structure compatible with plotting using calculated values
                    comp_results.append({
                        'salinity_wt': salinity,
                        'mo_ce_ratio': mo_ce_ratio,
                        'u_th_ratio': u_th_ratio,
                        'w_th_ratio': w_th_ratio,
                        'nb_la_ratio': nb_la_ratio,
                        'mo_w_ratio': mo_w_ratio,
                        'mo_fluid': row['MO_F_EM'],
                        'ce_fluid': row['CE_F_EM'],
                        'w_fluid': row['W_F_EM'],
                        'u_fluid': row['U_F_EM'],
                        'th_fluid': row['TH_F_EM'],
                        'nb_fluid': row['NB_F_EM'],
                        'la_fluid': row['LA_F_EM'],
                        'composition': comp_name,
                        'initial_concentrations': CONSTANTS.mantle_compositions[comp_name]
                    })
        
        if comp_results:
            all_composition_results[comp_name] = comp_results

# ========== MAIN DASHBOARD LAYOUT ==========
# Unified model conditions display
st.subheader("🌡️ Model Conditions Summary")

# Create comprehensive conditions box
log10_fo2_abs = compute_base_fmq_log10(temperature_k, pressure_bar) + dfmq

# Create 4 columns for the parameters
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("🌡️ Temperature", f"{temperature_c:.2f}°C", f"{temperature_k:.0f} K")

with col2:
    st.metric("⚡ Pressure", f"{pressure_bar/1000:.1f} kbar", f"{pressure_bar:,} bar")

with col3:
    st.metric("🔥 ΔFMQ", f"{dfmq:.2f}", f"log₁₀(fO₂) = {log10_fo2_abs:.2f}")

with col4:
    st.metric("⚙️ Modal Mineralogy", f"CPX: {mode_cpx_norm:.2f}", f"GRT: {mode_grt_norm:.2f} | RUT: {mode_rut_norm:.2f}")

# Salinity Selection
st.subheader("🧂 Salinity Selection")

col_s1, col_s2, col_s3, col_s4, col_s5 = st.columns(5)

with col_s1:
    st.checkbox("0.001% NaCl", key="sal_0001")
with col_s2:
    st.checkbox("5% NaCl", key="sal_5")
with col_s3:
    st.checkbox("10% NaCl", key="sal_10")
with col_s4:
    st.checkbox("15% NaCl", key="sal_15")
with col_s5:
    st.checkbox("20% NaCl", key="sal_20")

st.markdown("---")

# ========== DATA PROCESSING (Calculate first, display later) ==========
if use_excel_data and selected_data is not None:
    # Process selected row with all salinities using the core code's vectorized function
    selected_row_df = input_df.iloc[[selected_excel_row - 2]]  # Convert to 0-based index, keep as DataFrame
    all_salinities = [0.001, 5.0, 10.0, 15.0, 20.0]
    
    # Use the core code's process_all_vectorized function
    results_df = process_all_vectorized(selected_row_df, all_salinities)
    
    # Store results globally for plotting
    excel_results_df = results_df
    
    if results_df is not None and not results_df.empty:
        # Calculate DM ratios from CONSTANTS.dm_concentrations (following core script exactly)
        DM = CONSTANTS.dm_concentrations
        
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
        
        # Generate mixing data for all salinities
        f_values = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.9]
        
        # Function for safe mixing calculation (following core script exactly)
        def mix_val(f_em, dm_const, frac):
            try:
                if not math.isfinite(f_em):
                    return float('nan')
                return (float(f_em) * frac) + (float(dm_const) * (1.0 - frac))
            except Exception:
                return float('nan')
        
        def _safe_div_val(a, b):
            try:
                if not (math.isfinite(a) and math.isfinite(b)):
                    return float('nan')
                if float(b) == 0.0:
                    return float('nan')
                return float(a) / float(b)
            except Exception:
                return float('nan')
        
        # Generate mixing data for each salinity
        all_mixing_data = []
        
        for _, result_row in results_df.iterrows():
            salinity = result_row['NaCl_wt_pct']
            
            # Get fluid endmembers
            mo_em = result_row.get('MO_F_EM', float('nan'))
            ce_em = result_row.get('CE_F_EM', float('nan'))
            w_em = result_row.get('W_F_EM', float('nan'))
            u_em = result_row.get('U_F_EM', float('nan'))
            th_em = result_row.get('TH_F_EM', float('nan'))
            nb_em = result_row.get('NB_F_EM', float('nan'))
            la_em = result_row.get('LA_F_EM', float('nan'))
            
            # Create mixing table for this salinity
            for frac in f_values:
                # Calculate mixing using DM constants (exactly like core script)
                mo_mix = mix_val(mo_em, DM_MO, frac)
                ce_mix = mix_val(ce_em, DM_CE, frac)
                w_mix = mix_val(w_em, DM_W, frac)
                u_mix = mix_val(u_em, DM_U, frac)
                th_mix = mix_val(th_em, DM_TH, frac)
                nb_mix = mix_val(nb_em, DM_NB, frac)
                la_mix = mix_val(la_em, DM_LA, frac)
                
                # Calculate elemental ratios using safe division (exactly like core script)
                mo_ce = _safe_div_val(mo_mix, ce_mix)
                u_th = _safe_div_val(u_mix, th_mix)
                w_th = _safe_div_val(w_mix, th_mix)
                mo_w = _safe_div_val(mo_mix, w_mix)
                nb_la = _safe_div_val(nb_mix, la_mix)
                
                # Add to combined data for plotting
                all_mixing_data.append({
                    'NaCl_wt_pct': salinity,
                    'f': frac,
                    'MO_MIX': mo_mix,
                    'CE_MIX': ce_mix,
                    'W_MIX': w_mix,
                    'U_MIX': u_mix,
                    'TH_MIX': th_mix,
                    'NB_MIX': nb_mix,
                    'LA_MIX': la_mix,
                    'MO_CE_MIX': mo_ce,
                    'U_TH_MIX': u_th,
                    'W_TH_MIX': w_th,
                    'MO_W_MIX': mo_w,
                    'NB_LA_MIX': nb_la,
                })
        
        # Store mixing data for plotting
        excel_mixing_data = pd.DataFrame(all_mixing_data)
    else:
        excel_mixing_data = None
else:
    excel_mixing_data = None

# ========== ENHANCED BALI MODEL PLOTS ==========
# Create four separate plots showing mixing table data
plot_configs = [
    {
        'x_ratio': 'u_th', 'y_ratio': 'mo_ce', 
        'title': 'Mo/Ce vs U/Th', 'x_label': 'U/Th', 'y_label': 'Mo/Ce',
        'x_range': [0.1, 10], 'y_range': [0.01, 10]
    },
    {
        'x_ratio': 'u_th', 'y_ratio': 'w_th', 
        'title': 'W/Th vs U/Th', 'x_label': 'U/Th', 'y_label': 'W/Th',
        'x_range': [0.1, 10], 'y_range': [0.01, 10]
    },
    {
        'x_ratio': 'u_th', 'y_ratio': 'mo_w', 
        'title': 'Mo/W vs U/Th', 'x_label': 'U/Th', 'y_label': 'Mo/W',
        'x_range': [0.1, 10], 'y_range': [0.1, 100]
    },
    {
        'x_ratio': 'nb_la', 'y_ratio': 'mo_w', 
        'title': 'Mo/W vs Nb/La', 'x_label': 'Nb/La', 'y_label': 'Mo/W',
        'x_range': [0.1, 10], 'y_range': [0.1, 100]
    }
]

# Okabe-Ito colorblind-safe palette for different salinities
salinity_color_map = {
    0.001: '#0173B2',  # Blue
    5.0: '#DE8F05',    # Orange
    10.0: '#029E73',   # Bluish green
    15.0: '#CC78BC',   # Reddish purple
    20.0: '#CA3542'    # Vermillion
}

salinity_name_map = {
    0.001: '0.001%',
    5.0: '5%',
    10.0: '10%',
    15.0: '15%',
    20.0: '20%'
}

# Download Configuration UI
st.markdown("### 📈 Interactive Plots")
dl_col1, dl_col2 = st.columns([3, 1])
with dl_col1:
    pass
with dl_col2:
    download_format = st.selectbox(
        "Download format:",
        options=["PNG (300 dpi)", "SVG (Vector)"],
        index=0,
        help="Select the format for the camera download button on the plots."
    )

# Configure Plotly download options
# Display size: 4in x 3in -> 384px x 288px (at 96 DPI)
PLOT_WIDTH = 384
PLOT_HEIGHT = 288

if "PNG" in download_format:
    dl_config = {
        'toImageButtonOptions': {
            'format': 'png',
            'filename': 'bali_model_plot',
            'height': PLOT_HEIGHT,
            'width': PLOT_WIDTH,
            'scale': 4 
        },
        'displayModeBar': True,
        'displaylogo': False
    }
else:
    dl_config = {
        'toImageButtonOptions': {
            'format': 'svg',
            'filename': 'bali_model_plot',
            'height': PLOT_HEIGHT,
            'width': PLOT_WIDTH,
            'scale': 1
        },
        'displayModeBar': True,
        'displaylogo': False
    }

# Create plots in 2x2 grid layout
for row_num in range(2):
    cols = st.columns(2)
    
    for col_num in range(2):
        current_plot_index = row_num * 2 + col_num
        
        if current_plot_index < len(plot_configs):
            config = plot_configs[current_plot_index]
            
            with cols[col_num]:
                fig = go.Figure()
                
                # Plot mixing table data: 5 lines (one per salinity), each with 8 fraction points
                if use_excel_data and excel_mixing_data is not None:
                    for salinity in selected_salinities:
                        # Get mixing data for this salinity
                        salinity_data = excel_mixing_data[excel_mixing_data['NaCl_wt_pct'] == salinity]
                        
                        if not salinity_data.empty:
                            # Sort by fraction f to ensure line is in correct order
                            salinity_data = salinity_data.sort_values('f')
                            
                            # Extract x and y values based on plot configuration
                            if config['x_ratio'] == 'u_th':
                                x_values = salinity_data['U_TH_MIX'].values
                            elif config['x_ratio'] == 'nb_la':
                                x_values = salinity_data['NB_LA_MIX'].values
                            else:
                                continue
                            
                            if config['y_ratio'] == 'mo_ce':
                                y_values = salinity_data['MO_CE_MIX'].values
                            elif config['y_ratio'] == 'w_th':
                                y_values = salinity_data['W_TH_MIX'].values
                            elif config['y_ratio'] == 'mo_w':
                                y_values = salinity_data['MO_W_MIX'].values
                            else:
                                continue
                            
                            # Filter out NaN values
                            valid_mask = np.isfinite(x_values) & np.isfinite(y_values)
                            if np.any(valid_mask):
                                x_clean = x_values[valid_mask]
                                y_clean = y_values[valid_mask]
                                f_values = salinity_data['f'].values[valid_mask]
                                
                                salinity_color = salinity_color_map[salinity]
                                salinity_name = salinity_name_map[salinity]
                                
                                # Add line connecting all 8 fraction points for this salinity
                                fig.add_trace(go.Scatter(
                                    x=x_clean,
                                    y=y_clean,
                                    mode='lines+markers',
                                    name=f'{salinity_name} NaCl',
                                    line=dict(color=salinity_color, width=1.0),
                                    marker=dict(size=5, color=salinity_color,
                                              line=dict(color='black', width=0.3)),
                                    showlegend=True,
                                    hovertemplate=f'<b>{salinity_name} NaCl</b><br>' +
                                     f'f = %{{text}}<br>' +
                                     f'{config["x_label"]}: %{{x:.3f}}<br>' +
                                     f'{config["y_label"]}: %{{y:.3f}}<br>' +
                                     '<extra></extra>',
                                    text=[f'{f:.2f}' for f in f_values]
                                ))
                                
                                # Highlight the 1% mixing point (f=0.01) if it exists
                                # Find index where f is closest to 0.01 (allowing for floating point tolerance)
                                idx_1pct = np.where(np.isclose(f_values, 0.01, atol=1e-4))[0]
                                if len(idx_1pct) > 0:
                                    idx = idx_1pct[0]
                                    fig.add_trace(go.Scatter(
                                        x=[x_clean[idx]],
                                        y=[y_clean[idx]],
                                        mode='text',
                                        text=['1%'],
                                        textposition="top center",
                                        textfont=dict(color=salinity_color, size=10, family="Arial, sans-serif", weight="bold"),
                                        showlegend=False,
                                        hoverinfo='skip'
                                    ))
                
                # Layout configuration
                fig.update_layout(
                    title=dict(
                        text=f"{config['title']}",
                        font=dict(color='#FFFFFF', size=15, family='Arial, sans-serif', weight='bold'),
                        x=0.5,
                        xanchor='center',
                        yanchor='top'
                    ),
                    xaxis_title=dict(text=config['x_label'], font=dict(color='#FFFFFF', size=13, family='Arial, sans-serif', weight='bold')),
                    yaxis_title=dict(text=config['y_label'], font=dict(color='#FFFFFF', size=13, family='Arial, sans-serif', weight='bold')),
                    xaxis_type="log",
                    yaxis_type="log",
                    width=384,   # 4 inches (at 96 dpi)
                    height=288,  # 3 inches (at 96 dpi)
                    showlegend=True,
                    legend=dict(
                        font=dict(color='#FFFFFF', size=8),
                        bgcolor='rgba(38, 39, 48, 0.8)',
                        bordercolor='#FFFFFF',
                        borderwidth=1
                    ),
                    xaxis=dict(
                        type="log",
                        range=[np.log10(config['x_range'][0]), np.log10(config['x_range'][1])],
                        showgrid=True,
                        gridwidth=0.5,
                        gridcolor='rgba(128,128,128,0.4)',
                        griddash='dash',
                        zeroline=False,
                        tickfont=dict(color='#FFFFFF', family='Arial, sans-serif', size=11),
                        title=dict(font=dict(color='#FFFFFF', family='Arial, sans-serif')),
                        mirror=True,
                        linecolor='#000000',
                        linewidth=1,
                        tickvals=[0.1, 1, 10],
                        ticktext=['0.1', '1', '10']
                    ),
                    yaxis=dict(
                        type="log",
                        range=[np.log10(config['y_range'][0]), np.log10(config['y_range'][1])],
                        showgrid=True,
                        gridwidth=0.5,
                        gridcolor='rgba(128,128,128,0.4)',
                        griddash='dash',
                        zeroline=False,
                        tickfont=dict(color='#FFFFFF', family='Arial, sans-serif', size=11),
                        title=dict(font=dict(color='#FFFFFF', family='Arial, sans-serif')),
                        mirror=True,
                        linecolor='#000000',
                        linewidth=1,
                        tickvals=[config['y_range'][0], 1] + ([10] if config['y_range'][1] >= 10 else []) + ([100] if config['y_range'][1] >= 100 else []),
                        ticktext=[str(config['y_range'][0]), '1'] + (['10'] if config['y_range'][1] >= 10 else []) + (['100'] if config['y_range'][1] >= 100 else [])
                    ),
                    plot_bgcolor='white',
                    paper_bgcolor='#262730',
                    font=dict(color='#FFFFFF', family='Arial, sans-serif'),
                    margin=dict(t=60, b=40, l=60, r=20)
                )
                
                # Display the plot
                st.plotly_chart(fig, config=dl_config)

st.markdown("---")

# ========== CALCULATION TABLES (Display tables using already-calculated data) ==========
if use_excel_data and selected_data is not None and results_df is not None and not results_df.empty:
    with st.expander("📋 Detailed Calculation Tables", expanded=False):
        st.write("Click to view detailed calculation results for all salinities")
        st.success(f"✅ Processed {len(results_df)} salinity configurations")
        # Table 1 is hidden as inputs are already shown at the top of the dashboard
        # st.markdown("#### Table 1 — Inputs, logs and input element concentrations")
        # table1_cols = [
        #     "rock_id", "P_bar", "T_K", "dFMQ", "NaCl_wt_pct", "NaCl_m", "log10_NaCl",
        #     "log10_fO2_FMQ", "log10_fO2_abs",
        #     "C0_MO", "C0_Ce", "C0_W", "C0_U", "C0_Th", "C0_Nb", "C0_La",
        # ]
        # available_cols1 = [c for c in table1_cols if c in results_df.columns]
        # if available_cols1:
        #     st.dataframe(results_df[available_cols1], width="stretch")
        
        # Table 2: Per-mineral partition coefficients
        st.markdown("#### Table 2 — Per-mineral partition coefficients")
        table2_cols = [
            "NaCl_wt_pct", "MO_D_CPX", "MO_D_GRT", "MO_D_RUT",
            "W_D_CPX", "W_D_GRT", "U_D_CPX", "U_D_GRT", "U_D_RUT",
        ]
        available_cols2 = [c for c in table2_cols if c in results_df.columns]
        if available_cols2:
            st.dataframe(results_df[available_cols2], use_container_width=True)
        else:
            st.warning("Table 2 columns not found in results")
        
        # Table 3: Bulk (modal-weighted) partition coefficients and fluid endmembers
        st.markdown("#### Table 3 — Bulk partition coefficients and fluid endmembers")
        table3_cols = [
            "NaCl_wt_pct", "MO_Dbulk", "CE_Dbulk", "W_Dbulk", "U_Dbulk",
            "TH_Dbulk", "NB_Dbulk", "LA_Dbulk",
            "MO_F_EM", "CE_F_EM", "W_F_EM", "U_F_EM", "TH_F_EM", "NB_F_EM", "LA_F_EM",
        ]
        available_cols3 = [c for c in table3_cols if c in results_df.columns]
        if available_cols3:
            st.dataframe(results_df[available_cols3], use_container_width=True)
        else:
            st.warning("Table 3 columns not found in results")
        
        # Table 4: Endmember ratios (following exact core script structure)
        st.markdown("#### Table 4 — Endmember ratios")
        
        # Calculate DM ratios from CONSTANTS.dm_concentrations (following core script exactly)
        DM = CONSTANTS.dm_concentrations
        
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
        
        DM_MO_CE = safe_div(DM_MO, DM_CE)
        DM_U_TH = safe_div(DM_U, DM_TH) 
        DM_W_TH = safe_div(DM_W, DM_TH)
        DM_MO_W = safe_div(DM_MO, DM_W)
        DM_NB_LA = safe_div(DM_NB, DM_LA)
        
        # Build Table 4 exactly as in core script
        table4_rows = []
        for _, row in results_df.iterrows():
            # Calculate fluid endmember ratios
            F_MO_CE = safe_div(row.get('MO_F_EM', float('nan')), row.get('CE_F_EM', float('nan')))
            F_U_TH = safe_div(row.get('U_F_EM', float('nan')), row.get('TH_F_EM', float('nan')))
            F_W_TH = safe_div(row.get('W_F_EM', float('nan')), row.get('TH_F_EM', float('nan')))
            F_MO_W = safe_div(row.get('MO_F_EM', float('nan')), row.get('W_F_EM', float('nan')))
            F_NB_LA = safe_div(row.get('NB_F_EM', float('nan')), row.get('LA_F_EM', float('nan')))
            
            table4_rows.append({
                "NaCl_wt_pct": row.get("NaCl_wt_pct"),
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
        st.dataframe(table4_df, use_container_width=True)
        
        # Table 5: DM-F Mixing Model Results (use already-calculated mixing data)
        st.markdown("#### Table 5 — DM-F Mixing Model Results")
        
        # Display separate subtable for each salinity using excel_mixing_data
        if excel_mixing_data is not None:
            salinity_counter = 1
            for salinity in [0.001, 5.0, 10.0, 15.0, 20.0]:
                # Filter data for this salinity
                salinity_mix_data = excel_mixing_data[excel_mixing_data['NaCl_wt_pct'] == salinity]
                
                if not salinity_mix_data.empty:
                    # Create display dataframe
                    mix_df = pd.DataFrame({
                        "f": salinity_mix_data['f'],
                        "MO_DM_F_MIX": salinity_mix_data['MO_MIX'],
                        "CE_DM_F_MIX": salinity_mix_data['CE_MIX'],
                        "W_DM_F_MIX": salinity_mix_data['W_MIX'],
                        "U_DM_F_MIX": salinity_mix_data['U_MIX'],
                        "TH_DM_F_MIX": salinity_mix_data['TH_MIX'],
                        "NB_DM_F_MIX": salinity_mix_data['NB_MIX'],
                        "LA_DM_F_MIX": salinity_mix_data['LA_MIX'],
                        "MO_CE_MIX": salinity_mix_data['MO_CE_MIX'],
                        "U_TH_MIX": salinity_mix_data['U_TH_MIX'],
                        "W_TH_MIX": salinity_mix_data['W_TH_MIX'],
                        "MO_W_MIX": salinity_mix_data['MO_W_MIX'],
                        "NB_LA_MIX": salinity_mix_data['NB_LA_MIX'],
                    })
                    
                    sal_label = f"NaCl wt% = {salinity:.3g}"
                    st.markdown(f"**Table 5.{salinity_counter} — DM-F Mixing Table ({sal_label}):**")
                    st.dataframe(mix_df, width="stretch")
                    salinity_counter += 1

# ========== DATA TABLES ==========
# st.subheader("📋 Complete Results Table (All Selected Compositions & Salinities)")

# # Create comprehensive results table
# table_data = []
# for comp_name, results in all_composition_results.items():
#     comp_short = composition_styles[comp_name]['name']
#     for result in results:
#         table_data.append({
#             'Composition': comp_short,
#             'NaCl_wt_%': result['salinity_wt'],
#             'NaCl_molality': f"{result['nacl_molality']:.3f}",
#             'Mo/Ce': f"{result['mo_ce_ratio']:.4f}" if np.isfinite(result['mo_ce_ratio']) else "NaN",
#             'U/Th': f"{result['u_th_ratio']:.4f}" if np.isfinite(result['u_th_ratio']) else "NaN",
#             'W/Th': f"{result['w_th_ratio']:.4f}" if np.isfinite(result['w_th_ratio']) else "NaN",
#             'Nb/La': f"{result['nb_la_ratio']:.4f}" if np.isfinite(result['nb_la_ratio']) else "NaN",
#             'Mo_fluid_ppm': f"{result['mo_fluid']:.6f}" if np.isfinite(result['mo_fluid']) else "NaN",
#             'Ce_fluid_ppm': f"{result['ce_fluid']:.6f}" if np.isfinite(result['ce_fluid']) else "NaN"
#         })

