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
import matplotlib.pyplot as plt
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

# ========== EXCEL DATA LOADING ==========
@st.cache_data
def load_input_data():
    """Load and process Excel input data"""
    try:
        # Try to load the Excel file
        excel_path = "sample_input_data.xlsx"
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
    
    # Display the Excel data table
    st.dataframe(input_df_display, width="stretch", height=200)
    
    # Composition selection
    st.subheader("✅ Source selection")
    col1, col2 = st.columns(2)
    
    with col1:
        with st.container():
            st.markdown("""
            <div style="background-color: #2E2E2E; padding: 15px; border-radius: 8px; margin-bottom: 10px;">
            """, unsafe_allow_html=True)
            st.markdown("**Select Excel Row for Initial Composition:**")
            selected_excel_row = st.radio(
                "Choose Excel row (starting from row 2):",
                options=input_df_display["Excel_row"].tolist(),
                format_func=lambda x: f"Row {x}: {input_df_display[input_df_display['Excel_row']==x]['rock_id'].iloc[0]}",
                key="excel_row_selection"
            )
            st.markdown("</div>", unsafe_allow_html=True)
    
    with col2:
        with st.container():
            st.markdown("""
            <div style="background-color: #2E2E2E; padding: 15px; border-radius: 8px; margin-bottom: 10px;">
            """, unsafe_allow_html=True)
            st.markdown("**Selected Row Details:**")
            if selected_excel_row:
                selected_row_data = input_df_display[input_df_display['Excel_row'] == selected_excel_row].iloc[0]
                st.write(f"**Rock ID:** {selected_row_data['rock_id']}")
                st.write(f"**Temperature:** {selected_row_data['T_K']-273.15:.1f}°C ({selected_row_data['T_K']:.1f} K)")
                st.write(f"**Pressure:** {selected_row_data['P_bar']/1000:.1f} kbar")
                st.write(f"**ΔFMQ:** {selected_row_data.get('dFMQ', 'N/A')}")
                st.write(f"**Modal Composition:** CPX={selected_row_data.get('mode_cpx', 0):.2f}, GRT={selected_row_data.get('mode_grt', 0):.2f}, RUT={selected_row_data.get('mode_rut', 0):.2f}")
            st.markdown("</div>", unsafe_allow_html=True)
    
    st.markdown("---")
    
    # Use selected row data to override parameters
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

# ========== SIDEBAR: INPUT PARAMETERS ==========
st.sidebar.header("📊 Model Parameters")

# Core thermodynamic conditions
st.sidebar.subheader("Thermodynamic Conditions")

if use_excel_data and selected_data is not None:
    # Use Excel data with display
    temperature_c = float(selected_data['T_K']) - 273.15
    pressure_bar = float(selected_data['P_bar'])
    dfmq = float(selected_data.get('dFMQ', 0.0))
    
    st.sidebar.write(f"**Temperature:** {temperature_c:.1f}°C (from Excel)")
    st.sidebar.write(f"**Pressure:** {pressure_bar/1000:.1f} kbar (from Excel)")
    st.sidebar.write(f"**ΔFMQ:** {dfmq:.2f} (from Excel)")
else:
    # Manual parameter controls
    temperature_c = st.sidebar.slider("Temperature (°C)", 400, 1200, 1000, 25)
    pressure_bar = FIXED_PRESSURE_BAR  # Fixed pressure
    dfmq = st.sidebar.slider("ΔFMQ (log units)", -4.0, 4.0, 0.0, 0.1)

# Convert temperature
temperature_k = temperature_c + 273.15

# ========== SALINITY CONFIGURATION ==========
# This will be configured after sidebar checkboxes are defined

# ========== MANTLE COMPOSITION SELECTION ==========
st.sidebar.subheader("🌍 Initial Mantle Composition")

if use_excel_data and selected_data is not None:
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
    
    # Use Excel composition as default
    selected_compositions = [custom_composition_name]
    
    st.sidebar.write(f"**Using Excel Composition:** {selected_data['rock_id']}")
    st.sidebar.write("**Concentrations (ppm):**")
    for element, conc in excel_composition.items():
        st.sidebar.write(f"  {element}: {conc:.6f}")
    
    # Option to add standard compositions
    st.sidebar.markdown("**Additional Standard Compositions:**")
    use_dm = st.sidebar.checkbox("✅ Add Depleted Mantle (DM)", value=False, help="Workman & Hart (2005)")
    use_pm = st.sidebar.checkbox("✅ Add Primitive Mantle (PM)", value=False, help="Palme & O'Neill (2014)")
    
    if use_dm:
        selected_compositions.append("Depleted Mantle (DM)")
    if use_pm:
        selected_compositions.append("Primitive Mantle (PM)")
        
else:
    # Standard composition selection when no Excel data
    use_dm = st.sidebar.checkbox("✅ Depleted Mantle (DM)", value=True, help="Workman & Hart (2005)")
    use_pm = st.sidebar.checkbox("✅ Primitive Mantle (PM)", value=False, help="Palme & O'Neill (2014)")
    
    # Display selected compositions
    selected_compositions = []
    if use_dm:
        selected_compositions.append("Depleted Mantle (DM)")
    if use_pm:
        selected_compositions.append("Primitive Mantle (PM)")
    
    # If no composition selected, default to DM
    if not selected_compositions:
        selected_compositions = ["Depleted Mantle (DM)"]
        use_dm = True
        st.sidebar.warning("⚠️ No composition selected - defaulting to DM")

# ========== SALINITY PARAMETERS ==========
st.sidebar.subheader("🧂 Salinity")

# Checkboxes for controlling mixing line visibility for each salinity
show_mixing_0001 = st.sidebar.checkbox("0.001% Salinity", value=True)
show_mixing_05 = st.sidebar.checkbox("5% Salinity", value=True)
show_mixing_10 = st.sidebar.checkbox("10% Salinity", value=True)
show_mixing_15 = st.sidebar.checkbox("15% Salinity", value=True)
show_mixing_20 = st.sidebar.checkbox("20% Salinity", value=True)

# Create mapping for easy access
mixing_visibility = {
    0.001: show_mixing_0001,
    5.0: show_mixing_05,
    10.0: show_mixing_10,
    15.0: show_mixing_15,
    20.0: show_mixing_20
}

# Configure selected salinities based on checkboxes
selected_salinities = []
if show_mixing_0001:
    selected_salinities.append(0.001)
if show_mixing_05:
    selected_salinities.append(5.0)
if show_mixing_10:
    selected_salinities.append(10.0)
if show_mixing_15:
    selected_salinities.append(15.0)
if show_mixing_20:
    selected_salinities.append(20.0)

# If no salinity selected, default to 0.001% and 5%
if not selected_salinities:
    selected_salinities = [0.001, 5.0]

# Fixed fluid fraction steps as specified
fluid_fractions = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
highlight_fraction = 0.01  # 1% to be highlighted

# Modal mineralogy
st.sidebar.subheader("Modal Mineralogy")

if use_excel_data and selected_data is not None:
    # Use Excel data with display
    mode_cpx = float(selected_data.get('mode_cpx', 0.3))
    mode_grt = float(selected_data.get('mode_grt', 0.7))
    mode_rut = float(selected_data.get('mode_rut', 0.0))
    
    st.sidebar.write(f"**Clinopyroxene:** {mode_cpx:.3f} (from Excel)")
    st.sidebar.write(f"**Garnet:** {mode_grt:.3f} (from Excel)")
    st.sidebar.write(f"**Rutile:** {mode_rut:.3f} (from Excel)")
else:
    # Manual parameter controls
    mode_cpx = st.sidebar.slider("Clinopyroxene", 0.0, 1.0, 0.3, 0.01)
    mode_grt = st.sidebar.slider("Garnet", 0.0, 1.0, 0.7, 0.01)
    mode_rut = st.sidebar.slider("Rutile", 0.0, 1.0, 0.0, 0.01)

# Normalize modal proportions
total_modes = mode_cpx + mode_grt + mode_rut
if total_modes > 0:
    mode_cpx_norm = mode_cpx / total_modes
    mode_grt_norm = mode_grt / total_modes  
    mode_rut_norm = mode_rut / total_modes
else:
    mode_cpx_norm = mode_grt_norm = mode_rut_norm = 0

st.sidebar.markdown(f"**Normalized:** CPX={mode_cpx_norm:.2f}, GRT={mode_grt_norm:.2f}, RUT={mode_rut_norm:.2f}")

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
# Top row: Summary metrics
st.subheader("🌡️ Current Conditions Summary")
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("Temperature", f"{temperature_c}°C", f"{temperature_k:.0f} K")
with col2:
    st.metric("Pressure", f"{pressure_bar/1000:.1f} kbar", f"{pressure_bar:,} bar (FIXED)")
with col3:
    st.metric("ΔFMQ", f"{dfmq:.2f}")
with col4:
    log10_fo2_abs = compute_base_fmq_log10(temperature_k, pressure_bar) + dfmq
    st.metric("log₁₀(fO₂)", f"{log10_fo2_abs:.2f}")

st.markdown("---")

# ========== MODEL CONDITIONS ==========
# Display unified model conditions
st.markdown(f"""
<div style='text-align: center; background-color: #262730; padding: 8px; border-radius: 4px; margin-bottom: 15px; border: 1px solid #404040;'>
    <h5 style='margin: 0; color: #FFFFFF; font-family: Arial, sans-serif; font-weight: bold;'>Model Conditions</h5>
    <p style='margin: 3px 0; color: #CCCCCC; font-size: 13px; font-family: Arial, sans-serif;'>
        <b>Temperature:</b> {temperature_c}°C | <b>Pressure:</b> {pressure_bar/1000:.1f} kbar | <b>ΔFMQ:</b> {dfmq:.1f}<br>
        <b>Modal Composition:</b> CPX={mode_cpx_norm:.2f}, GRT={mode_grt_norm:.2f}, RUT={mode_rut_norm:.2f}
    </p>
</div>
""", unsafe_allow_html=True)

# ========== ENHANCED BALI MODEL PLOTS WITH MIXING ==========
# Create four separate plots as requested with specific axis limits
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

# Store all figures to create a unified legend later
all_figures = []

# Okabe-Ito colorblind-safe palette for different salinities
salinity_color_map = {
    0.001: '#0173B2',  # Blue (Okabe-Ito)
    5.0: '#DE8F05',    # Orange (Okabe-Ito)
    10.0: '#029E73',   # Bluish green (Okabe-Ito)
    15.0: '#CC78BC',   # Reddish purple (Okabe-Ito)
    20.0: '#CA3542'    # Vermillion (Okabe-Ito)
}

salinity_name_map = {
    0.001: '0.001%',
    5.0: '5%',
    10.0: '10%',
    15.0: '15%',
    20.0: '20%'
}

# Create plots in 2x2 grid layout - centered
# Add centering container
st.markdown('<div style="display: flex; justify-content: center;">', unsafe_allow_html=True)

for row_num in range(2):
    cols = st.columns(2)  # Create 2 columns for each row
    
    for col_num in range(2):
        current_plot_index = row_num * 2 + col_num  # Calculate the plot index
        
        # Ensure plot_configs is properly defined and current_plot_index is an integer
        if isinstance(plot_configs, list) and isinstance(current_plot_index, int) and current_plot_index < len(plot_configs):
            config = plot_configs[current_plot_index]
            
            with cols[col_num]:  # Display plot in the appropriate column
                fig = go.Figure()
                
                # Plot mixing lines from Table 5 data for each salinity
                if use_excel_data and excel_mixing_data is not None:
                    for salinity in selected_salinities:
                        if mixing_visibility.get(salinity, True):  # Only plot if checkbox is checked
                            # Get mixing data for this salinity from Table 5
                            salinity_data = excel_mixing_data[excel_mixing_data['NaCl_wt_pct'] == salinity]
                            
                            if not salinity_data.empty:
                                # Extract x and y values based on plot configuration
                                if config['x_ratio'] == 'u_th':
                                    x_values = salinity_data['U_TH_MIX'].values
                                elif config['x_ratio'] == 'nb_la':
                                    x_values = salinity_data['NB_LA_MIX'].values
                                
                                if config['y_ratio'] == 'mo_ce':
                                    y_values = salinity_data['MO_CE_MIX'].values
                                elif config['y_ratio'] == 'w_th':
                                    y_values = salinity_data['W_TH_MIX'].values
                                elif config['y_ratio'] == 'mo_w':
                                    y_values = salinity_data['MO_W_MIX'].values
                                
                                # Filter out NaN values
                                valid_mask = np.isfinite(x_values) & np.isfinite(y_values)
                                if np.any(valid_mask):
                                    x_clean = x_values[valid_mask]
                                    y_clean = y_values[valid_mask]
                                    f_values = salinity_data['f'].values[valid_mask]
                                    
                                    salinity_color = salinity_color_map[salinity]
                                    salinity_name = salinity_name_map[salinity]
                                    
                                    # Add mixing line (solid line connecting all mixing points)
                                    fig.add_trace(go.Scatter(
                                        x=x_clean,
                                        y=y_clean,
                                        mode='lines',
                                        name=f'Mix {salinity_name}',
                                        line=dict(color=salinity_color, width=1.5),
                                        showlegend=True,
                                        hovertemplate=f'<b>{salinity_name} NaCl Mixing</b><br>' +
                                         f'{config["x_label"]}: %{{x:.3f}}<br>' +
                                         f'{config["y_label"]}: %{{y:.3f}}<br>' +
                                         '<extra></extra>',
                                        opacity=0.9
                                    ))
                                    
                                    # Add 1% marker as X
                                    one_percent_mask = np.isclose(f_values, 0.01)
                                    if np.any(one_percent_mask):
                                        one_x = x_clean[one_percent_mask]
                                        one_y = y_clean[one_percent_mask]
                                        fig.add_trace(go.Scatter(
                                            x=one_x,
                                            y=one_y,
                                            mode='markers',
                                            name=f'1% {salinity_name}',
                                            marker=dict(size=5, color=salinity_color, symbol='x',
                                                      line=dict(color='black', width=1)),
                                            showlegend=False,
                                            hovertemplate=f'<b>1% {salinity_name} NaCl</b><br>' +
                                             f'{config["x_label"]}: %{{x:.3f}}<br>' +
                                             f'{config["y_label"]}: %{{y:.3f}}<br>' +
                                             '<extra></extra>'
                                        ))
                                    
                                    # Add 5% marker as X
                                    five_percent_mask = np.isclose(f_values, 0.05)
                                    if np.any(five_percent_mask):
                                        five_x = x_clean[five_percent_mask]
                                        five_y = y_clean[five_percent_mask]
                                        fig.add_trace(go.Scatter(
                                            x=five_x,
                                            y=five_y,
                                            mode='markers',
                                            name=f'5% {salinity_name}',
                                            marker=dict(size=5, color=salinity_color, symbol='x',
                                                      line=dict(color='black', width=1)),
                                            showlegend=False,
                                            hovertemplate=f'<b>5% {salinity_name} NaCl</b><br>' +
                                             f'{config["x_label"]}: %{{x:.3f}}<br>' +
                                             f'{config["y_label"]}: %{{y:.3f}}<br>' +
                                             '<extra></extra>'
                                        ))
                                    
                                    # Add DM point (f=0, mantle source) as star
                                    dm_mask = np.isclose(f_values, 0.0)
                                    if np.any(dm_mask):
                                        dm_x = x_clean[dm_mask]
                                        dm_y = y_clean[dm_mask]
                                        fig.add_trace(go.Scatter(
                                            x=dm_x,
                                            y=dm_y,
                                            mode='markers',
                                            name=f'DM',
                                            marker=dict(size=12, color='white', symbol='star',
                                                      line=dict(color='black', width=1)),
                                            showlegend=(salinity == selected_salinities[0]),  # Show legend only once
                                            hovertemplate=f'<b>DM (Depleted Mantle)</b><br>' +
                                             f'{config["x_label"]}: %{{x:.3f}}<br>' +
                                             f'{config["y_label"]}: %{{y:.3f}}<br>' +
                                             '<extra></extra>'
                                        ))
                
                # Enhanced layout with light gray plot background, logarithmic axes, and specific ranges
                            else:
                                # Original approach using fluid_fractions
                                if highlight_fraction in fluid_fractions:
                                    idx_1_percent = fluid_fractions.index(highlight_fraction)
                                    if idx_1_percent < len(mixing_x):
                                        fig.add_trace(go.Scatter(
                                            x=[mixing_x[idx_1_percent]],
                                            y=[mixing_y[idx_1_percent]],
                                            mode='markers+text',
                                            name=f'{comp_style["name"]} 1% {salinity_name}',
                                            marker=dict(size=5, color=salinity_color, symbol='x',
                                                      line=dict(color='black', width=1)),
                                            text=['1%'],
                                            textposition='top center',
                                            textfont=dict(color='black', size=9, family='Arial, sans-serif'),
                                            showlegend=False,
                                            hovertemplate=f'<b>{comp_style["name"]} - {salinity_name} NaCl</b><br>' +
                                             f'1% Fluid Fraction<br>' +
                                             f'{config["x_label"]}: %{{x:.3f}}<br>' +
                                             f'{config["y_label"]}: %{{y:.3f}}<br>' +
                                             '<extra></extra>'
                                        ))
                            
                            # Add 5% marker as x (use Excel data if available, otherwise use mixing_x/mixing_y)
                            if use_excel_data and excel_mixing_data is not None:
                                # Find 5% point in Excel mixing data
                                five_percent_data = excel_mixing_data[
                                    (excel_mixing_data['NaCl_wt_pct'] == salinity) & 
                                    (excel_mixing_data['f'] == 0.05)
                                ]
                                if not five_percent_data.empty:
                                    row = five_percent_data.iloc[0]
                                    if config['x_ratio'] == 'u_th':
                                        five_x = row['U_TH_MIX']
                                    elif config['x_ratio'] == 'nb_la':
                                        five_x = row['NB_LA_MIX']
                                    else:
                                        five_x = np.nan
                                        
                                    if config['y_ratio'] == 'mo_ce':
                                        five_y = row['MO_CE_MIX']
                                    elif config['y_ratio'] == 'w_th':
                                        five_y = row['W_TH_MIX']
                                    elif config['y_ratio'] == 'mo_w':
                                        five_y = row['MO_W_MIX']
                                    else:
                                        five_y = np.nan
                                    
                                    if np.isfinite(five_x) and np.isfinite(five_y):
                                        fig.add_trace(go.Scatter(
                                            x=[five_x],
                                            y=[five_y],
                                            mode='markers+text',
                                            name=f'{comp_style["name"]} 5% {salinity_name}',
                                            marker=dict(size=5, color=salinity_color, symbol='x',
                                                      line=dict(color='black', width=1)),
                                            text=['5%'],
                                            textposition='top center',
                                            textfont=dict(color='black', size=9, family='Arial, sans-serif'),
                                            showlegend=False,
                                            hovertemplate=f'<b>{comp_style["name"]} - {salinity_name} NaCl</b><br>' +
                                             f'5% Fluid Fraction (Excel)<br>' +
                                             f'{config["x_label"]}: %{{x:.3f}}<br>' +
                                             f'{config["y_label"]}: %{{y:.3f}}<br>' +
                                             '<extra></extra>'
                                        ))
                            else:
                                # Original approach using fluid_fractions
                                if 0.05 in fluid_fractions:
                                    idx_5_percent = fluid_fractions.index(0.05)
                                    if idx_5_percent < len(mixing_x):
                                        fig.add_trace(go.Scatter(
                                            x=[mixing_x[idx_5_percent]],
                                            y=[mixing_y[idx_5_percent]],
                                            mode='markers+text',
                                            name=f'{comp_style["name"]} 5% {salinity_name}',
                                            marker=dict(size=5, color=salinity_color, symbol='x',
                                                      line=dict(color='black', width=1)),
                                            text=['5%'],
                                            textposition='top center',
                                            textfont=dict(color='black', size=9, family='Arial, sans-serif'),
                                            showlegend=False,
                                            hovertemplate=f'<b>{comp_style["name"]} - {salinity_name} NaCl</b><br>' +
                                             f'5% Fluid Fraction<br>' +
                                             f'{config["x_label"]}: %{{x:.3f}}<br>' +
                                             f'{config["y_label"]}: %{{y:.3f}}<br>' +
                                             '<extra></extra>'
                                        ))

                # Add initial composition markers (brought to front)
                for comp_name in selected_compositions:
                    comp_style = composition_styles[comp_name]
                    comp_data = CONSTANTS.mantle_compositions[comp_name]
                    
                    # Calculate initial composition ratios
                    if config['x_ratio'] == 'u_th':
                        initial_x = comp_data['U'] / comp_data['Th']
                    elif config['x_ratio'] == 'nb_la':
                        initial_x = comp_data['Nb'] / comp_data['La']
                    else:
                        continue  # Skip if ratio not recognized
                        
                    if config['y_ratio'] == 'mo_ce':
                        initial_y = comp_data['MO'] / comp_data['Ce']
                    elif config['y_ratio'] == 'w_th':
                        initial_y = comp_data['W'] / comp_data['Th']
                    elif config['y_ratio'] == 'mo_w':
                        initial_y = comp_data['MO'] / comp_data['W']
                    else:
                        continue  # Skip if ratio not recognized
                    
                    # Add initial composition marker
                    fig.add_trace(go.Scatter(
                        x=[initial_x],
                        y=[initial_y],
                        mode='markers',
                        name=f'{comp_style["name"]} Source',
                        marker=dict(size=12, color=comp_style['color'], symbol='star',
                                    line=dict(color='black', width=1)),
                        showlegend=True,
                        hovertemplate=f'<b>{comp_style["name"]} Initial Composition</b><br>' +
                         f'{config["x_label"]}: %{{x:.3f}}<br>' +
                         f'{config["y_label"]}: %{{y:.3f}}<br>' +
                         '<extra></extra>'
                    ))

                # Enhanced layout with light gray plot background, logarithmic axes, and specific ranges
                fig.update_layout(
                    title=dict(
                        text=f"{config['title']}",
                        font=dict(color='#FFFFFF', size=15, family='Arial, sans-serif', weight='bold'),
                        x=0.5
                    ),
                    xaxis_title=dict(text=config['x_label'], font=dict(color='#FFFFFF', size=13, family='Arial, sans-serif', weight='bold')),
                    yaxis_title=dict(text=config['y_label'], font=dict(color='#FFFFFF', size=13, family='Arial, sans-serif', weight='bold')), 
                    xaxis_type="log",
                    yaxis_type="log",
                    height=450,  # Reduced height to make space for unified legend
                    showlegend=False,  # Hide individual legends
                    xaxis=dict(
                        type="log",
                        range=[np.log10(config['x_range'][0]), np.log10(config['x_range'][1])],  # Set specific range
                        showgrid=True,
                        gridwidth=1,
                        gridcolor='rgba(128,128,128,0.4)',
                        griddash='dash',
                        zeroline=False,
                        tickfont=dict(color='#FFFFFF', family='Arial, sans-serif', size=11),  # White text for better contrast
                        title=dict(font=dict(color='#FFFFFF', family='Arial, sans-serif')),
                        mirror=True,
                        linecolor='#000000',  # Black axis lines
                        linewidth=2,
                        tickvals=[0.1, 1, 10],  # Only show these tick values
                        ticktext=['0.1', '1', '10']
                    ),
                    yaxis=dict(
                        type="log",
                        range=[np.log10(config['y_range'][0]), np.log10(config['y_range'][1])],  # Set specific range
                        showgrid=True,
                        gridwidth=1,
                        gridcolor='rgba(128,128,128,0.4)',
                        griddash='dash',
                        zeroline=False,
                        tickfont=dict(color='#FFFFFF', family='Arial, sans-serif', size=11),  # White text for better contrast
                        title=dict(font=dict(color='#FFFFFF', family='Arial, sans-serif')),
                        mirror=True,
                        linecolor='#000000',  # Black axis lines
                        linewidth=2,
                        tickvals=[config['y_range'][0], 1] + ([10] if config['y_range'][1] >= 10 else []) + ([100] if config['y_range'][1] >= 100 else []),
                        ticktext=[str(config['y_range'][0]), '1'] + (['10'] if config['y_range'][1] >= 10 else []) + (['100'] if config['y_range'][1] >= 100 else [])
                    ),
                    plot_bgcolor='white',  # White background inside plot area
                    paper_bgcolor='#262730',  # Match Model Conditions container background
                    font=dict(color='#FFFFFF', family='Arial, sans-serif'),
                    margin=dict(t=60, b=40, l=60, r=20)  # Reduced margins
                )

                # Store figure for unified legend
                all_figures.append(fig)
                
                # Display the plot in this column
                st.plotly_chart(fig, width="stretch")

# Close centering container
st.markdown('</div>', unsafe_allow_html=True)

st.markdown("---")


# ========== LEGEND ==========
st.markdown("### 🗂️ Legend")

# Create legend using HTML for better control - clean white theme
legend_html = f"""
<div style='background-color: #f8f9fa; padding: 15px; border-radius: 8px; border: 2px solid #dee2e6; margin: 10px 0; font-family: Arial, sans-serif; box-shadow: 0 2px 4px rgba(0,0,0,0.1);'>
    <div style='display: flex; flex-wrap: wrap; justify-content: space-around; align-items: center;'>"""

# Add salinity colors dynamically based on selection (Okabe-Ito colorblind-safe palette)
salinity_info = {
    0.001: ('#0173B2', '0.001% NaCl'),
    5.0: ('#DE8F05', '5% NaCl'),
    10.0: ('#029E73', '10% NaCl'),
    15.0: ('#CC78BC', '15% NaCl'),
    20.0: ('#CA3542', '20% NaCl')
}

for salinity in selected_salinities:
    color, label = salinity_info[salinity]
    legend_html += f"""
        <div style='margin: 6px 12px;'>
            <span style='display: inline-block; width: 18px; height: 18px; background-color: {color}; margin-right: 10px; border-radius: 3px; border: 1px solid #333;'></span>
            <strong style='color: #333333; font-size: 14px;'>{label}</strong>
        </div>"""

legend_html += """
    </div>
    <hr style='margin: 10px 0; border: 0; border-top: 2px solid #dee2e6;'>
    <div style='display: flex; flex-wrap: wrap; justify-content: space-around; align-items: center;'>"""

# Add composition symbols dynamically based on selection (Okabe-Ito colors)
if "Depleted Mantle (DM)" in selected_compositions:
    legend_html += """
        <div style='margin: 6px 12px;'>
            <span style='color: #E69F00; font-size: 18px; margin-right: 10px; text-shadow: 1px 1px 1px rgba(0,0,0,0.3);'>★</span>
            <strong style='color: #333333; font-size: 14px;'>Depleted Mantle (DM)</strong>
        </div>"""

if "Primitive Mantle (PM)" in selected_compositions:
    legend_html += """
        <div style='margin: 6px 12px;'>
            <span style='color: #009E73; font-size: 18px; margin-right: 10px; text-shadow: 1px 1px 1px rgba(0,0,0,0.3);'>★</span>
            <strong style='color: #333333; font-size: 14px;'>Primitive Mantle (PM)</strong>
        </div>"""

legend_html += """
        <div style='margin: 6px 12px;'>
            <span style='color: #666666; font-style: italic; font-size: 13px;'>X markers = 1% & 5% fluid fractions</span>
        </div>
    </div>
</div>
"""

st.markdown(legend_html, unsafe_allow_html=True)

st.markdown("---")

# ========== CALCULATION TABLES (like core code) ==========
if use_excel_data and selected_data is not None:
    st.subheader("📋 Detailed Calculation Tables")
    
    # Process selected row with all salinities using the core code's vectorized function
    selected_row_df = input_df.iloc[[selected_excel_row - 2]]  # Convert to 0-based index, keep as DataFrame
    all_salinities = [0.001, 5.0, 10.0, 15.0, 20.0]
    
    # Use the core code's process_all_vectorized function
    results_df = process_all_vectorized(selected_row_df, all_salinities)
    
    # Store results globally for plotting
    excel_results_df = results_df
    
    if results_df is not None and not results_df.empty:
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
            st.dataframe(results_df[available_cols2], width="stretch")
        
        # Table 3: Bulk (modal-weighted) partition coefficients and fluid endmembers
        st.markdown("#### Table 3 — Bulk partition coefficients and fluid endmembers")
        table3_cols = [
            "NaCl_wt_pct", "MO_Dbulk", "CE_Dbulk", "W_Dbulk", "U_Dbulk",
            "TH_Dbulk", "NB_Dbulk", "LA_Dbulk",
            "MO_F_EM", "CE_F_EM", "W_F_EM", "U_F_EM", "TH_F_EM", "NB_F_EM", "LA_F_EM",
        ]
        available_cols3 = [c for c in table3_cols if c in results_df.columns]
        if available_cols3:
            st.dataframe(results_df[available_cols3], width="stretch")
        
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
        st.dataframe(table4_df, width="stretch")
        
        # Table 5: DM-F Mixing Model Results (separate subtable for each salinity - following exact core script)
        st.markdown("#### Table 5 — DM-F Mixing Model Results")
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
        
        # Generate separate subtable for each salinity (exactly like core script)
        salinity_counter = 1
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
            mix_rows = []
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
                
                # Also add to combined data for plotting
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
            
            # Display subtable for this salinity (exactly like core script)
            mix_df = pd.DataFrame(mix_rows)
            sal_label = f"NaCl wt% = {salinity:.3g}" if (salinity is not None and not pd.isna(salinity)) else "Unknown salinity"
            st.markdown(f"**Table 5.{salinity_counter} — DM-F Mixing Table ({sal_label}):**")
            # st.text("Columns: f, MO_DM_F_MIX, CE_DM_F_MIX, W_DM_F_MIX, U_DM_F_MIX, TH_DM_F_MIX, NB_DM_F_MIX, LA_DM_F_MIX, MO_CE_MIX, U_TH_MIX, W_TH_MIX, MO_W_MIX, NB_LA_MIX")
            st.dataframe(mix_df, width="stretch")
            salinity_counter += 1
        
        # Store mixing data for plotting
        excel_mixing_data = pd.DataFrame(all_mixing_data)
    else:
        st.error("Failed to process calculation data")
        excel_mixing_data = None
else:
    excel_mixing_data = None

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

