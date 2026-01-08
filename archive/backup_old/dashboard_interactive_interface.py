#!/usr/bin/env python3
"""
Enhanced Bali et al. 2012 Model Dashboard v5 - Enhanced with Mixing Model

Complete interactive implementation with:
- Fixed mantle compositions (DM, PM)
- Multi-salinity analysis with individual control
- fO2 variation and modal mineralogy exploration
- Fixed pressure at 26100 bar
- Mixing model visualization
- Enhanced plot aesthetics with white grid lines

Usage: streamlit run enhanced_dashboard_v5.py
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
from Bali_v5_dashboard import (
    BaliConstants, 
    calculate_nacl_properties,
    calculate_molybdenum_partitioning,
    calculate_tungsten_partitioning, 
    calculate_uranium_partitioning,
    calculate_bulk_partitioning,
    calculate_fluid_endmember,
    compute_base_fmq_log10
)

# ========== PAGE CONFIGURATION ==========
st.set_page_config(
    page_title="Bali et al. 2012 Enhanced Model v5",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Set dark theme
st.markdown("""
<style>
    .stApp {
        background-color: #0E1117;
        color: white;
    }
    .stSidebar {
        background-color: #262730;
    }
    .stDataFrame {
        background-color: #1E1E1E;
    }
</style>
""", unsafe_allow_html=True)

# ========== CONSTANTS ==========
CONSTANTS = BaliConstants()
FIXED_PRESSURE_BAR = 26100  # Fixed pressure as requested

# ========== DASHBOARD TITLE ==========
st.title("🧪 Bali et al. 2012 Enhanced Interactive Model v5")
st.markdown("**Complete Trace Element Partitioning Analysis with Mixing Model • Fixed Pressure • Selectable Salinities • Enhanced Visualization**")
st.markdown("---")

# ========== SIDEBAR: INPUT PARAMETERS ==========
st.sidebar.header("📊 Model Parameters")

# Core thermodynamic conditions
st.sidebar.subheader("Thermodynamic Conditions")
temperature_c = st.sidebar.slider("Temperature (°C)", 400, 1200, 1000, 25)
pressure_bar = FIXED_PRESSURE_BAR  # Fixed pressure
dfmq = st.sidebar.slider("ΔFMQ (log units)", -4.0, 4.0, 0.0, 0.1)

# Convert temperature
temperature_k = temperature_c + 273.15

# Display fixed pressure value
st.sidebar.info(f"**Fixed Pressure:** {pressure_bar:,} bar ({pressure_bar/1000:.1f} kbar)")

# ========== SALINITY CONFIGURATION ==========
# Always use all 5 salinities - no selection needed
selected_salinities = [0.001, 5.0, 10.0, 15.0, 20.0]

st.sidebar.markdown("### 📊 **Plotting All 5 Salinities**")
st.sidebar.markdown("🔹 Freshwater (0.001% NaCl)")
st.sidebar.markdown("🔸 5% NaCl") 
st.sidebar.markdown("🔹 10% NaCl")
st.sidebar.markdown("🔸 15% NaCl")
st.sidebar.markdown("🔹 20% NaCl")

# ========== MANTLE COMPOSITION SELECTION ==========
st.sidebar.subheader("🌍 Initial Mantle Composition")

# Checkbox system for mantle selection (hidden details as requested)
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

# ========== MIXING MODEL PARAMETERS ==========
st.sidebar.subheader("🔄 Mixing Model")
st.sidebar.markdown("*Control mixing line visibility:*")

# Fixed fluid fraction steps as specified
fluid_fractions = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
highlight_fraction = 0.01  # 1% to be highlighted

# Modal mineralogy
st.sidebar.subheader("Modal Mineralogy")
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
            
            composition_results.append({
                'salinity_wt': sal_wt,
                'nacl_molality': nacl_mol,
                'log10_nacl': log10_nacl,
                'mo_ce_ratio': mo_ce_ratio,
                'u_th_ratio': u_th_ratio,
                'w_th_ratio': w_th_ratio,
                'nb_la_ratio': nb_la_ratio,
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
all_composition_results = calculate_selected_compositions_salinities(
    temperature_k, pressure_bar, dfmq, mode_cpx_norm, mode_grt_norm, mode_rut_norm,
    selected_compositions, selected_salinities
)

# ========== MAIN DASHBOARD LAYOUT ==========
# Top row: Summary metrics
st.subheader("🌡️ Current Conditions Summary")
col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    st.metric("Temperature", f"{temperature_c}°C", f"{temperature_k:.0f} K")
with col2:
    st.metric("Pressure", f"{pressure_bar/1000:.1f} kbar", f"{pressure_bar:,} bar (FIXED)")
with col3:
    st.metric("ΔFMQ", f"{dfmq:.2f}")
with col4:
    log10_fo2_abs = compute_base_fmq_log10(temperature_k, pressure_bar) + dfmq
    st.metric("log₁₀(fO₂)", f"{log10_fo2_abs:.2f}")
with col5:
    st.metric("Salinities", f"{len(selected_salinities)}", "Selected")

st.markdown("---")

# ========== ENHANCED BALI MODEL PLOT WITH MIXING ==========
st.subheader("📈 Bali et al. 2012 Model: Mo/Ce vs U/Th with Mixing Model")

# Create the signature Mo/Ce vs U/Th plot with enhanced features
fig_model = go.Figure()

# Colors for different salinities and compositions
salinity_color_map = {
    0.001: '#1f77b4',  # Blue
    5.0: '#ff7f0e',    # Orange
    10.0: '#2ca02c',   # Green
    15.0: '#d62728',   # Red
    20.0: '#9467bd'    # Purple
}

salinity_name_map = {
    0.001: '0.001%',
    5.0: '5%',
    10.0: '10%',
    15.0: '15%',
    20.0: '20%'
}

# Composition colors and markers
composition_styles = {
    "Depleted Mantle (DM)": {'color': '#FF0000', 'symbol': 'circle', 'name': 'DM'},
    "Primitive Mantle (PM)": {'color': '#00AA00', 'symbol': 'square', 'name': 'PM'}
}

# Plot ALL salinities together (always visible)
all_compositions_data = []

for comp_name, results in all_composition_results.items():
    comp_style = composition_styles[comp_name]
    
    # Collect all points from ALL salinities for this composition
    all_x = []
    all_y = []
    all_colors = []
    all_sizes = []
    all_hover_text = []
    all_symbols = []
    
    # Process ALL available salinities (not just selected ones)
    all_available_salinities = [0.001, 5.0, 10.0, 15.0, 20.0]
    
    for salinity in all_available_salinities:
        salinity_results = [r for r in results if r['salinity_wt'] == salinity]
        if salinity_results:
            result = salinity_results[0]
            if np.isfinite(result['mo_ce_ratio']) and np.isfinite(result['u_th_ratio']) and result['mo_ce_ratio'] > 0 and result['u_th_ratio'] > 0:
                all_x.append(result['u_th_ratio'])
                all_y.append(result['mo_ce_ratio'])
                all_colors.append(salinity_color_map[salinity])
                all_sizes.append(14)
                all_symbols.append(comp_style['symbol'])
                all_hover_text.append(f"{comp_style['name']}: {salinity_name_map[salinity]} NaCl<br>U/Th: {result['u_th_ratio']:.3f}<br>Mo/Ce: {result['mo_ce_ratio']:.3f}")
    
    # Plot all salinities for this composition as a single trace with multiple colors
    if all_x:
        fig_model.add_trace(go.Scatter(
            x=all_x,
            y=all_y,
            mode='markers',
            name=f"{comp_style['name']} All Salinities",
            marker=dict(
                size=all_sizes, 
                color=all_colors,
                symbol=all_symbols[0],  # Use composition symbol
                line=dict(color='black', width=2)
            ),
            text=all_hover_text,
            hovertemplate='%{text}<extra></extra>',
            showlegend=True
        ))
    
    # Add evolution path connecting all fluid compositions
    if len(all_x) > 1:
        # Sort points by U/Th for smooth connection
        combined_points = list(zip(all_x, all_y))
        combined_points.sort(key=lambda p: p[0])
        
        fig_model.add_trace(go.Scatter(
            x=[p[0] for p in combined_points],
            y=[p[1] for p in combined_points],
            mode='lines',
            name=f'{comp_style["name"]} Evolution Path',
            line=dict(color=comp_style['color'], width=2, dash='dash'),
            showlegend=True,
            hoverinfo='skip'
        ))

# Add salinity color legend traces
salinity_legend_added = False
for salinity in selected_salinities:
    fig_model.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=12, color=salinity_color_map[salinity], symbol='circle',
                   line=dict(color='black', width=1)),
        name=f'{salinity_name_map[salinity]} NaCl',
        showlegend=True,
        hoverinfo='skip'
    ))

# Add mantle reference points
for comp_name in selected_compositions:
    comp_data = CONSTANTS.mantle_compositions[comp_name]
    comp_style = composition_styles[comp_name]
    
    # Calculate mantle ratios
    mantle_mo_ce = comp_data['MO'] / comp_data['Ce']
    mantle_u_th = comp_data['U'] / comp_data['Th']
    
    fig_model.add_trace(go.Scatter(
        x=[mantle_u_th],
        y=[mantle_mo_ce],
        mode='markers',
        name=f'{comp_style["name"]} Source',
        marker=dict(size=20, color=comp_style['color'], symbol='star',
                    line=dict(color='white', width=3)),
        showlegend=True
    ))

# Add enhanced mixing lines controlled by individual salinity checkboxes
# Add mixing lines controlled by individual salinity checkboxes
all_available_salinities = [0.001, 5.0, 10.0, 15.0, 20.0]

# Process all available compositions (not hardcoded)
for comp_name in all_composition_results.keys():  # Use actual available compositions
    comp_style = composition_styles[comp_name]
    comp_data = CONSTANTS.mantle_compositions[comp_name]
    
    # Get mantle source ratios
    mantle_mo_ce = comp_data['MO'] / comp_data['Ce']
    mantle_u_th = comp_data['U'] / comp_data['Th']
    
    for salinity in all_available_salinities:
        # ALWAYS show all mixing lines - no checkbox dependency
        # Get fluid endmember for this specific salinity
        salinity_results = [r for r in all_composition_results[comp_name] if r['salinity_wt'] == salinity]
        if salinity_results:
                fluid_result = salinity_results[0]
                fluid_mo_ce = fluid_result['mo_ce_ratio']
                fluid_u_th = fluid_result['u_th_ratio']
                
                if np.isfinite(fluid_mo_ce) and np.isfinite(fluid_u_th) and np.isfinite(mantle_mo_ce) and np.isfinite(mantle_u_th):
                    # Create mixing line from mantle source to fluid endmember using custom fluid fractions
                    mixing_x = []
                    mixing_y = []
                    mixing_sizes = []
                    
                    for f_fluid in fluid_fractions:
                        f_mantle = 1 - f_fluid  # Fraction of mantle source
                        
                        # Linear mixing between mantle source and fluid endmember
                        mixed_u_th = f_fluid * fluid_u_th + f_mantle * mantle_u_th
                        mixed_mo_ce = f_fluid * fluid_mo_ce + f_mantle * mantle_mo_ce
                        
                        mixing_x.append(mixed_u_th)
                        mixing_y.append(mixed_mo_ce)
                        
                        # Highlight the 1% (0.01) fraction point
                        if f_fluid == highlight_fraction:
                            mixing_sizes.append(16)  # Much larger size for highlighting (was 12)
                        else:
                            mixing_sizes.append(6)   # Normal size
                    
                    # Add this salinity's mixing line
                    salinity_color = salinity_color_map[salinity]
                    salinity_name = salinity_name_map[salinity]
                    
                    fig_model.add_trace(go.Scatter(
                        x=mixing_x,
                        y=mixing_y,
                        mode='lines+markers',
                        name=f'{comp_style["name"]} Mix {salinity_name}',
                        line=dict(color=salinity_color, width=2, dash='dot'),
                        marker=dict(size=mixing_sizes, color=salinity_color, symbol='diamond',
                                   line=dict(color='black', width=0.5)),
                        showlegend=True,
                        hovertemplate=f'<b>{comp_style["name"]} - {salinity_name} NaCl Mixing</b><br>' +
                                     'U/Th: %{x:.3f}<br>' +
                                     'Mo/Ce: %{y:.3f}<br>' +
                                     '<extra></extra>',
                        opacity=0.8
                    ))

# Enhanced layout with dark mode and improved aesthetics
fig_model.update_layout(
    title=dict(
        text=f"Model Conditions: T={temperature_c}°C, P={pressure_bar/1000:.1f} kbar (FIXED), ΔFMQ={dfmq:.1f}<br>" +
             f"Modal: CPX={mode_cpx_norm:.2f}, GRT={mode_grt_norm:.2f}, RUT={mode_rut_norm:.2f}",
        font=dict(color='white', size=16)
    ),
    xaxis_title=dict(text="U/Th", font=dict(color='white', size=14)),
    yaxis_title=dict(text="Mo/Ce", font=dict(color='white', size=14)), 
    xaxis_type="log",
    yaxis_type="log",
    height=800,
    showlegend=True,
    legend=dict(
        x=0.02, y=0.98, 
        bgcolor="rgba(40,40,40,0.9)",
        bordercolor="white",
        borderwidth=1,
        font=dict(color='white')
    ),
    # FIXED AXIS RANGES AS REQUESTED
    xaxis=dict(
        range=[np.log10(0.1), np.log10(10)],
        showgrid=True,
        gridwidth=1,
        gridcolor='white',
        zeroline=True,
        zerolinewidth=2,
        zerolinecolor='white',
        tickfont=dict(color='white'),
        title=dict(font=dict(color='white'))
    ),
    yaxis=dict(
        range=[np.log10(0.01), np.log10(1)],
        showgrid=True,
        gridwidth=1,
        gridcolor='white',
        zeroline=True,
        zerolinewidth=2,
        zerolinecolor='white',
        tickfont=dict(color='white'),
        title=dict(font=dict(color='white'))
    ),
    plot_bgcolor='#2E2E2E',  # Dark gray background
    paper_bgcolor='#1E1E1E', # Darker paper background
    font=dict(color='white')
)

# Add reference line at Mo/Ce = 0.1
fig_model.add_hline(y=0.1, line_dash="dot", line_color="cyan", 
                   annotation_text="Mo/Ce = 0.1 Reference",
                   annotation=dict(font=dict(color='white', size=12)))

st.plotly_chart(fig_model, use_container_width=True)

# ========== DETAILED ANALYSIS SECTION ==========
st.markdown("---")

# Create three columns for detailed analysis
col1, col2, col3 = st.columns([1, 1, 1])

with col1:
    st.subheader("🪨 Modal Mineralogy")
    
    if total_modes > 0:
        # Pie chart
        fig_pie = go.Figure(data=[go.Pie(
            labels=['Clinopyroxene', 'Garnet', 'Rutile'], 
            values=[mode_cpx_norm, mode_grt_norm, mode_rut_norm],
            hole=0.4,
            marker_colors=['#FF6B6B', '#4ECDC4', '#45B7D1'],
            textinfo='label+percent',
            textposition='outside'
        )])
        fig_pie.update_layout(
            height=300, 
            showlegend=False, 
            margin=dict(t=20, b=20, l=20, r=20),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            font=dict(color='white')
        )
        st.plotly_chart(fig_pie, use_container_width=True)
        
        # Display normalized values
        st.write("**Normalized Modal Proportions:**")
        st.write(f"• CPX: {mode_cpx_norm:.3f}")
        st.write(f"• GRT: {mode_grt_norm:.3f}")  
        st.write(f"• RUT: {mode_rut_norm:.3f}")
    else:
        st.warning("⚠️ No modal mineralogy specified")

with col2:
    st.subheader("⚛️ Partition Coefficients")
    
    # Show partition coefficients for the first composition and highest salinity
    if all_composition_results:
        primary_composition = selected_compositions[0]
        highest_salinity = max(selected_salinities)
        current_result = next((r for r in all_composition_results[primary_composition] if r['salinity_wt'] == highest_salinity), None)
        
        if current_result:
            st.write(f"**{composition_styles[primary_composition]['name']} at {salinity_name_map[highest_salinity]} NaCl**")
            
            # Create partition coefficient comparison
            elements = ['Mo', 'W', 'U']
            cpx_values = [current_result['mo_results']['D_CPX'], 
                          current_result['w_results']['D_CPX'], 
                          current_result['u_results']['D_CPX']]
            grt_values = [current_result['mo_results']['D_GRT'],
                          current_result['w_results']['D_GRT'], 
                          current_result['u_results']['D_GRT']]
            rut_values = [current_result['mo_results']['D_RUT'],
                          CONSTANTS.partition_constants['W']['RUT'], 
                          current_result['u_results']['D_RUT']]
            
            fig_partition = go.Figure()
            fig_partition.add_trace(go.Bar(name='CPX', x=elements, y=cpx_values, marker_color='#FF6B6B'))
            fig_partition.add_trace(go.Bar(name='GRT', x=elements, y=grt_values, marker_color='#4ECDC4'))
            fig_partition.add_trace(go.Bar(name='RUT', x=elements, y=rut_values, marker_color='#45B7D1'))
            
            fig_partition.update_layout(
                title=dict(text=f"D(mineral/fluid)", font=dict(color='white')),
                xaxis_title=dict(text="Element", font=dict(color='white')), 
                yaxis_title=dict(text="Partition Coefficient", font=dict(color='white')),
                yaxis_type="log",
                height=300,
                barmode='group',
                margin=dict(t=40, b=40, l=40, r=40),
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(30,30,30,0.8)',
                font=dict(color='white'),
                xaxis=dict(tickfont=dict(color='white')),
                yaxis=dict(tickfont=dict(color='white'))
            )
            st.plotly_chart(fig_partition, use_container_width=True)

with col3:
    st.subheader("🧪 Fluid Endmembers")
    
    if all_composition_results:
        primary_composition = selected_compositions[0]
        highest_salinity = max(selected_salinities)
        current_result = next((r for r in all_composition_results[primary_composition] if r['salinity_wt'] == highest_salinity), None)
        
        if current_result:
            st.write(f"**{composition_styles[primary_composition]['name']} at {salinity_name_map[highest_salinity]} NaCl**")
            
            # Fluid endmember concentrations
            fluid_elements = ['Mo', 'Ce', 'W', 'U', 'Th', 'Nb', 'La']
            fluid_values = [current_result['mo_fluid'], current_result['ce_fluid'], 
                           current_result['w_fluid'], current_result['u_fluid'],
                           current_result['th_fluid'], current_result['nb_fluid'],
                           current_result['la_fluid']]
            
            fig_fluid = go.Figure()
            fig_fluid.add_trace(go.Bar(
                x=fluid_elements,
                y=fluid_values,
                marker_color=['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F', '#BB8FCE'],
                text=[f'{v:.2e}' if np.isfinite(v) else 'NaN' for v in fluid_values],
                textposition='outside'
            ))
            
            fig_fluid.update_layout(
                title=dict(text=f"Fluid Concentrations (ppm)", font=dict(color='white')),
                xaxis_title=dict(text="Element", font=dict(color='white')),
                yaxis_title=dict(text="Concentration (ppm)", font=dict(color='white')),
                yaxis_type="log",
                height=300,
                showlegend=False,
                margin=dict(t=40, b=40, l=40, r=40),
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(30,30,30,0.8)',
                font=dict(color='white'),
                xaxis=dict(tickfont=dict(color='white')),
                yaxis=dict(tickfont=dict(color='white'))
            )
            st.plotly_chart(fig_fluid, use_container_width=True)

# ========== COMPOSITION COMPARISON TABLE ==========
st.markdown("---")
st.subheader("🌍 Mantle Composition Database")

# Create comparison table of initial concentrations
composition_table_data = []
for comp_name in CONSTANTS.mantle_compositions.keys():
    comp_data = CONSTANTS.mantle_compositions[comp_name]
    selected_marker = "✅" if comp_name in selected_compositions else "⬜"
    composition_table_data.append({
        'Selected': selected_marker,
        'Composition': comp_name,
        'Reference': "Workman & Hart (2005)" if "DM" in comp_name else "Palme & O'Neill (2014)" if "PM" in comp_name else "Representative",
        'Mo_ppm': f"{comp_data['MO']:.4f}",
        'Ce_ppm': f"{comp_data['Ce']:.4f}",
        'W_ppm': f"{comp_data['W']:.4f}",
        'U_ppm': f"{comp_data['U']:.4f}",
        'Th_ppm': f"{comp_data['Th']:.4f}",
        'Nb_ppm': f"{comp_data['Nb']:.4f}",
        'La_ppm': f"{comp_data['La']:.4f}",
        'Mo/Ce': f"{comp_data['MO']/comp_data['Ce']:.4f}",
        'U/Th': f"{comp_data['U']/comp_data['Th']:.3f}"
    })

composition_df = pd.DataFrame(composition_table_data)
st.dataframe(composition_df, use_container_width=True)

# ========== DATA TABLE ==========
st.markdown("---")
st.subheader("📋 Complete Results Table (All Selected Compositions & Salinities)")

# Create comprehensive results table
table_data = []
for comp_name, results in all_composition_results.items():
    comp_short = composition_styles[comp_name]['name']
    for result in results:
        table_data.append({
            'Composition': comp_short,
            'NaCl_wt_%': result['salinity_wt'],
            'NaCl_molality': f"{result['nacl_molality']:.3f}",
            'Mo/Ce': f"{result['mo_ce_ratio']:.4f}" if np.isfinite(result['mo_ce_ratio']) else "NaN",
            'U/Th': f"{result['u_th_ratio']:.4f}" if np.isfinite(result['u_th_ratio']) else "NaN",
            'W/Th': f"{result['w_th_ratio']:.4f}" if np.isfinite(result['w_th_ratio']) else "NaN",
            'Nb/La': f"{result['nb_la_ratio']:.4f}" if np.isfinite(result['nb_la_ratio']) else "NaN",
            'Mo_fluid_ppm': f"{result['mo_fluid']:.2e}" if np.isfinite(result['mo_fluid']) else "NaN",
            'Ce_fluid_ppm': f"{result['ce_fluid']:.2e}" if np.isfinite(result['ce_fluid']) else "NaN"
        })

results_df = pd.DataFrame(table_data)
st.dataframe(results_df, use_container_width=True)

# ========== FOOTER ==========
st.markdown("---")
st.markdown("**Model Reference:** Bali et al. (2012) - Trace element partitioning between clinopyroxene, garnet, rutile and silicate melts")
st.markdown("**Dashboard v5:** Enhanced implementation with mixing model • Fixed pressure: 26.1 kbar • Selectable salinities • White grid visualization")
st.markdown("**Features:** Fluid-melt mixing curves • All selected salinities displayed • Enhanced plot aesthetics • Real-time parameter variation")
