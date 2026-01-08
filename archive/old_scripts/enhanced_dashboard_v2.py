#!/usr/bin/env python3
"""
Enhanced Bali et al. 2012 Model Dashboard v2

Complete interactive implementation with fixed mantle compositions,
multi-salinity analysis, fO2 variation, and modal mineralogy exploration.

Usage: streamlit run enhanced_dashboard_v2.py
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
    page_title="Bali et al. 2012 Enhanced Model v2",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ========== CONSTANTS ==========
CONSTANTS = BaliConstants()

# ========== DASHBOARD TITLE ==========
st.title("🧪 Bali et al. 2012 Enhanced Interactive Model v2")
st.markdown("**Complete Trace Element Partitioning Analysis: Fixed Mantle Compositions • Multi-Salinity • fO₂ • Modal Variations**")
st.markdown("---")

# ========== SIDEBAR: INPUT PARAMETERS ==========
st.sidebar.header("📊 Model Parameters")

# Core thermodynamic conditions
st.sidebar.subheader("Thermodynamic Conditions")
temperature_c = st.sidebar.slider("Temperature (°C)", 400, 1200, 800, 25)
pressure_bar = st.sidebar.slider("Pressure (bar)", 1000, 50000, 20000, 1000)
dfmq = st.sidebar.slider("ΔFMQ (log units)", -4.0, 4.0, 0.0, 0.1)

# Convert temperature
temperature_k = temperature_c + 273.15

# Salinity (for current point highlighting)
st.sidebar.subheader("Current Salinity Focus")
salinity_options = [0.001, 5.0, 10.0, 15.0, 20.0]
salinity_wt_pct = st.sidebar.selectbox("Highlight NaCl (wt%)", salinity_options, index=2)

# ========== MANTLE COMPOSITION SELECTION ==========
st.sidebar.subheader("🌍 Initial Mantle Composition")
st.sidebar.markdown("*Select one or more compositions to compare:*")

# Checkbox system for mantle selection
use_dm = st.sidebar.checkbox("✅ Depleted Mantle (DM)", value=True, help="Workman & Hart (2005)")
use_pm = st.sidebar.checkbox("✅ Primitive Mantle (PM)", value=False, help="Palme & O'Neill (2014)")
use_em = st.sidebar.checkbox("✅ Enriched Mantle (EM)", value=False, help="Representative enriched composition")

# Display selected compositions
selected_compositions = []
if use_dm:
    selected_compositions.append("Depleted Mantle (DM)")
if use_pm:
    selected_compositions.append("Primitive Mantle (PM)")
if use_em:
    selected_compositions.append("Enriched Mantle (EM)")

# If no composition selected, default to DM
if not selected_compositions:
    selected_compositions = ["Depleted Mantle (DM)"]
    use_dm = True
    st.sidebar.warning("⚠️ No composition selected - defaulting to DM")

# Show composition details
st.sidebar.markdown("**Selected Compositions:**")
for comp_name in selected_compositions:
    comp_data = CONSTANTS.mantle_compositions[comp_name]
    
    # Calculate key ratios
    mo_ce = comp_data['MO'] / comp_data['Ce']
    u_th = comp_data['U'] / comp_data['Th']
    
    st.sidebar.markdown(f"• **{comp_name.split('(')[1].split(')')[0]}**")
    st.sidebar.markdown(f"  - Mo/Ce: {mo_ce:.4f}")
    st.sidebar.markdown(f"  - U/Th: {u_th:.3f}")
    
    with st.sidebar.expander(f"View detailed concentrations"):
        col1, col2 = st.columns(2)
        with col1:
            st.write(f"Mo: {comp_data['MO']:.4f} ppm")
            st.write(f"Ce: {comp_data['Ce']:.4f} ppm")
            st.write(f"W: {comp_data['W']:.4f} ppm")
            st.write(f"U: {comp_data['U']:.4f} ppm")
        with col2:
            st.write(f"Th: {comp_data['Th']:.4f} ppm")
            st.write(f"Nb: {comp_data['Nb']:.4f} ppm")
            st.write(f"La: {comp_data['La']:.4f} ppm")

# Modal mineralogy
st.sidebar.subheader("Modal Mineralogy")
mode_cpx = st.sidebar.slider("Clinopyroxene", 0.0, 1.0, 0.3, 0.01)
mode_grt = st.sidebar.slider("Garnet", 0.0, 1.0, 0.2, 0.01)
mode_rut = st.sidebar.slider("Rutile", 0.0, 1.0, 0.05, 0.01)

# Normalize modal proportions
total_modes = mode_cpx + mode_grt + mode_rut
if total_modes > 0:
    mode_cpx_norm = mode_cpx / total_modes
    mode_grt_norm = mode_grt / total_modes  
    mode_rut_norm = mode_rut / total_modes
else:
    mode_cpx_norm = mode_grt_norm = mode_rut_norm = 0

st.sidebar.markdown(f"**Normalized:** CPX={mode_cpx_norm:.2f}, GRT={mode_grt_norm:.2f}, RUT={mode_rut_norm:.2f}")

# ========== CALCULATE ALL SALINITIES FOR ALL COMPOSITIONS ==========
@st.cache_data
def calculate_all_compositions_salinities(temp_k, press_bar, dfmq_val, cpx_norm, grt_norm, rut_norm, 
                                        composition_names):
    """Calculate results for all salinity values and all selected mantle compositions - cached for performance."""
    
    salinity_values = [0.001, 5.0, 10.0, 15.0, 20.0]
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
        
        for sal_wt in salinity_values:
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
all_composition_results = calculate_all_compositions_salinities(
    temperature_k, pressure_bar, dfmq, mode_cpx_norm, mode_grt_norm, mode_rut_norm,
    selected_compositions
)

# Get current salinity result for the first selected composition (for detailed display)
primary_composition = selected_compositions[0]
current_result = next(r for r in all_composition_results[primary_composition] if r['salinity_wt'] == salinity_wt_pct)

# ========== MAIN DASHBOARD LAYOUT ==========
# Top row: Summary metrics
st.subheader("🌡️ Current Conditions Summary")
col1, col2, col3, col4, col5, col6 = st.columns(6)

with col1:
    st.metric("Temperature", f"{temperature_c}°C", f"{temperature_k:.0f} K")
with col2:
    st.metric("Pressure", f"{pressure_bar/1000:.1f} kbar", f"{pressure_bar:,} bar")
with col3:
    st.metric("ΔFMQ", f"{dfmq:.2f}")
with col4:
    log10_fo2_abs = compute_base_fmq_log10(temperature_k, pressure_bar) + dfmq
    st.metric("log₁₀(fO₂)", f"{log10_fo2_abs:.2f}")
with col5:
    st.metric("Modal Total", f"{total_modes:.2f}", "Sum of phases")
with col6:
    st.metric("Compositions", f"{len(selected_compositions)}", "Selected")

st.markdown("---")

# ========== SIGNATURE BALI MODEL PLOT ==========
st.subheader("📈 Bali et al. 2012 Model: Mo/Ce vs U/Th (All Salinities & Compositions)")

# Create the signature Mo/Ce vs U/Th plot with all salinity lines for all compositions
fig_model = go.Figure()

# Colors for different salinities and compositions
salinity_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
salinity_names = ['0.001%', '5%', '10%', '15%', '20%']

# Composition colors and markers
composition_styles = {
    "Depleted Mantle (DM)": {'color': '#FF0000', 'symbol': 'circle', 'name': 'DM'},
    "Primitive Mantle (PM)": {'color': '#00AA00', 'symbol': 'square', 'name': 'PM'},
    "Enriched Mantle (EM)": {'color': '#0000FF', 'symbol': 'diamond', 'name': 'EM'}
}

# Plot each composition
for comp_name, results in all_composition_results.items():
    comp_style = composition_styles[comp_name]
    
    # Collect valid points for connecting lines
    valid_points = []
    for i, result in enumerate(results):
        if np.isfinite(result['mo_ce_ratio']) and np.isfinite(result['u_th_ratio']) and result['mo_ce_ratio'] > 0 and result['u_th_ratio'] > 0:
            valid_points.append({
                'x': result['u_th_ratio'],
                'y': result['mo_ce_ratio'],
                'salinity': result['salinity_wt'],
                'color': salinity_colors[i % len(salinity_colors)],
                'salinity_name': salinity_names[i]
            })
    
    # Sort by U/Th ratio for connecting lines
    valid_points.sort(key=lambda p: p['x'])
    
    # Add connecting line for this composition
    if len(valid_points) > 1:
        fig_model.add_trace(go.Scatter(
            x=[p['x'] for p in valid_points],
            y=[p['y'] for p in valid_points],
            mode='lines',
            name=f'{comp_style["name"]} Trend',
            line=dict(color=comp_style['color'], width=3, dash='dash'),
            showlegend=True
        ))
    
    # Plot points for each salinity in this composition
    for i, result in enumerate(results):
        if np.isfinite(result['mo_ce_ratio']) and np.isfinite(result['u_th_ratio']):
            # Show individual points for focused composition only
            if comp_name == primary_composition:
                fig_model.add_trace(go.Scatter(
                    x=[result['u_th_ratio']],
                    y=[result['mo_ce_ratio']],
                    mode='markers',
                    name=f"{comp_style['name']}: {salinity_names[i]} NaCl",
                    marker=dict(size=12, color=salinity_colors[i % len(salinity_colors)], 
                               symbol=comp_style['symbol'],
                               line=dict(color='white', width=2)),
                    showlegend=True
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
                    line=dict(color='black', width=3)),
        showlegend=True
    ))

# Highlight current salinity for primary composition
if np.isfinite(current_result['mo_ce_ratio']) and np.isfinite(current_result['u_th_ratio']):
    fig_model.add_trace(go.Scatter(
        x=[current_result['u_th_ratio']],
        y=[current_result['mo_ce_ratio']],
        mode='markers',
        name=f'FOCUS: {salinity_wt_pct}% NaCl ({composition_styles[primary_composition]["name"]})',
        marker=dict(size=25, color='gold', symbol='star',
                    line=dict(color='black', width=3)),
        showlegend=True
    ))

fig_model.update_layout(
    title=f"Model Conditions: T={temperature_c}°C, P={pressure_bar/1000:.1f} kbar, ΔFMQ={dfmq:.1f}<br>" +
          f"Modal: CPX={mode_cpx_norm:.2f}, GRT={mode_grt_norm:.2f}, RUT={mode_rut_norm:.2f}<br>" +
          f"Selected: {', '.join([composition_styles[comp]['name'] for comp in selected_compositions])}",
    xaxis_title="U/Th",
    yaxis_title="Mo/Ce", 
    xaxis_type="log",
    yaxis_type="log",
    height=700,
    showlegend=True,
    legend=dict(x=0.02, y=0.98, bgcolor="rgba(255,255,255,0.8)"),
    xaxis=dict(range=[np.log10(0.05), np.log10(20)]),
    yaxis=dict(range=[np.log10(0.001), np.log10(2)]),
    template="plotly_white"
)

# Add reference line at Mo/Ce = 0.1
fig_model.add_hline(y=0.1, line_dash="dot", line_color="gray", 
                   annotation_text="Mo/Ce = 0.1 Reference")

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
        fig_pie.update_layout(height=300, showlegend=False, margin=dict(t=20, b=20, l=20, r=20))
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
    st.write(f"**Focus: {composition_styles[primary_composition]['name']} at {salinity_wt_pct}% NaCl**")
    
    # Create partition coefficient comparison for current salinity
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
        title=f"D(mineral/fluid)",
        xaxis_title="Element", 
        yaxis_title="Partition Coefficient",
        yaxis_type="log",
        height=300,
        barmode='group',
        margin=dict(t=40, b=40, l=40, r=40)
    )
    st.plotly_chart(fig_partition, use_container_width=True)

with col3:
    st.subheader("🧪 Fluid Endmembers")
    st.write(f"**Focus: {composition_styles[primary_composition]['name']} at {salinity_wt_pct}% NaCl**")
    
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
        title=f"Fluid Concentrations (ppm)",
        xaxis_title="Element",
        yaxis_title="Concentration (ppm)",
        yaxis_type="log",
        height=300,
        showlegend=False,
        margin=dict(t=40, b=40, l=40, r=40)
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
st.subheader("📋 Complete Results Table (Selected Compositions)")

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
st.markdown("**Dashboard v2:** Enhanced implementation with fixed mantle compositions • Multi-composition comparison • Interactive parameter variation")
st.markdown("**Features:** Fixed mantle compositions • Multi-composition analysis • Real-time parameter variation • Interactive visualizations")

# ========== HELP SECTION ==========
with st.expander("📚 How to Add More Mantle Compositions"):
    st.markdown("""
    ### Adding New Mantle Compositions
    
    To add additional mantle compositions (e.g., HIMU, EM1, EM2, OIB), follow these steps:
    
    #### 1. Update the BaliConstants class in `Bali_v5_dashboard.py`:
    
    ```python
    # Add new composition data (example: HIMU)
    self.himu_concentrations = {
        "MO": 0.040, "Ce": 2.200, "W": 0.015, "U": 0.025, 
        "Th": 0.065, "Nb": 1.200, "La": 1.100
    }
    
    # Add to the mantle compositions database
    self.mantle_compositions = {
        "Depleted Mantle (DM)": self.dm_concentrations,
        "Primitive Mantle (PM)": self.pm_concentrations,
        "Enriched Mantle (EM)": self.em_concentrations,
        "HIMU": self.himu_concentrations  # Add new composition here
    }
    ```
    
    #### 2. Update the composition styles in this dashboard:
    
    ```python
    composition_styles = {
        "Depleted Mantle (DM)": {'color': '#FF0000', 'symbol': 'circle', 'name': 'DM'},
        "Primitive Mantle (PM)": {'color': '#00AA00', 'symbol': 'square', 'name': 'PM'},
        "Enriched Mantle (EM)": {'color': '#0000FF', 'symbol': 'diamond', 'name': 'EM'},
        "HIMU": {'color': '#FF00FF', 'symbol': 'triangle-up', 'name': 'HIMU'}  # Add styling
    }
    ```
    
    #### 3. Add checkbox in sidebar:
    
    ```python
    use_himu = st.sidebar.checkbox("✅ HIMU", value=False, help="High μ (U/Pb) mantle")
    
    # Update selected_compositions logic
    if use_himu:
        selected_compositions.append("HIMU")
    ```
    
    ### Literature References for Mantle Compositions:
    - **DM**: Workman & Hart (2005) - Major and trace element composition of the depleted MORB mantle (DMM)
    - **PM**: Palme & O'Neill (2014) - Cosmochemical constraints on mantle composition  
    - **EM**: Zindler & Hart (1986) - Chemical geodynamics of mantle endmembers
    - **HIMU**: Weaver (1991) - The origin of ocean island basalt endmember compositions
    - **MORB**: Gale et al. (2013) - The mean composition of ocean ridge basalts
    
    ### Current Compositions Implemented:
    """)
    
    # Display current compositions in help
    for comp_name, comp_data in CONSTANTS.mantle_compositions.items():
        st.write(f"**{comp_name}:**")
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.write(f"Mo: {comp_data['MO']:.4f}")
            st.write(f"Ce: {comp_data['Ce']:.4f}")
        with col2:
            st.write(f"W: {comp_data['W']:.4f}")
            st.write(f"U: {comp_data['U']:.4f}")
        with col3:
            st.write(f"Th: {comp_data['Th']:.4f}")
            st.write(f"Nb: {comp_data['Nb']:.4f}")
        with col4:
            st.write(f"La: {comp_data['La']:.4f}")
            st.write(f"Mo/Ce: {comp_data['MO']/comp_data['Ce']:.4f}")
