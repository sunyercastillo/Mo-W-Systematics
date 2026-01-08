#!/usr/bin/env python3
"""
Enhanced Bali et al. 2012 Model Dashboard

Complete interactive implementation with multi-salinity analysis,
fO2 variation, and modal mineralogy exploration.

Usage: streamlit run enhanced_dashboard.py
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
    page_title="Bali et al. 2012 Enhanced Model",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ========== CONSTANTS ==========
CONSTANTS = BaliConstants()

# ========== DASHBOARD TITLE ==========
st.title("🧪 Bali et al. 2012 Enhanced Interactive Model")
st.markdown("**Complete Trace Element Partitioning Analysis: Multi-Salinity, fO₂, and Modal Variations**")
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

# Initial concentrations
st.sidebar.subheader("Initial Concentrations (ppm)")
c0_mo = st.sidebar.number_input("Mo", value=0.025, format="%.3f")
c0_ce = st.sidebar.number_input("Ce", value=0.772, format="%.3f")
c0_w = st.sidebar.number_input("W", value=0.0024, format="%.4f")
c0_u = st.sidebar.number_input("U", value=0.0047, format="%.4f")
c0_th = st.sidebar.number_input("Th", value=0.0137, format="%.4f")
c0_nb = st.sidebar.number_input("Nb", value=0.21, format="%.3f")
c0_la = st.sidebar.number_input("La", value=0.234, format="%.3f")

# ========== CALCULATE ALL SALINITIES ==========
@st.cache_data
def calculate_all_salinities(temp_k, press_bar, dfmq_val, cpx_norm, grt_norm, rut_norm, 
                           c0_mo_val, c0_ce_val, c0_w_val, c0_u_val, c0_th_val, c0_nb_val, c0_la_val):
    """Calculate results for all salinity values - cached for performance."""
    
    salinity_values = [0.001, 5.0, 10.0, 15.0, 20.0]
    results = []
    
    # Base calculations (same for all salinities)
    log10_fo2_fmq_base = compute_base_fmq_log10(temp_k, press_bar)
    log10_fo2_abs_base = log10_fo2_fmq_base + dfmq_val
    
    modal_props_base = {'cpx': cpx_norm, 'grt': grt_norm, 'rut': rut_norm}
    
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
        
        results.append({
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
            'u_results': u_res
        })
    
    return results

# Calculate results for all salinities
all_salinity_results = calculate_all_salinities(
    temperature_k, pressure_bar, dfmq, mode_cpx_norm, mode_grt_norm, mode_rut_norm,
    c0_mo, c0_ce, c0_w, c0_u, c0_th, c0_nb, c0_la
)

# Get current salinity result for detailed display
current_result = next(r for r in all_salinity_results if r['salinity_wt'] == salinity_wt_pct)

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
    st.metric("Focus Salinity", f"{salinity_wt_pct}%", "NaCl")

st.markdown("---")

# ========== SIGNATURE BALI MODEL PLOT ==========
st.subheader("📈 Bali et al. 2012 Model: Mo/Ce vs U/Th (All Salinities)")

# Create the signature Mo/Ce vs U/Th plot with all salinity lines
fig_model = go.Figure()

# Colors for different salinities
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
salinity_names = ['0.001%', '5%', '10%', '15%', '20%']

# Collect valid points for connecting lines
valid_points = []
for i, result in enumerate(all_salinity_results):
    if np.isfinite(result['mo_ce_ratio']) and np.isfinite(result['u_th_ratio']) and result['mo_ce_ratio'] > 0 and result['u_th_ratio'] > 0:
        valid_points.append({
            'x': result['u_th_ratio'],
            'y': result['mo_ce_ratio'],
            'salinity': result['salinity_wt'],
            'color': colors[i % len(colors)],
            'name': salinity_names[i]
        })

# Sort by U/Th ratio for connecting lines
valid_points.sort(key=lambda p: p['x'])

# Add connecting line
if len(valid_points) > 1:
    fig_model.add_trace(go.Scatter(
        x=[p['x'] for p in valid_points],
        y=[p['y'] for p in valid_points],
        mode='lines',
        name='Salinity Trend',
        line=dict(color='black', width=2, dash='dash'),
        showlegend=True
    ))

# Plot points for each salinity
for i, result in enumerate(all_salinity_results):
    if np.isfinite(result['mo_ce_ratio']) and np.isfinite(result['u_th_ratio']):
        fig_model.add_trace(go.Scatter(
            x=[result['u_th_ratio']],
            y=[result['mo_ce_ratio']],
            mode='markers',
            name=f"{salinity_names[i]} NaCl",
            marker=dict(size=14, color=colors[i % len(colors)], 
                       line=dict(color='white', width=2)),
            showlegend=True
        ))

# Add depleted mantle point
dm_mo_ce = CONSTANTS.dm_concentrations['MO'] / CONSTANTS.dm_concentrations['Ce']
dm_u_th = CONSTANTS.dm_concentrations['U'] / CONSTANTS.dm_concentrations['Th']

fig_model.add_trace(go.Scatter(
    x=[dm_u_th],
    y=[dm_mo_ce],
    mode='markers',
    name='Depleted Mantle (DM)',
    marker=dict(size=18, color='grey', symbol='square', 
                line=dict(color='black', width=3)),
    showlegend=True
))

# Highlight current salinity
current_idx = next(i for i, r in enumerate(all_salinity_results) if r['salinity_wt'] == salinity_wt_pct)
if np.isfinite(current_result['mo_ce_ratio']) and np.isfinite(current_result['u_th_ratio']):
    fig_model.add_trace(go.Scatter(
        x=[current_result['u_th_ratio']],
        y=[current_result['mo_ce_ratio']],
        mode='markers',
        name=f'FOCUS: {salinity_wt_pct}% NaCl',
        marker=dict(size=25, color='red', symbol='star',
                    line=dict(color='black', width=3)),
        showlegend=True
    ))

fig_model.update_layout(
    title=f"Model Conditions: T={temperature_c}°C, P={pressure_bar/1000:.1f} kbar, ΔFMQ={dfmq:.1f}<br>" +
          f"Modal: CPX={mode_cpx_norm:.2f}, GRT={mode_grt_norm:.2f}, RUT={mode_rut_norm:.2f}",
    xaxis_title="U/Th",
    yaxis_title="Mo/Ce", 
    xaxis_type="log",
    yaxis_type="log",
    height=600,
    showlegend=True,
    legend=dict(x=0.02, y=0.98, bgcolor="rgba(255,255,255,0.8)"),
    xaxis=dict(range=[np.log10(0.1), np.log10(10)]),
    yaxis=dict(range=[np.log10(0.01), np.log10(1)]),
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
        title=f"D(mineral/fluid) at {salinity_wt_pct}% NaCl",
        xaxis_title="Element", 
        yaxis_title="Partition Coefficient",
        yaxis_type="log",
        height=350,
        barmode='group',
        margin=dict(t=40, b=40, l=40, r=40)
    )
    st.plotly_chart(fig_partition, use_container_width=True)

with col3:
    st.subheader("🧪 Fluid Endmembers")
    
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
        title=f"Fluid Endmember (ppm) at {salinity_wt_pct}% NaCl",
        xaxis_title="Element",
        yaxis_title="Concentration (ppm)",
        yaxis_type="log",
        height=350,
        showlegend=False,
        margin=dict(t=40, b=40, l=40, r=40)
    )
    st.plotly_chart(fig_fluid, use_container_width=True)

# ========== DATA TABLE ==========
st.markdown("---")
st.subheader("📋 Complete Results Table")

# Create comprehensive results table
table_data = []
for result in all_salinity_results:
    table_data.append({
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
st.markdown("**Dashboard:** Enhanced interactive implementation with full multi-salinity analysis")
st.markdown("**Features:** Real-time parameter variation • Multi-salinity modeling • Interactive visualizations")
