#!/usr/bin/env python3
"""
Interactive Bali et al. 2012 Model Dashboard

A Streamlit-based interactive dashboard for exploring trace element partitioning
under different P-T-fO2-salinity conditions.

Usage: streamlit run streamlit_dashboard.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

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
    page_title="Bali et al. 2012 Model",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ========== CONSTANTS ==========
CONSTANTS = BaliConstants()

# ========== DASHBOARD TITLE ==========
st.title("🧪 Bali et al. 2012 Interactive Model")
st.markdown("**Trace Element Partitioning: Mo, W, U in Fluid-Mineral Systems**")
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

# Salinity
st.sidebar.subheader("Fluid Composition")
salinity_options = [0.001, 5.0, 10.0, 15.0, 20.0]
salinity_wt_pct = st.sidebar.selectbox("NaCl (wt%)", salinity_options, index=2)

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

# ========== MAIN CALCULATIONS ==========
# Calculate fundamental properties
log10_fo2_fmq = compute_base_fmq_log10(temperature_k, pressure_bar)
log10_fo2_abs = log10_fo2_fmq + dfmq

# NaCl properties
nacl_props = calculate_nacl_properties(salinity_wt_pct)
nacl_molality = nacl_props['molality']
log10_nacl = nacl_props['log10_nacl']

# Element partitioning
mo_results = calculate_molybdenum_partitioning(log10_fo2_abs, log10_nacl, temperature_k)
w_results = calculate_tungsten_partitioning(log10_fo2_abs, temperature_k)
u_results = calculate_uranium_partitioning(log10_fo2_abs, nacl_molality)

# Modal properties for bulk calculations
modal_props = {'cpx': mode_cpx_norm, 'grt': mode_grt_norm, 'rut': mode_rut_norm}

# Bulk partition coefficients
mo_bulk = calculate_bulk_partitioning(
    {'CPX': mo_results['D_CPX'], 'GRT': mo_results['D_GRT'], 'RUT': mo_results['D_RUT']},
    modal_props
)

w_bulk = calculate_bulk_partitioning(
    {'CPX': w_results['D_CPX'], 'GRT': w_results['D_GRT'], 'RUT': CONSTANTS.partition_constants['W']['RUT']},
    modal_props
)

u_bulk = calculate_bulk_partitioning(
    {'CPX': u_results['D_CPX'], 'GRT': u_results['D_GRT'], 'RUT': u_results['D_RUT']},
    modal_props
)

# Constant bulk coefficients
ce_bulk = calculate_bulk_partitioning(CONSTANTS.partition_constants['CE'], modal_props)
th_bulk = calculate_bulk_partitioning(CONSTANTS.partition_constants['TH'], modal_props)
nb_bulk = calculate_bulk_partitioning(CONSTANTS.partition_constants['NB'], modal_props)
la_bulk = calculate_bulk_partitioning(CONSTANTS.partition_constants['LA'], modal_props)

# Fluid endmembers
mo_fluid = calculate_fluid_endmember(c0_mo, mo_bulk)
ce_fluid = calculate_fluid_endmember(c0_ce, ce_bulk)
w_fluid = calculate_fluid_endmember(c0_w, w_bulk)
u_fluid = calculate_fluid_endmember(c0_u, u_bulk)
th_fluid = calculate_fluid_endmember(c0_th, th_bulk)
nb_fluid = calculate_fluid_endmember(c0_nb, nb_bulk)
la_fluid = calculate_fluid_endmember(c0_la, la_bulk)

# ========== MAIN DASHBOARD LAYOUT ==========
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("🌡️ Thermodynamic Summary")
    
    # Create summary metrics
    metrics_col1, metrics_col2, metrics_col3 = st.columns(3)
    
    with metrics_col1:
        st.metric("Temperature", f"{temperature_c}°C", f"{temperature_k:.1f} K")
        st.metric("Pressure", f"{pressure_bar:,} bar", f"{pressure_bar/1000:.1f} kbar")
    
    with metrics_col2:
        st.metric("ΔFMQ", f"{dfmq:.2f}", "log units")
        st.metric("log₁₀(fO₂)ₐᵦₛ", f"{log10_fo2_abs:.2f}")
    
    with metrics_col3:
        st.metric("NaCl", f"{salinity_wt_pct}%", f"{nacl_molality:.3f} m")
        st.metric("log₁₀(NaCl)", f"{log10_nacl:.2f}")

with col2:
    st.subheader("⚛️ Partition Coefficients")
    
    # Create partition coefficient comparison
    elements = ['Mo', 'W', 'U']
    cpx_values = [mo_results['D_CPX'], w_results['D_CPX'], u_results['D_CPX']]
    grt_values = [mo_results['D_GRT'], w_results['D_GRT'], u_results['D_GRT']]
    rut_values = [mo_results['D_RUT'], CONSTANTS.partition_constants['W']['RUT'], u_results['D_RUT']]
    
    fig_partition = go.Figure()
    fig_partition.add_trace(go.Bar(name='CPX', x=elements, y=cpx_values, marker_color='#FF6B6B'))
    fig_partition.add_trace(go.Bar(name='GRT', x=elements, y=grt_values, marker_color='#4ECDC4'))
    fig_partition.add_trace(go.Bar(name='RUT', x=elements, y=rut_values, marker_color='#45B7D1'))
    
    fig_partition.update_layout(
        title="Mineral/Fluid Partition Coefficients",
        xaxis_title="Element", 
        yaxis_title="D (mineral/fluid)",
        yaxis_type="log",
        height=350,
        barmode='group'
    )
    st.plotly_chart(fig_partition, use_container_width=True)

# ========== FOOTER ==========
st.markdown("---")
st.markdown("**Model Reference:** Bali et al. (2012) - Trace element partitioning model")
st.markdown("**Dashboard:** Interactive implementation with modular calculations")
