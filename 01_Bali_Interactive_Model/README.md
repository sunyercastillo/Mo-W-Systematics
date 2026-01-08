# 🎯 01_Bali_Interactive_Model - Bali et al. 2012 Dashboard

## 📁 **File Structure Guide - What Each File Does**

### **🚀 MAIN ACTIVE FILES:**
- **`MAIN_streamlit_web_dashboard.py`** - Interactive web dashboard interface (USE THIS)
- **`CORE_scientific_calculations_engine.py`** - Core Bali et al. 2012 calculations engine
- **`LAUNCHER_python_script.py`** - Python launcher script (recommended way to start)
- **`LAUNCHER_bash_script.sh`** - Bash launcher script (alternative way to start)
- **`README.md`** - This documentation file

### **📚 REFERENCE FILES (Old Versions):**
- **`OLD_bali_dashboard_version5.py`** - Original Bali dashboard version 5
- **`OLD_enhanced_dashboard_v3_fixed_pressure.py`** - Enhanced version with fixed pressure
- **`OLD_enhanced_dashboard_v5.py`** - Enhanced dashboard version 5

### **⚙️ SYSTEM FILES:**
- **`.venv/`** - Python virtual environment with all dependencies
- **`.github/`** - GitHub configuration files
- **`__pycache__/`** - Python compiled cache files

## 🚀 Quick Start Guide

### Option 1: Python Launcher (Recommended)
```bash
python LAUNCHER_python_script.py
```

### Option 2: Bash Script
```bash
./LAUNCHER_bash_script.sh
```

### Option 3: Direct Streamlit
```bash
streamlit run MAIN_streamlit_web_dashboard.py
```

## 🧪 **What This Dashboard Does:**
- **Scientific Model**: Implements Bali et al. 2012 trace element partitioning calculations
- **Interactive Visualization**: Real-time parameter adjustment with immediate plot updates
- **Multi-Salinity Analysis**: 5 different NaCl concentrations (0.001%, 5%, 10%, 15%, 20%)
- **Custom Fluid Fractions**: [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
- **Dark Mode Visualization**: Professional scientific plotting theme

## 📊 **Key Features:**
- Molybdenum and Uranium partitioning calculations
- Interactive mixing line visualization
- Temperature and pressure parameter control
- Modal mineralogy adjustment sliders
- Real-time plot updates with parameter changes

## 🔧 **Technical Details:**
- **Framework**: Streamlit + Plotly for interactive web interface
- **Python Version**: 3.9.6 (isolated virtual environment)
- **Pressure Setting**: Fixed at 26,100 bar (following Bali et al. 2012)
- **Starting Conditions**: 1000°C temperature, 30/70/0 modal mineralogy

## 📈 **File Purposes Explained:**
- **MAIN_** prefix = Primary active file to use
- **CORE_** prefix = Essential calculation engine
- **LAUNCHER_** prefix = Scripts to start the dashboard
- **OLD_** prefix = Historical versions for reference

---
*For questions about the scientific model, refer to Bali et al. 2012 publication*
*For technical issues, check the OLD_ files for comparison*
