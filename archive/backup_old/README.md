# Bali et al. 2012 Trace Element Partitioning Dashboard

## Quick Start 🚀

### Option 1: Python Launcher (Recommended)
```bash
python3 run_dashboard.py
```

### Option 2: Bash Script
```bash
./launch_dashboard.sh
```

### Option 3: Direct Streamlit
```bash
streamlit run enhanced_dashboard_v5.py
```

## Features ✨

### **Fixed Parameters:**
- **Pressure**: 26,100 bar (fixed)
- **Starting Temperature**: 1000°C (adjustable)
- **Starting Modal Mineralogy**: 30/70/0 (CPX/GRT/RUT)

### **Automatic Multi-Salinity Plotting:**
- ✅ **All 5 salinities plotted automatically**
- 🔹 Freshwater (0.001% NaCl) - Blue
- 🔸 5% NaCl - Orange  
- 🔹 10% NaCl - Green
- 🔸 15% NaCl - Red
- 🔹 20% NaCl - Purple

### **Enhanced Visualization:**
- **Dark mode** with white grid lines
- **Custom fluid fractions**: [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
- **1% highlighting**: Larger markers for 0.01 fluid fraction
- **Mixing lines**: DM and PM mantle sources to fluid endmembers
- **Interactive controls**: Adjust temperature and modal mineralogy

## Files 📁

- `enhanced_dashboard_v5.py` - Main dashboard application
- `Bali_v5_dashboard.py` - Core calculation functions  
- `run_dashboard.py` - Python launcher script
- `launch_dashboard.sh` - Bash launcher script

## Usage Notes 📝

1. **No salinity selection needed** - all 5 salinities are always plotted
2. **1% fluid fraction highlighted** with larger markers on all mixing lines
3. **Starting conditions**: 1000°C, 30/70/0 modal mineralogy
4. **Fixed pressure** at 26,100 bar following Bali et al. 2012
5. **Dashboard auto-opens** in your default browser

## Controls 🎛️

### Sidebar Controls:
- **Temperature**: 400-1200°C (starts at 1000°C)
- **ΔFMQ**: Oxygen fugacity relative to FMQ buffer
- **Modal Mineralogy**: CPX/GRT/RUT proportions
- **Composition**: DM and PM mantle sources

### Main Plot:
- **Log-scale axes**: U/Th vs Mo/Ce
- **Fixed ranges**: Optimized for trace element data
- **Interactive legend**: Click to hide/show traces
- **Hover information**: Detailed data on mouse over

## Troubleshooting 🔧

### Dashboard won't start:
```bash
# Kill existing processes
pkill -f streamlit

# Try again
python3 run_dashboard.py
```

### Browser doesn't open automatically:
- Manually go to: http://localhost:8501

### Module import errors:
```bash
# Activate virtual environment
source .venv/bin/activate

# Install requirements
pip install streamlit plotly pandas numpy
```

---
**Based on:** Bali et al. (2012) trace element partitioning model  
**Dashboard Version:** v5.0  
**Last Updated:** August 2025
