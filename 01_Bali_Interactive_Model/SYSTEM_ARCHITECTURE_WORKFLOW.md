# 🏗️ Complete System Architecture & Workflow Tree

## 📋 **SYSTEM OVERVIEW**
Interactive scientific dashboard implementing the Bali et al. (2012) trace element partitioning model with real-time parameter exploration and comprehensive data visualization.

---

## 🧮 **1. LOGIC MODULES: Scientific Calculation Engine**

### **1.1 Core Constants & Configuration (`BaliConstants` Class)**
```
📦 Constants
├── 🔬 Physical Constants
│   ├── molar_masses: {MO: 95.95, W: 183.84, NACL: 58.44}
│   └── log_scales: Pre-computed logarithmic scales for performance
├── 🌍 Mantle Compositions (ppm)
│   ├── dm_concentrations: Depleted Mantle (Workman & Hart 2005)
│   │   └── {MO: 0.025, Ce: 0.772, W: 0.0024, U: 0.0047, Th: 0.0137, Nb: 0.21, La: 0.234}
│   ├── pm_concentrations: Primitive Mantle (Palme & O'Neill 2014)
│   │   └── {MO: 0.047, Ce: 1.7529, W: 0.012, U: 0.0229, Th: 0.0849, Nb: 0.595, La: 0.6832} X
│   └── mantle_compositions: Combined database X
├── ⚖️ Partition Coefficients (mineral/fluid)
│   ├── partition_constants: Fixed coefficients for CE, TH, NB, LA, W
│   └── partition_denominators: Temperature-dependent denominators for MO, W, U
└── 📐 Model Equation Coefficients
    ├── FMQ: {A: -25096.3, B: 8.735, C: 0.11}
    ├── MO: {fO2: 0.435, NaCl: 0.42, temp: -1.8, const: 4.8}
    ├── W: {fO2: 0.07, temp: -4.7236, const: 4.4271}
    └── U: {fO2: 0.1433, NaCl: 0.594, const_cg: 2.681, const_rut: 1.7954}
```

### **1.2 Thermodynamic Foundation**
```
🌡️ Thermodynamic Calculations
├── compute_base_fmq_log10(T_K, P_bar)
│   ├── INPUT: Temperature (K), Pressure (bar)
│   ├── EQUATION: log10(fO2_FMQ) = A/T + B + C*P/T
│   ├── CONSTANTS: A=-25096.3, B=8.735, C=0.11
│   └── OUTPUT: Base FMQ log10(fO2) reference
└── Absolute fO2 Calculation
    ├── INPUT: Base FMQ + ΔFMQ offset
    ├── EQUATION: log10(fO2_abs) = log10(fO2_FMQ) + ΔFMQ
    └── OUTPUT: Absolute oxygen fugacity for partitioning calculations
```

### **1.3 NaCl Solution Properties**
```
🧂 calculate_nacl_properties(salinity_wt_pct)
├── INPUT: NaCl weight percentage (0.001% to 20%)
├── PROCESS: Molality Conversion
│   ├── EQUATION: molality = 1000 * wt% / (M_NaCl * (100 - wt%))
│   ├── M_NaCl = 58.44 g/mol
│   └── Guard against division by zero at 100% salinity
├── PROCESS: Logarithmic Transformation
│   ├── EQUATION: log10_NaCl = log10(molality)
│   └── Guard against log(0) for ultra-low salinities
└── OUTPUT: {molality: float, log10_nacl: float}
```

### **1.4 Element-Specific Partitioning Calculations**

#### **1.4.1 Molybdenum Partitioning**
```
⚛️ calculate_molybdenum_partitioning(log10_fo2_abs, log10_nacl, T_K)
├── INPUT VALIDATION: Check finite values for all inputs
├── PROCESS: LogMO Calculation
│   ├── EQUATION: LogMO = 0.435*log10(fO2) + 0.42*log10(NaCl) - 1.8*(1000/T) + 4.8
│   └── Temperature dependence: -1.8*(1000/T_K) term
├── PROCESS: E_MO Calculation (with overflow protection)
│   ├── Combined log = LogMO + log10(M_MO * 1000)
│   ├── Overflow guard: -308 ≤ combined_log ≤ 308
│   └── EQUATION: E_MO = (10^LogMO) * M_MO * 1000
├── PROCESS: Partition Coefficients
│   ├── D_CPX = 40.0 / E_MO
│   ├── D_GRT = 12.0 / E_MO
│   └── D_RUT = 87670.0 / E_MO
└── OUTPUT: {LogMO, E_MO, D_CPX, D_GRT, D_RUT}
```

#### **1.4.2 Tungsten Partitioning**
```
⚛️ calculate_tungsten_partitioning(log10_fo2_abs, T_K)
├── INPUT VALIDATION: Check finite values
├── PROCESS: LogW Calculation
│   ├── EQUATION: LogW = 0.07*log10(fO2) - 4.7236*(1000/T) + 4.4271
│   └── Strong temperature dependence: -4.7236*(1000/T_K)
├── PROCESS: E_W Calculation (with overflow protection)
│   ├── Combined log = LogW + log10(M_W * 1000)
│   ├── Overflow guard: -308 ≤ combined_log ≤ 308
│   └── EQUATION: E_W = (10^LogW) * M_W * 1000
├── PROCESS: Partition Coefficients
│   ├── D_CPX = 60.0 / E_W
│   ├── D_GRT = 12.0 / E_W
│   └── D_RUT = 1.250 (constant)
└── OUTPUT: {LogW, E_W, D_CPX, D_GRT, D_RUT}
```

#### **1.4.3 Uranium Partitioning**
```
⚛️ calculate_uranium_partitioning(log10_fo2_abs, nacl_molality)
├── INPUT VALIDATION: Check finite values
├── PROCESS: Dual LogU Calculations
│   ├── LogU_cg = 2.681 + 0.1433*log10(fO2) + 0.594*NaCl_molality
│   └── LogU_rut = 1.7954 + 0.1433*log10(fO2) + 0.594*NaCl_molality
├── PROCESS: E_U Calculations (with overflow protection)
│   ├── E_U_cg = 10^LogU_cg (for CPX/GRT)
│   └── E_U_rut = 10^LogU_rut (for RUT)
├── PROCESS: Partition Coefficients
│   ├── D_CPX = 11.0 / E_U_cg
│   ├── D_GRT = 40.0 / E_U_cg
│   └── D_RUT = 94.0 / E_U_rut
└── OUTPUT: {LogU_cg, LogU_rut, E_U_cg, E_U_rut, D_CPX, D_GRT, D_RUT}
```

### **1.5 Bulk Partitioning & Fluid Endmembers**
```
⚖️ calculate_bulk_partitioning(partition_coeffs, modal_props)
├── INPUT: Element-specific partition coefficients + modal mineralogy
├── PROCESS: Modal Weighted Average
│   ├── EQUATION: D_bulk = Σ(D_mineral * mode_mineral)
│   ├── D_bulk = (D_CPX * f_CPX) + (D_GRT * f_GRT) + (D_RUT * f_RUT)
│   └── Modal fractions normalized: f_CPX + f_GRT + f_RUT = 1
└── OUTPUT: Bulk partition coefficient (float)

💧 calculate_fluid_endmember(initial_concentration, bulk_partition_coeff)
├── INPUT: Initial mantle concentration (ppm) + bulk D value
├── PROCESS: Fluid Endmember Calculation
│   ├── EQUATION: C_fluid = C_initial / (D_bulk + (1 - D_bulk) * f)
│   ├── At f=1 (pure fluid): C_fluid = C_initial / 1 = C_initial
│   └── Physical meaning: Concentration in equilibrium fluid
└── OUTPUT: Fluid endmember concentration (ppm)
```

### **1.6 Mixing Model Implementation**
```
🌊 Mixing Model Calculations
├── INPUT: Fluid endmember concentrations + DM constants + fluid fractions
├── PROCESS: Binary Mixing Equation
│   ├── EQUATION: C_mix = (C_fluid * f) + (C_DM * (1-f))
│   ├── f = fluid fraction (0 to 1)
│   ├── C_fluid = fluid endmember concentration
│   └── C_DM = depleted mantle concentration
├── PROCESS: Elemental Ratio Calculations
│   ├── Mo/Ce = C_Mo_mix / C_Ce_mix
│   ├── U/Th = C_U_mix / C_Th_mix
│   ├── W/Th = C_W_mix / C_Th_mix
│   ├── Mo/W = C_Mo_mix / C_W_mix
│   └── Nb/La = C_Nb_mix / C_La_mix
└── OUTPUT: Complete mixing curves for all fluid fractions
```

### **1.7 Master Calculation Engine**
```
🚀 process_all_vectorized(df, sal_wts)
├── INPUT PREPARATION
│   ├── DataFrame with rock compositions and conditions
│   ├── List of salinity values [0.001, 5.0, 10.0, 15.0, 20.0]
│   └── Expected columns validation and default values
├── THERMODYNAMIC FOUNDATION
│   ├── Compute base FMQ for all rows: log10_fO2_FMQ
│   ├── Add ΔFMQ offset: log10_fO2_abs = log10_fO2_FMQ + dFMQ
│   └── Handle NaN values in ΔFMQ gracefully
├── SALINITY PROCESSING
│   ├── Replicate DataFrame for each salinity value
│   ├── Calculate NaCl molality and log10_NaCl for each row
│   └── Maintain pandas dtypes throughout processing
├── ELEMENT-SPECIFIC CALCULATIONS (Vectorized)
│   ├── Molybdenum: LogMO → E_MO → D_CPX/GRT/RUT
│   ├── Tungsten: LogW → E_W → D_CPX/GRT/RUT
│   ├── Uranium: LogU_cg/rut → E_U_cg/rut → D_CPX/GRT/RUT
│   └── Overflow protection for all exponential calculations
├── BULK PARTITIONING
│   ├── Modal-weighted averages for all elements
│   ├── Fixed coefficients: CE, TH, NB, LA
│   └── Temperature-dependent: MO, W, U
├── FLUID ENDMEMBERS
│   ├── Calculate equilibrium fluid concentrations
│   ├── Safe division with NaN handling
│   └── Physical validation of results
└── OUTPUT: Complete results DataFrame with all intermediate and final values
```

---

## 📊 **2. DATA MODULES: Input/Output Processing**

### **2.1 Excel Data Input Pipeline**
```
📋 Excel Data Loading & Processing
├── load_excel(path, sheet)
│   ├── INPUT: Excel file path and sheet number
│   ├── PROCESS: pandas.read_excel()
│   └── OUTPUT: Raw DataFrame
├── normalize_headers(df)
│   ├── INPUT: Raw DataFrame with various header formats
│   ├── PROCESS: Header Normalization Pipeline
│   │   ├── Convert to lowercase and strip whitespace
│   │   ├── Remove special characters: °, space, tabs, parentheses
│   │   ├── Normalize delta symbols: Δ → d, δ → d
│   │   └── Map to standardized column names via HEADER_MAP
│   ├── HEADER_MAP Mappings:
│   │   ├── Identifiers: "input"/"id"/"sample" → "rock_id"
│   │   ├── Pressure: "pbar"/"p(bar)" → "P_bar"
│   │   ├── Temperature: "temperaturek"/"tk" → "T_K"
│   │   ├── ΔFMQ: "dfmq"/"Δfmq" → "dFMQ"
│   │   ├── Modal: "modegrt"/"grt" → "mode_grt"
│   │   └── Elements: "mo" → "C0_MO", "ce" → "C0_Ce", etc.
│   └── OUTPUT: Standardized DataFrame with consistent headers
└── Data Validation
    ├── Check for required columns
    ├── Handle missing values with defaults
    └── Type conversion and error handling
```

### **2.2 Results Table Generation**
```
📈 Results Processing Pipeline
├── TABLE 1: Input Parameters (Hidden in Dashboard)
│   ├── rock_id, P_bar, T_K, dFMQ
│   ├── NaCl_wt_pct, NaCl_m, log10_NaCl
│   ├── log10_fO2_FMQ, log10_fO2_abs
│   └── Initial concentrations: C0_MO, C0_Ce, C0_W, C0_U, C0_Th, C0_Nb, C0_La
├── TABLE 2: Per-Mineral Partition Coefficients
│   ├── Temperature-dependent: MO_D_CPX/GRT/RUT, W_D_CPX/GRT, U_D_CPX/GRT/RUT
│   └── Salinity breakdown for each NaCl concentration
├── TABLE 3: Bulk Partition Coefficients & Fluid Endmembers
│   ├── Modal-weighted D_bulk values: MO_Dbulk, CE_Dbulk, W_Dbulk, etc.
│   └── Equilibrium fluid concentrations: MO_F_EM, CE_F_EM, W_F_EM, etc.
├── TABLE 4: Endmember Ratios
│   ├── DM ratios: DM_MO_CE, DM_U_TH, DM_W_TH, DM_MO_W, DM_NB_LA
│   └── Fluid ratios: F_MO_CE, F_U_TH, F_W_TH, F_MO_W, F_NB_LA
└── TABLE 5.x: DM-F Mixing Model Results (Separate for each salinity)
    ├── Table 5.1: 0.001% NaCl mixing curves
    ├── Table 5.2: 5% NaCl mixing curves
    ├── Table 5.3: 10% NaCl mixing curves
    ├── Table 5.4: 15% NaCl mixing curves
    ├── Table 5.5: 20% NaCl mixing curves
    └── Each table contains:
        ├── Fluid fractions: f = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.9]
        ├── Mixed concentrations: MO_DM_F_MIX, CE_DM_F_MIX, etc.
        └── Mixed ratios: MO_CE_MIX, U_TH_MIX, W_TH_MIX, MO_W_MIX, NB_LA_MIX
```

### **2.3 Data Flow Architecture**
```
🔄 Data Flow Control
├── PRIMARY PATH: Excel Row Processing
│   ├── User selects Excel row → Extract source parameters
│   ├── Apply interactive T and ΔFMQ overrides
│   ├── Process through core engine: process_all_vectorized()
│   └── Generate all tables and mixing data
├── FALLBACK PATH: Manual Composition Selection
│   ├── User selects composition type (DM/PM)
│   ├── Use interactive T, ΔFMQ, and fixed modal mineralogy
│   ├── Create temporary DataFrame → Process through core engine
│   └── Generate results for selected compositions
└── CACHING STRATEGY
    ├── @st.cache_data for Excel loading functions
    ├── @st.cache_data for expensive calculations
    └── Cache invalidation on parameter changes
```

---

## 🖥️ **3. UI MODULES: Dashboard Interface Components**

### **3.1 Page Configuration & Styling**
```
🎨 Dashboard Setup
├── Page Configuration
│   ├── Title: "Bali et al. 2012 Enhanced Model v5"
│   ├── Icon: 🧪, Layout: wide, Sidebar: expanded
│   └── Dark theme optimized for scientific visualization
├── CSS Styling (Dark Theme)
│   ├── Background: #0E1117 (main), #262730 (sidebar)
│   ├── Text: #FFFFFF with contrast optimization
│   ├── DataFrames: #1E1E1E background
│   ├── Metrics: #262730 with border styling
│   └── Interactive elements: Optimized for visibility
└── Global Helper Functions
    ├── safe_div(): Division with zero protection
    └── safe_ratio(): Pandas-aware ratio calculations
```

### **3.2 Source Data Management**
```
📊 Source Parameters Section
├── Excel Data Display Table
│   ├── INPUT: Raw Excel data with row numbers (starting from row 2)
│   ├── COLUMNS: All source parameters including T_K, P_bar, dFMQ, compositions
│   ├── FEATURES: Scrollable, full-width display
│   └── PURPOSE: Complete visibility of available data
└── Source Selection Interface
    ├── Radio Button Selection
    │   ├── OPTIONS: All Excel rows with rock_id labels
    │   ├── FORMAT: "Row X: Rock_ID" for easy identification
    │   └── KEY: "excel_row_selection" for state management
    └── Selected Row Details Panel
        ├── Rock ID display
        ├── Temperature: °C and K with precision
        ├── Pressure: kbar conversion
        ├── ΔFMQ: Oxidation state
        └── Modal composition: CPX, GRT, RUT fractions
```

### **3.3 Interactive Parameter Controls**
```
🎛️ Sidebar Parameter Controls
├── Temperature Control Section
│   ├── Dynamic Slider Setup
│   │   ├── Range: 400-1200°C, Step: 25°C
│   │   ├── Default: Excel T_K value (when available)
│   │   ├── Key: f"temp_slider_{selected_excel_row}" (dynamic reset)
│   │   └── Help text: Shows Excel default value
│   ├── Parameter Preservation
│   │   ├── Pressure: Maintained from Excel (read-only display)
│   │   ├── ΔFMQ: Initial value from Excel (before interactive override)
│   │   └── Modal composition: Preserved from source
│   └── Change Indicators
│       ├── Temperature change: ±XX°C from Excel (when modified)
│       └── Real-time calculation triggering
├── ΔFMQ Control Section
│   ├── Dynamic Slider Setup
│   │   ├── Range: -3.0 to +3.0, Step: 0.05
│   │   ├── Default: Excel dFMQ value (when available)
│   │   ├── Key: f"dfmq_slider_{selected_excel_row}" (dynamic reset)
│   │   └── Help text: Shows Excel default value
│   ├── Parameter Integration
│   │   ├── Source parameters display (Pressure from Excel)
│   │   ├── Temperature change indicator (if modified)
│   │   └── ΔFMQ change indicator (if modified)
│   └── Real-time Updates
│       ├── All calculations update immediately on slider change
│       └── Both sliders reset when Excel row selection changes
└── Salinity Controls
    ├── Checkbox Array: [0.001%, 5%, 10%, 15%, 20%]
    ├── Default selections: All salinities enabled
    ├── Purpose: Control mixing line visibility in plots
    └── Dynamic legend generation based on selections
```

### **3.4 Current Conditions Summary**
```
📊 Conditions Summary Dashboard
├── Metrics Row (4 columns)
│   ├── Temperature: Current °C and K values
│   ├── Pressure: kbar and bar with (FIXED) indicator
│   ├── ΔFMQ: Current oxidation state (2 decimal precision)
│   └── log₁₀(fO₂): Calculated absolute oxygen fugacity
├── Model Conditions Box
│   ├── Styled container with border and background
│   ├── All current parameters in compact format
│   ├── Modal composition display: CPX, GRT, RUT fractions
│   └── Centered layout for professional appearance
└── Real-time Updates
    ├── All values update automatically with slider changes
    ├── log₁₀(fO₂) recalculates: base_FMQ + interactive_ΔFMQ
    └── Visual feedback for parameter modifications
```

### **3.5 Visualization System**
```
📈 Enhanced Plotting System
├── Plot Configuration Array
│   ├── Plot 1: Mo/Ce vs U/Th, Range: [0.1,10] x [0.01,10]
│   ├── Plot 2: W/Th vs U/Th, Range: [0.1,10] x [0.01,10]  
│   ├── Plot 3: Mo/W vs U/Th, Range: [0.1,10] x [0.1,100]
│   └── Plot 4: Mo/W vs Nb/La, Range: [0.1,10] x [0.1,100]
├── Data Sources & Styling
│   ├── Table 5.x Mixing Data: One curve per salinity
│   ├── Okabe-Ito Colorblind-Safe Palette:
│   │   ├── 0.001% NaCl: #0173B2 (Blue)
│   │   ├── 5% NaCl: #DE8F05 (Orange)
│   │   ├── 10% NaCl: #029E73 (Bluish green)
│   │   ├── 15% NaCl: #CC78BC (Reddish purple)
│   │   └── 20% NaCl: #CA3542 (Vermillion)
│   └── Special Markers:
│       ├── Mixing lines: markers+lines, size=6, width=1.5
│       ├── 1% fluid fraction: X markers, size=6
│       └── DM (f=0): White stars, size=12, black outline
├── Plot Layout & Styling
│   ├── Logarithmic axes with specific ranges
│   ├── White plot background for readability
│   ├── Grid: Dashed gray lines, moderate opacity
│   ├── Axis styling: Black lines, white text, custom tick values
│   └── Dark theme integration: Paper background matches UI
├── Interactive Features
│   ├── Hover templates: Detailed information on mouseover
│   ├── Custom data: Fluid fraction values for each point
│   ├── Unified legend: Separate from individual plots
│   └── Real-time updates: All plots refresh with parameter changes
└── Layout Management
    ├── 2x2 grid layout using Streamlit columns
    ├── Individual plots: 450px height for optimal viewing
    ├── Centered display with responsive design
    └── Legend below plots: Comprehensive symbol explanation
```

### **3.6 Legend & Documentation**
```
🗂️ Legend System
├── Dynamic Legend Generation
│   ├── Salinity colors: Only show selected salinities
│   ├── Composition symbols: Based on active selections
│   ├── Special markers: X for 1%, star for DM/PM
│   └── Clean white theme with professional styling
├── Legend Components
│   ├── Salinity Section: Color-coded entries for each selected salinity
│   ├── Separator: Visual division between sections
│   ├── Special Markers: DM star, 1% X explanation
│   └── Composition entries: Dynamic based on selections
└── Styling Features
    ├── Background: #f8f9fa with border and shadow
    ├── Flex layout: Responsive wrapping of legend items
    ├── Color swatches: 18px squares with borders
    └── Typography: Bold labels, readable text sizes
```

### **3.7 Detailed Calculation Tables**
```
📋 Table Display System
├── Table Organization
│   ├── Table 2: Per-Mineral Partition Coefficients
│   │   ├── Columns: NaCl_wt_pct, MO_D_CPX/GRT/RUT, W_D_CPX/GRT, U_D_CPX/GRT/RUT
│   │   └── Purpose: Show temperature-dependent partitioning
│   ├── Table 3: Bulk Partition Coefficients & Fluid Endmembers
│   │   ├── Columns: D_bulk values + F_EM concentrations for all elements
│   │   └── Purpose: Modal-weighted results and equilibrium fluids
│   ├── Table 4: Endmember Ratios
│   │   ├── DM ratios: Depleted mantle elemental ratios
│   │   ├── Fluid ratios: Equilibrium fluid elemental ratios
│   │   └── Purpose: Reference values for mixing model endpoints
│   └── Tables 5.1-5.5: DM-F Mixing Model Results
│       ├── Separate subtable for each salinity
│       ├── Columns: f, concentrations, elemental ratios
│       └── Purpose: Complete mixing curves for all fluid fractions
├── Table Features
│   ├── Full-width display: width="stretch"
│   ├── Interactive scrolling for large datasets
│   ├── Automatic formatting: Precision appropriate for each column
│   └── Real-time updates: All tables refresh with parameter changes
├── Data Processing
│   ├── Safe division functions: Handle NaN and zero values
│   ├── DM constants integration: From BaliConstants class
│   ├── Mixing calculations: Exact replication of core engine logic
│   └── Validation: Cross-check with plotting data
└── Integration with Interactive Controls
    ├── Temperature changes: Affect all partition coefficients
    ├── ΔFMQ changes: Modify all fO2-dependent calculations
    ├── Row selection: Complete recalculation with new source parameters
    └── Salinity selection: Filter displayed mixing tables
```

### **3.8 State Management & Performance**
```
⚡ Performance Optimization
├── Caching Strategy
│   ├── @st.cache_data: Excel loading and processing functions
│   ├── Cache invalidation: On parameter changes
│   └── Memory management: Efficient DataFrame operations
├── State Management
│   ├── Dynamic keys: Slider reset when row selection changes
│   ├── Session state: Preserve user selections
│   ├── Reactive updates: Automatic recalculation on input changes
│   └── Error handling: Graceful degradation for invalid inputs
├── Data Flow Optimization
│   ├── Single calculation engine: Avoid code duplication
│   ├── Vectorized operations: Efficient pandas processing
│   ├── Minimal recomputation: Only recalculate when necessary
│   └── Progressive loading: Display results as they become available
└── User Experience
    ├── Real-time feedback: Immediate visual response to changes
    ├── Progress indicators: For long calculations
    ├── Error messages: Clear guidance for invalid inputs
    └── Responsive design: Optimal viewing on different screen sizes
```

```
🚀 USER INTERACTION → CALCULATION → VISUALIZATION PIPELINE

1. USER INPUT
   ├── Select Excel Row → Extract source parameters
   ├── Adjust Temperature → Override Excel T_K
   ├── Adjust ΔFMQ → Override Excel dFMQ
   └── Select Salinities → Control mixing line visibility

2. PARAMETER PROCESSING
   ├── Dynamic key generation → Reset sliders on row change
   ├── Parameter validation → Ensure finite, reasonable values
   ├── Source preservation → Keep Excel P, modal composition, concentrations
   └── Interactive overrides → Apply T and ΔFMQ modifications

3. CORE CALCULATIONS (process_all_vectorized)
   ├── Thermodynamic foundation → FMQ + ΔFMQ → Absolute fO2
   ├── NaCl properties → Weight% → Molality → log10(NaCl)
   ├── Element partitioning → Temperature & fO2 dependent coefficients
   ├── Bulk partitioning → Modal-weighted averaging
   ├── Fluid endmembers → Equilibrium concentrations
   └── Mixing models → DM-fluid binary mixing for all fractions

4. RESULTS GENERATION
   ├── Table 2 → Per-mineral partition coefficients
   ├── Table 3 → Bulk coefficients and fluid endmembers  
   ├── Table 4 → Endmember ratios (DM vs Fluid)
   └── Tables 5.1-5.5 → Mixing curves for each salinity

5. VISUALIZATION PIPELINE
   ├── Plot data extraction → Table 5.x mixing results
   ├── Color/symbol mapping → Okabe-Ito palette + special markers
   ├── Plot generation → 4 elemental ratio plots with logarithmic axes
   ├── Legend creation → Dynamic based on selections
   └── Real-time updates → All components refresh on parameter changes

6. USER FEEDBACK
   ├── Summary metrics → Current T, P, ΔFMQ, log₁₀(fO₂)
   ├── Change indicators → Show modifications from Excel defaults
   ├── Interactive plots → Hover information, mixing curves
   ├── Detailed tables → Complete calculation transparency
   └── Professional layout → Dark theme, responsive design
```

---

