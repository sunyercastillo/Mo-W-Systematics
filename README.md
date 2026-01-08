# Bali et al. (2012) Interactive Model Dashboard

This repository contains an interactive dashboard for modeling trace element partitioning (Mo, W, U, Th) in subduction zones, based on the scientific framework of Bali et al. (2012).

## 📊 Features
- **Real-time Interaction:** Adjust Temperatures, fO₂, and mineralogy instantly.
- **Scientific Visualization:** Log-Log plots of element ratios (Mo/Ce, W/Th, etc.) and mixing lines for various salinities.
- **Publication Ready:** Export high-resolution plots (PNG/SVG) and access detailed data tables.

## 📁 Project Structure

```
├── 01_Bali_Interactive_Model/    # Main Application Source Code
│   ├── MAIN_streamlit_web_dashboard.py   # Entry point
│   ├── CORE_scientific_calculations_engine.py  # Physics/Calc Logic
│   └── sample_input_data.xlsx    # Default input parameters
├── requirements.txt              # Python dependencies
└── README.md                     # Documentation
```

## 🚀 Getting Started

### 1. Prerequisites
- Python 3.8+
- Git

### 2. Installation

Clone the repository:
```bash
git clone https://github.com/YOUR_USERNAME/Bali_Model_Project.git
cd Bali_Model_Project
```

Create a virtual environment (recommended):
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows use: .venv\Scripts\activate
```

Install dependencies:
```bash
pip install -r requirements.txt
```

### 3. Running the Dashboard
Navigate to the application folder and run Streamlit:

```bash
cd 01_Bali_Interactive_Model
streamlit run MAIN_streamlit_web_dashboard.py
```

The dashboard will open in your browser at `http://localhost:8501`.

## 🤝 Collaboration
This project is structured for GitHub.
- **Master Branch:** Contains the stable, working version.
- **Updates:** Pull the latest changes with `git pull`.
- **Contribution:** Create a new branch for features, then create a Pull Request.

## 🛠 Tech Stack
- **Frontend:** Streamlit
- **Visualization:** Plotly
- **Data:** Pandas, NumPy, OpenPyXL
