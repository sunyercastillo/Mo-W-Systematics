#!/bin/bash
# Simple bash script to launch the Bali Dashboard
# Usage: ./launch_dashboard.sh

echo "🚀 Launching Bali et al. 2012 Dashboard..."
cd "$(dirname "$0")"
python3 run_dashboard.py
