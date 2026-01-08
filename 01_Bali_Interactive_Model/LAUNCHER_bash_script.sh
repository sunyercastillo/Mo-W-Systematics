#!/bin/bash
# Simple bash script to launch the Bali Dashboard
# Usage: ./launch_dashboard.sh

cd "$(dirname "$0")"

# Check for virtual environment in parent directory
if [ -d "../.venv" ]; then
    PYTHON_CMD="../.venv/bin/python"
# Check for local venv
elif [ -d ".venv" ]; then
    PYTHON_CMD=".venv/bin/python"
elif [ -d "venv" ]; then
    PYTHON_CMD="venv/bin/python"
else
    PYTHON_CMD="python3"
fi

echo "🚀 Launching Bali et al. 2012 Dashboard using $PYTHON_CMD..."
"$PYTHON_CMD" LAUNCHER_python_script.py
