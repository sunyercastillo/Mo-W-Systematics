#!/bin/bash
# This is a clickable launcher for macOS
# Double-click this file to start the dashboard

# Get the directory where this script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "------------------------------------------------------------------"
echo "   Opening Bali et al 2012 Interactive Dashboard"
echo "------------------------------------------------------------------"

# Ensure the sub-script is executable
chmod +x "01_Bali_Interactive_Model/LAUNCHER_bash_script.sh"

# Run the main launcher
./01_Bali_Interactive_Model/LAUNCHER_bash_script.sh
