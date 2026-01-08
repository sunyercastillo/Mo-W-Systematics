#!/usr/bin/env python3
"""
Simple launcher script for the Bali et al. 2012 Dashboard
This script automatically opens the dashboard in your browser
"""

import subprocess
import webbrowser
import time
import sys
import os

def run_dashboard():
    """Launch the Streamlit dashboard and open it in browser"""
    
    print("🚀 Starting Bali et al. 2012 Dashboard...")
    print("=" * 50)
    
    # Change to the correct directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    # Start Streamlit in background
    try:
        # Kill any existing streamlit processes
        print("🧹 Cleaning up old processes...")
        subprocess.run(["pkill", "-f", "streamlit"], capture_output=True)
        time.sleep(2)  # Wait longer for port to free
        
        print("📊 Launching Streamlit server...")
        # Use sys.executable to run streamlit module correctly from current environment
        proc = subprocess.Popen([
            sys.executable, "-m", "streamlit", "run", "MAIN_streamlit_web_dashboard.py",
            "--server.headless", "true",
            "--server.port", "8501",
            "--server.address", "localhost"
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Wait a moment for server to start
        print("⏳ Waiting for server to initialize...")
        time.sleep(5)
        
        # Open browser
        url = "http://localhost:8501"
        print(f"🌐 Opening dashboard in browser: {url}")
        webbrowser.open(url)
        
        print("✅ Dashboard is now running!")
        print("📝 Press Ctrl+C to stop the dashboard")
        print("=" * 50)
        
        # Keep the process running
        proc.wait()
        
    except KeyboardInterrupt:
        print("\n🛑 Stopping dashboard...")
        proc.terminate()
        subprocess.run(["pkill", "-f", "streamlit"], capture_output=True)
        print("✅ Dashboard stopped!")
    except Exception as e:
        print(f"❌ Error starting dashboard: {e}")
        sys.exit(1)

if __name__ == "__main__":
    run_dashboard()
