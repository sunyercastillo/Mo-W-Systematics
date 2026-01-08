#!/usr/bin/env python3
"""
Simple wrapper for automated DM-fluid mixing analysis.
Just run this script and it will do everything automatically!
"""

import os
import sys
from datetime import datetime
from pathlib import Path

def main():
    print("🔬 DM-FLUID MIXING AUTOMATED ANALYSIS")
    print("="*50)
    print("Enhanced Bali et al. 2012 model with full automation")
    print("="*50)
    
    # Default paths
    script_dir = Path(__file__).parent
    input_file = script_dir / "input.xlsx"
    
    # Create timestamp-based output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = script_dir / f"results_{timestamp}"
    
    # Check if input file exists
    if not input_file.exists():
        print(f"❌ Input file not found: {input_file}")
        print("\nPlease ensure you have an 'input.xlsx' file in this directory")
        print("with your rock composition data.")
        return
    
    print(f"📂 Input file: {input_file}")
    print(f"📁 Output directory: {output_dir}")
    print("\n🚀 Starting automated analysis...")
    
    # Run the analysis
    cmd = f'python "{script_dir}/Bali_v3.py" --auto --in_xlsx "{input_file}" --output_dir "{output_dir}"'
    
    try:
        result = os.system(cmd)
        if result == 0:
            print("\n✅ ANALYSIS COMPLETED SUCCESSFULLY!")
            print(f"\n📊 Results saved to: {output_dir}")
            print("\nGenerated files:")
            print("  🖼️  Publication-ready plots (PNG + PDF)")
            print("  📋 Organized data tables (CSV)")
            print("  📈 Mixing analysis for all salinity conditions")
            print("  📊 Summary statistics")
            
            # Show key output files
            if output_dir.exists():
                plot_png = output_dir / "DM_fluid_mixing_model.png"
                plot_pdf = output_dir / "DM_fluid_mixing_model.pdf"
                main_results = output_dir / "complete_results.csv"
                
                print("\n🎯 Key files:")
                if plot_png.exists():
                    print(f"  📊 Main plot: {plot_png}")
                if plot_pdf.exists():
                    print(f"  📄 PDF version: {plot_pdf}")
                if main_results.exists():
                    print(f"  📋 Complete data: {main_results}")
            
            print("\n🎉 Ready for publication!")
            
        else:
            print("❌ Analysis failed. Please check the error messages above.")
            
    except Exception as e:
        print(f"❌ Error running analysis: {e}")

if __name__ == "__main__":
    main()
