#!/usr/bin/env python3
"""
Quick test script for the enhanced automation features.
"""

import sys
import os
from pathlib import Path

# Add the directory containing Bali_v3.py to the Python path
sys.path.insert(0, str(Path(__file__).parent))

try:
    from Bali_v3 import automated_batch_analysis, load_excel, normalize_headers
    print("✓ Successfully imported enhanced Bali_v3 module")
except ImportError as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

def test_automation():
    """Test the automated workflow with sample data."""
    
    print("🧪 TESTING AUTOMATED WORKFLOW")
    print("="*50)
    
    # Default input file path
    input_file = "/Users/polsuka/Desktop/Bali et al 2012 model/input.xlsx"
    
    if not os.path.exists(input_file):
        print(f"⚠️  Input file not found: {input_file}")
        print("Please ensure your input.xlsx file exists in the expected location.")
        return False
    
    try:
        # Load and normalize data
        df_raw = load_excel(input_file, 0)
        df = normalize_headers(df_raw)
        print(f"✓ Loaded {len(df)} samples from input file")
        
        # Test the automated analysis
        results = automated_batch_analysis(
            df, 
            output_dir="test_output",
            create_plots=True
        )
        
        if results:
            print("\n✅ AUTOMATION TEST PASSED!")
            print(f"Generated {len(results.get('tables', {}))} data tables")
            print(f"Generated {len(results.get('mixing_tables', []))} mixing tables")
            if results.get('plots'):
                print(f"Generated {len(results['plots'])} plots")
            return True
        else:
            print("❌ AUTOMATION TEST FAILED!")
            return False
            
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_automation()
    if success:
        print("\n🎉 Ready for automated DM-fluid mixing analysis!")
        print("\nTo run full automation, use:")
        print("  python Bali_v3.py --auto --output_dir results")
    else:
        print("\n🔧 Please check the configuration and try again.")
