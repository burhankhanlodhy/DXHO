#!/usr/bin/env python3
"""
Example usage of SIMBAD enhancement functionality
This script demonstrates how to enhance astronomical data with SIMBAD information.
"""

import pandas as pd
import numpy as np
import os
from simbad_simple import SIMBADEnhancer

def create_sample_data():
    """Create sample astronomical data similar to your CSV format."""
    np.random.seed(42)  # For reproducible results
    
    # Sample coordinates around NGC 7380 region (your test data area)
    base_ra = 341.0  # degrees
    base_dec = 57.3  # degrees
    n_objects = 20
    
    # Generate sample data
    data = {
        'filename': ['test.fits'] * n_objects,
        'id': range(1, n_objects + 1),
        'ra': np.random.normal(base_ra, 0.5, n_objects),  # Small spread around base
        'dec': np.random.normal(base_dec, 0.3, n_objects),
        'flux': np.random.exponential(1000, n_objects),
        'mag_calibrated': np.random.normal(16, 2, n_objects),
        'gaia_source_id': [''] * n_objects,  # Mostly empty for demo
        'simbad_main_id': [''] * n_objects,  # All empty - need to fill
        'otype': [''] * n_objects,
        'sp_type': [''] * n_objects, 
        'rv_value': [np.nan] * n_objects,
        'distance_result': [np.nan] * n_objects
    }
    
    # Add a few known objects with existing SIMBAD data
    data['simbad_main_id'][0] = 'HD 216763'
    data['otype'][0] = 'Star'
    data['sp_type'][0] = 'G5V'
    
    return pd.DataFrame(data)

def demonstrate_enhancement():
    """Demonstrate the SIMBAD enhancement process."""
    
    print("=" * 60)
    print("SIMBAD Data Enhancement Example")
    print("=" * 60)
    print()
    
    # Step 1: Create sample data
    print("1. Creating sample astronomical data...")
    df = create_sample_data()
    sample_csv = 'sample_data.csv'
    df.to_csv(sample_csv, index=False)
    print(f"   + Created {sample_csv} with {len(df)} astronomical objects")
    
    # Show initial state
    missing_count = (df['simbad_main_id'] == '').sum()
    print(f"   + {missing_count} objects need SIMBAD data")
    print()
    
    # Step 2: Initialize enhancer
    print("2. Setting up SIMBAD enhancer...")
    enhancer = SIMBADEnhancer()
    
    # Configure for demo (smaller batches, shorter delays)
    enhancer.batch_size = 10
    enhancer.delay_between_batches = 1.0
    
    print("   + Configured for gentle API usage")
    print("   + Batch size: 10 objects")
    print("   + Delay between batches: 1.0 seconds")
    print()
    
    # Step 3: Preview what will be processed
    print("3. Analyzing data for SIMBAD enhancement...")
    df_full, missing_rows = enhancer.find_missing_simbad_rows(sample_csv)
    
    print("   Sample objects needing SIMBAD data:")
    for idx, row in missing_rows.head(5).iterrows():
        ra_hms = f"{int(row['ra']/15):02d}:{int((row['ra']/15 % 1)*60):02d}:{((row['ra']/15 % 1)*3600 % 60):05.2f}"
        dec_dms = f"{int(abs(row['dec'])):02d}:{int((abs(row['dec']) % 1)*60):02d}:{((abs(row['dec']) % 1)*3600 % 60):05.2f}"
        if row['dec'] < 0:
            dec_dms = "-" + dec_dms
        else:
            dec_dms = "+" + dec_dms
        print(f"   • Object {row['id']:2d}: RA {row['ra']:9.5f}° ({ra_hms}) Dec {row['dec']:8.5f}° ({dec_dms})")
    print()
    
    # Step 4: Run enhancement
    print("4. Enhancing data with SIMBAD information...")
    print("   This will:")
    print("   • Query SIMBAD for each object (3\" radius first, then 5\")")
    print("   • Extract: main_id, object_type, spectral_type, radial_velocity, distance")  
    print("   • Generate friendly source names")
    print("   • Log objects without matches")
    print("   • Cache results to avoid duplicate queries")
    print()
    
    try:
        enhanced_csv = enhancer.enhance_csv_with_simbad(sample_csv, 'sample_enhanced.csv')
        print(f"   + Enhancement completed!")
        print(f"   + Results saved to: {enhanced_csv}")
        
        # Step 5: Show results
        print()
        print("5. Results Summary:")
        enhanced_df = pd.read_csv(enhanced_csv)
        
        # Count matches
        matches = enhanced_df[enhanced_df['simbad_main_id'] != '']
        
        print(f"   • Total objects processed: {len(enhanced_df)}")
        print(f"   • SIMBAD matches found: {len(matches)}")
        print(f"   • Success rate: {len(matches)/len(enhanced_df)*100:.1f}%")
        
        if len(matches) > 0:
            print()
            print("   Sample matches found:")
            for idx, row in matches.head(3).iterrows():
                name = row['source_name'] if row['source_name'] else row['simbad_main_id']
                otype = f" ({row['otype']})" if row['otype'] else ""
                print(f"   • {name}{otype}")
                if not np.isnan(row.get('rv_value', np.nan)):
                    print(f"     Radial velocity: {row['rv_value']:.1f} km/s")
                if not np.isnan(row.get('distance_result', np.nan)):
                    print(f"     Distance: {row['distance_result']:.1f} pc")
        
        # Show misses
        if enhancer.misses:
            print()
            print(f"   • Objects without SIMBAD matches: {len(enhancer.misses)}")
            misses_file = 'sample_enhanced.misses.csv'
            if os.path.exists(misses_file):
                print(f"   • Miss details logged to: {misses_file}")
        
        # Step 6: Database integration
        print()
        print("6. Database Integration:")
        print("   The enhancer can:")
        print("   • Add 'source_name' columns to your SQLite database")
        print("   • Sync all enhanced data to database tables")
        print("   • Update existing photometry_sources records")
        print("   • Create/update simbad_matches records")
        print()
        print("   To sync to database, use:")
        print("   enhancer.sync_csv_to_database('sample_enhanced.csv')")
        
        # Step 7: Show the enhanced CSV structure
        print()
        print("7. Enhanced CSV Structure:")
        print("   New/updated columns added:")
        cols_info = [
            ('simbad_main_id', 'SIMBAD main identifier (e.g., HD 216763)'),
            ('otype', 'Object type (e.g., Star, Galaxy, Nebula)'),
            ('sp_type', 'Spectral type (e.g., G5V, K2III)'),
            ('rv_value', 'Radial velocity in km/s'),
            ('distance_result', 'Distance in parsecs'),
            ('source_name', 'Clean, readable source name')
        ]
        
        for col, desc in cols_info:
            print(f"   • {col:16s}: {desc}")
        
    except KeyboardInterrupt:
        print("\n   Enhancement interrupted by user")
    except Exception as e:
        print(f"\n   Enhancement failed: {e}")
        print("   This might be due to:")
        print("   • Network connectivity issues")
        print("   • SIMBAD service temporarily unavailable")
        print("   • Rate limiting by SIMBAD")
    
    # Cleanup
    print()
    print("8. Cleaning up...")
    cleanup_files = [sample_csv, 'sample_enhanced.csv', 'sample_enhanced.misses.csv']
    for file in cleanup_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"   + Removed {file}")
    
    print()
    print("=" * 60)
    print("Example completed!")
    print("=" * 60)
    print()
    print("To use with your actual data:")
    print()
    print("python simbad_simple.py v1.csv --batch-size 50 --delay 2.0")
    print("python simbad_simple.py v1.csv --output v1_enhanced.csv --sync-db")
    print()

if __name__ == '__main__':
    demonstrate_enhancement()