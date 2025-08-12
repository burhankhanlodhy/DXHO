"""
Demo script for SIMBAD enhancement functionality
Shows how the system works with mock data when SIMBAD API is not available
"""

import pandas as pd
import numpy as np
import os
from typing import Optional, Dict
from simbad_enhancement_astroquery import SIMBADEnhancer

def create_test_csv():
    """Create a small test CSV with some missing SIMBAD data."""
    test_data = {
        'filename': ['test.fits'] * 10,
        'id': range(1, 11),
        'ra': [341.80877, 342.31682, 339.99572, 340.00722, 340.03693, 
               340.05423, 340.22059, 340.35298, 340.54525, 340.72633],
        'dec': [57.181023, 57.143944, 57.295562, 57.295166, 57.293586,
                57.292723, 57.283001, 57.274634, 57.263864, 57.252521],
        'flux': [1000, 1500, 800, 2000, 1200, 1800, 900, 1100, 1300, 1600],
        'mag_calibrated': [15.5, 15.2, 15.8, 14.9, 15.4, 14.8, 15.7, 15.6, 15.3, 15.1],
        'gaia_source_id': ['2007132539214702976', '', '', '2007275853684844416', '',
                          '2007275681886171520', '2007277743470543872', '', '', ''],
        'simbad_main_id': ['', '', '', 'ZTF J224001.61+571742.9', '',
                          '', '', '', '', ''],
        'otype': ['', '', '', 'EB*', '', '', '', '', '', ''],
        'sp_type': ['', '', '', '', '', '', '', '', '', ''],
        'rv_value': [np.nan, np.nan, np.nan, np.nan, np.nan,
                     np.nan, np.nan, np.nan, np.nan, np.nan],
        'distance_result': [np.nan, np.nan, np.nan, np.nan, np.nan,
                           np.nan, np.nan, np.nan, np.nan, np.nan]
    }
    
    df = pd.DataFrame(test_data)
    return df

class MockSIMBADEnhancer(SIMBADEnhancer):
    """Mock version of SIMBAD enhancer for demo purposes."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # Mock database of known objects
        self.mock_simbad_data = {
            # RA/Dec ranges for mock objects
            (341.8, 57.18): {
                'main_id': 'HD 216763',
                'otype': 'Star',
                'sp_type': 'G5V',
                'rv_value': -12.5,
                'distance_result': 45.2,
                'simbad_ra': 341.80876,
                'simbad_dec': 57.181024,
                'match_distance_arcsec': 0.5
            },
            (342.3, 57.14): {
                'main_id': 'TYC 3956-1234-1',
                'otype': 'Star',
                'sp_type': 'K2V',
                'rv_value': 8.3,
                'distance_result': 78.9,
                'simbad_ra': 342.31680,
                'simbad_dec': 57.143945,
                'match_distance_arcsec': 1.2
            },
            (340.22, 57.28): {
                'main_id': '2MASS J22404893+5716588',
                'otype': 'Star',
                'sp_type': 'M3V',
                'rv_value': np.nan,
                'distance_result': 156.3,
                'simbad_ra': 340.22058,
                'simbad_dec': 57.283002,
                'match_distance_arcsec': 0.8
            }
        }
    
    def query_simbad_cone(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """Mock SIMBAD query that uses local data."""
        
        # Check cache first
        cached_result = self.cache.get(ra, dec, radius)
        if cached_result is not None:
            return cached_result
        
        # Look for nearby mock objects
        for (mock_ra, mock_dec), data in self.mock_simbad_data.items():
            # Simple distance check (not perfect spherical geometry but good for demo)
            ra_diff = abs(ra - mock_ra)
            dec_diff = abs(dec - mock_dec)
            
            # Convert to approximate arcseconds
            distance_arcsec = np.sqrt((ra_diff * 3600 * np.cos(np.radians(dec)))**2 + 
                                    (dec_diff * 3600)**2)
            
            if distance_arcsec <= radius:
                result = data.copy()
                result['match_distance_arcsec'] = distance_arcsec
                
                # Cache and return
                self.cache.set(ra, dec, radius, result)
                return result
        
        # No match found
        self.cache.set(ra, dec, radius, None)
        return None

def demo_simbad_enhancement():
    """Demonstrate the SIMBAD enhancement process."""
    
    print("=== SIMBAD Enhancement Demo ===")
    print()
    
    # Create test data
    print("1. Creating test CSV with missing SIMBAD data...")
    df = create_test_csv()
    test_csv = 'test_data.csv'
    df.to_csv(test_csv, index=False)
    print(f"   Created {test_csv} with {len(df)} rows")
    
    # Show initial state
    missing_count = (df['simbad_main_id'].isna() | (df['simbad_main_id'] == '')).sum()
    print(f"   {missing_count} rows have missing SIMBAD data")
    print()
    
    # Initialize mock enhancer
    print("2. Initializing SIMBAD enhancer...")
    enhancer = MockSIMBADEnhancer()
    enhancer.batch_size = 5  # Small batches for demo
    enhancer.delay_between_batches = 0.5  # Faster for demo
    print("   Using mock SIMBAD data for demo purposes")
    print()
    
    # Show what we'll be looking for
    print("3. Preview of rows needing enhancement:")
    df_preview, missing_rows = enhancer.find_missing_simbad_rows(test_csv)
    for idx, row in missing_rows.head(5).iterrows():
        print(f"   Row {idx}: RA={row['ra']:.5f}, Dec={row['dec']:.5f}")
    print()
    
    # Run enhancement
    print("4. Running SIMBAD enhancement...")
    output_csv = enhancer.enhance_csv_with_simbad(test_csv, 'test_data_enhanced.csv')
    print(f"   Enhanced data saved to: {output_csv}")
    print()
    
    # Show results
    print("5. Results summary:")
    enhanced_df = pd.read_csv(output_csv)
    
    # Count successful matches
    matches = enhanced_df[
        enhanced_df['simbad_main_id'].notna() & 
        (enhanced_df['simbad_main_id'] != '')
    ]
    
    print(f"   Total objects: {len(enhanced_df)}")
    print(f"   SIMBAD matches found: {len(matches)}")
    print(f"   Objects with source names: {(enhanced_df['source_name'] != '').sum()}")
    
    if len(matches) > 0:
        print("\n   Sample matches:")
        for idx, row in matches.iterrows():
            print(f"     {row['simbad_main_id']} ({row['otype']}) - {row['source_name']}")
    
    # Show misses
    if enhancer.misses:
        print(f"\n   Objects without SIMBAD matches: {len(enhancer.misses)}")
        misses_file = 'test_data_enhanced.misses.csv'
        if os.path.exists(misses_file):
            print(f"   Miss details saved to: {misses_file}")
    
    print("\n6. Database schema update capability:")
    print("   The enhancer can update your SQLite database schema to include:")
    print("   - source_name column in photometry_sources table")
    print("   - source_name column in simbad_matches table")
    print("   - Sync all enhanced CSV data back to database")
    
    # Show CSV structure
    print(f"\n7. Enhanced CSV structure:")
    print("   New/updated columns:")
    print("   - simbad_main_id: SIMBAD main identifier")
    print("   - otype: Object type")
    print("   - sp_type: Spectral type")  
    print("   - rv_value: Radial velocity (km/s)")
    print("   - distance_result: Distance (pc)")
    print("   - source_name: Cleaned source name")
    
    # Cleanup
    print(f"\n8. Cleaning up demo files...")
    for file in [test_csv, output_csv, 'test_data_enhanced.misses.csv']:
        if os.path.exists(file):
            os.remove(file)
            print(f"   Removed {file}")
    
    print("\n=== Demo Complete ===")
    print("\nTo use with real data:")
    print("python simbad_enhancement_astroquery.py your_file.csv --batch-size 50 --delay 2.0 --sync-db")

if __name__ == '__main__':
    demo_simbad_enhancement()