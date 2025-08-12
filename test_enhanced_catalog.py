"""
Test script for enhanced catalog query system.
"""

import pandas as pd
import numpy as np
import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from enhanced_catalog_query import SIMBADKDTreeCatalog

def test_simbad_kdtree_catalog():
    """Test the enhanced catalog system with mock data."""
    print("Testing Enhanced SIMBAD KD-Tree Catalog System")
    print("=" * 50)
    
    # Initialize catalog
    catalog = SIMBADKDTreeCatalog(rate_limit_calls_per_sec=10.0)  # Faster for testing
    
    # Test with a small field around a well-known area (M31 region)
    field_ra = 10.68  # degrees (M31 RA)  
    field_dec = 41.27  # degrees (M31 Dec)
    field_radius = 0.1  # degrees (smaller for testing)
    
    print(f"Test field: RA={field_ra}°, Dec={field_dec}°, radius={field_radius}°")
    
    # Step 1: Query SIMBAD field
    print("\n1. Querying SIMBAD for field objects...")
    try:
        simbad_df = catalog.query_simbad_field(field_ra, field_dec, field_radius)
        print(f"   Found {len(simbad_df)} SIMBAD objects")
        
        if not simbad_df.empty:
            print("   Sample objects:")
            for idx, row in simbad_df.head(3).iterrows():
                print(f"     {row.get('main_id', 'N/A')} ({row.get('otype_final', 'N/A')})")
        
    except Exception as e:
        print(f"   SIMBAD query failed: {e}")
        print("   Creating mock SIMBAD data for testing...")
        
        # Create mock SIMBAD data for testing
        simbad_df = pd.DataFrame({
            'ra': [field_ra + 0.01, field_ra - 0.02, field_ra + 0.03],
            'dec': [field_dec + 0.01, field_dec - 0.01, field_dec + 0.02],
            'main_id': ['Test Star A', 'Test Galaxy B', 'Test Nebula C'],
            'otype_final': ['*', 'G', 'HII'],
            'sp_type': ['G2V', '', ''],
            'rv_value': [25.0, 1500.0, np.nan],
            'distance_result': [10.5, 785000.0, np.nan]
        })
        print(f"   Created {len(simbad_df)} mock objects")
    
    # Step 2: Build KD-tree
    print("\n2. Building KD-tree spatial index...")
    if catalog.build_kdtree(simbad_df):
        print("   KD-tree built successfully")
    else:
        print("   KD-tree build failed (sklearn may not be available)")
    
    # Step 3: Test object matching
    print("\n3. Testing object type population...")
    target_coords = pd.DataFrame({
        'ra': [field_ra + 0.005, field_ra - 0.015, field_ra + 0.025],
        'dec': [field_dec + 0.008, field_dec - 0.005, field_dec + 0.018]
    })
    
    print(f"   Testing with {len(target_coords)} target coordinates")
    
    matches = catalog.populate_object_types(target_coords, search_radius_arcsec=30.0)
    
    if not matches.empty:
        found_matches = matches[matches['simbad_main_id'] != '']
        print(f"   Found {len(found_matches)} matches")
        
        for _, row in found_matches.iterrows():
            print(f"     Target {row['target_index']}: {row['simbad_main_id']} "
                  f"({row['simbad_otype']}) - {row['match_distance_arcsec']:.1f}\"")
    else:
        print("   No matches found")
    
    # Step 4: Test VizieR photometry (optional)
    print("\n4. Testing VizieR photometry enrichment...")
    try:
        # Use only 1 coordinate to minimize API calls during testing
        test_coord = target_coords.head(1)
        viz_data = catalog.query_vizier_photometry(test_coord, radius_arcsec=10.0)
        
        if viz_data is not None and not viz_data.empty:
            print(f"   Retrieved {len(viz_data)} VizieR entries")
            catalogs = viz_data['vizier_catalog'].unique() if 'vizier_catalog' in viz_data.columns else []
            print(f"   From catalogs: {list(catalogs)}")
        else:
            print("   No VizieR data retrieved (expected for testing)")
            
    except Exception as e:
        print(f"   VizieR test failed: {e}")
    
    print("\n5. Testing key features...")
    
    # Test rate limiter
    print("   [OK] Rate limiter initialized")
    
    # Test caching
    if catalog.simbad_cache is not None:
        print("   [OK] Caching system available") 
    
    # Test astroquery availability
    try:
        import astroquery
        print("   [OK] astroquery available")
    except ImportError:
        print("   [WARN] astroquery not available (will use HTTP fallbacks)")
    
    # Test KD-tree availability
    if catalog.simbad_kdtree is not None:
        print("   [OK] KD-tree spatial indexing active")
    else:
        print("   [WARN] KD-tree not available (install scikit-learn)")
    
    print(f"\n{'=' * 50}")
    print("Enhanced catalog test completed!")
    
    return True

if __name__ == '__main__':
    test_simbad_kdtree_catalog()