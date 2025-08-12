#!/usr/bin/env python3
"""
Simple test of VizieR-primary catalog enhancement
"""

from catalog_enhancement import CatalogEnhancer
import pandas as pd

def test_vizier_primary():
    """Test VizieR as primary service."""
    
    print("=" * 50)
    print("VizieR Primary Service Test")  
    print("=" * 50)
    print()
    
    # Initialize enhancer
    enhancer = CatalogEnhancer()
    
    # Test with a known bright star that should be in SIMBAD
    # Vega coordinates
    ra, dec = 279.234735, 38.783689
    print(f"Testing with Vega coordinates: RA={ra}, Dec={dec}")
    print()
    
    print("1. Testing VizieR Primary (astroquery)...")
    try:
        result = enhancer.query_vizier_primary(ra, dec, 5.0)
        if result:
            print(f"   SUCCESS: Found {result['main_id']}")
            print(f"   Object type: {result['otype']}")
            print(f"   Spectral type: {result['sp_type']}")
        else:
            print("   No match found")
    except Exception as e:
        print(f"   FAILED: {e}")
    
    print()
    print("2. Testing VizieR HTTP fallback...")
    try:
        result = enhancer.query_vizier_http_fallback(ra, dec, 5.0)
        if result:
            print(f"   SUCCESS: Found {result['main_id']}")
            print(f"   Object type: {result['otype']}")
        else:
            print("   No match found")
    except Exception as e:
        print(f"   FAILED: {e}")
    
    print()
    print("3. Testing Unified Query (VizieR B catalogs ONLY)...")
    try:
        result = enhancer.query_catalog_data(ra, dec, 5.0)
        if result:
            print(f"   SUCCESS: Found {result['main_id']}")
            print(f"   Object type: {result['otype']}")
            print(f"   Catalog source: {result.get('catalog_source', 'Unknown')}")
            source_name = enhancer.get_source_name(result)
            print(f"   Source name: {source_name}")
        else:
            print("   No match found")
    except Exception as e:
        print(f"   FAILED: {e}")
    
    print()
    print("4. Testing Gaia Query...")
    try:
        result = enhancer.query_gaia_data(ra, dec, 5.0)
        if result:
            print(f"   SUCCESS: Gaia ID {result['gaia_source_id']}")
            print(f"   G magnitude: {result.get('phot_g_mean_mag', 'N/A')}")
        else:
            print("   No match found")
    except Exception as e:
        print(f"   FAILED: {e}")
    
    print()
    print("=" * 50)
    print("Service Configuration:")
    print(f"VizieR B Catalogs: {'Available' if hasattr(enhancer, 'vizier') and enhancer.vizier else 'Not available'}")
    print(f"SIMBAD: COMPLETELY DISABLED")
    print(f"VizieR HTTP Fallbacks: Available")
    print("Catalogs: Gaia EDR3/DR3, Hipparcos, MK Types, GCVS, ASCC")
    print("=" * 50)

def create_test_csv():
    """Create a small test CSV."""
    data = {
        'ra': [279.234735, 341.80877, 340.00722],  # Vega, and your objects
        'dec': [38.783689, 57.181023, 57.295166],
        'id': [1, 2, 3],
        'gaia_source_id': ['', '', ''],
        'simbad_main_id': ['', '', ''],
        'otype': ['', '', ''],
        'source_name': ['', '', '']
    }
    
    df = pd.DataFrame(data)
    test_file = 'test_vizier.csv'
    df.to_csv(test_file, index=False)
    
    print(f"Created test file: {test_file}")
    print("To enhance:")
    print(f"python catalog_enhancement.py {test_file} --batch-size 3 --delay 1.0")
    
    return test_file

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == '--create-csv':
        create_test_csv()
    else:
        test_vizier_primary()