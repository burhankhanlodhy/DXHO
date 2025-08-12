#!/usr/bin/env python3
"""
Test script for catalog enhancement services
Tests SIMBAD, VizieR, and Gaia services with fallback strategies
"""

import pandas as pd
import numpy as np
from catalog_enhancement import CatalogEnhancer

def test_catalog_services():
    """Test catalog services with known coordinates."""
    
    print("=" * 60)
    print("Catalog Enhancement Service Test")
    print("=" * 60)
    print()
    
    # Initialize enhancer
    enhancer = CatalogEnhancer()
    
    # Test coordinates (from your CSV) - these should have data
    test_coords = [
        {'ra': 341.80877, 'dec': 57.181023, 'name': 'Object 1'},
        {'ra': 340.00722, 'dec': 57.295166, 'name': 'Object 2 (has ZTF data)'},
        {'ra': 340.22059, 'dec': 57.283001, 'name': 'Object 3'},
    ]
    
    print("Testing catalog services with known coordinates:")
    print()
    
    for coord in test_coords:
        ra, dec = coord['ra'], coord['dec']
        name = coord['name']
        
        print(f"Testing {name} (RA={ra:.5f}, Dec={dec:.5f}):")
        
        # Test SIMBAD services (VizieR PRIMARY, NO SIMBAD TAP)
        print("  SIMBAD Services (via VizieR - NO SIMBAD TAP):")
        
        # Try VizieR PRIMARY (astroquery)
        try:
            vizier_result = enhancer.query_vizier_primary(ra, dec, 5.0)
            if vizier_result:
                print(f"    ✓ VizieR PRIMARY: {vizier_result['main_id']} ({vizier_result['otype']})")
            else:
                print("    - VizieR PRIMARY: No match")
        except Exception as e:
            print(f"    ✗ VizieR PRIMARY: Failed ({e})")
        
        # Try VizieR HTTP fallback
        try:
            vizier_http_result = enhancer.query_vizier_http_fallback(ra, dec, 5.0)
            if vizier_http_result:
                print(f"    ✓ VizieR HTTP: {vizier_http_result['main_id']} ({vizier_http_result['otype']})")
            else:
                print("    - VizieR HTTP: No match")
        except Exception as e:
            print(f"    ✗ VizieR HTTP: Failed ({e})")
        
        # Try SIMBAD HTTP (last resort only)
        try:
            http_result = enhancer.query_simbad_http_fallback(ra, dec, 5.0)
            if http_result:
                print(f"    ✓ SIMBAD HTTP (last resort): {http_result['main_id']} ({http_result['otype']})")
            else:
                print("    - SIMBAD HTTP: No match")
        except Exception as e:
            print(f"    ✗ SIMBAD HTTP: Failed ({e})")
        
        # Test unified SIMBAD query (VizieR PRIMARY with fallbacks)
        try:
            unified_result = enhancer.query_simbad_data(ra, dec, 5.0)
            if unified_result:
                print(f"    ✓ UNIFIED (VizieR→HTTP): {unified_result['main_id']} ({unified_result['otype']})")
                source_name = enhancer.get_source_name(unified_result)
                if source_name:
                    print(f"      Source name: {source_name}")
            else:
                print("    - UNIFIED: No match found")
        except Exception as e:
            print(f"    ✗ UNIFIED: Failed ({e})")
        
        # Test Gaia services
        print("  Gaia Services:")
        try:
            gaia_result = enhancer.query_gaia_data(ra, dec, 5.0)
            if gaia_result:
                gaia_id = gaia_result.get('gaia_source_id', 'N/A')
                g_mag = gaia_result.get('phot_g_mean_mag', np.nan)
                print(f"    ✓ Gaia: {gaia_id}")
                if not np.isnan(g_mag):
                    print(f"      G magnitude: {g_mag:.2f}")
            else:
                print("    - Gaia: No match")
        except Exception as e:
            print(f"    ✗ Gaia: Failed ({e})")
        
        print()
    
    print("=" * 60)
    print("Service Availability Summary:")
    print("=" * 60)
    
    # Check service availability
    services_status = {
        'astroquery': 'Available' if enhancer.ASTROQUERY_AVAILABLE else 'Not Available',
        'requests': 'Available' if enhancer.REQUESTS_AVAILABLE else 'Not Available',
        'astropy': 'Available' if enhancer.ASTROPY_AVAILABLE else 'Not Available'
    }
    
    for service, status in services_status.items():
        print(f"{service:12s}: {status}")
    
    print()
    print("NEW Service Strategy (NO SIMBAD TAP):")
    print("1. VizieR PRIMARY (astroquery VizieR SIMBAD catalog)")
    print("2. Fall back to VizieR HTTP requests")  
    print("3. Fall back to direct SIMBAD HTTP (last resort)")
    print("4. Gaia queries use astroquery Gaia DR3")
    print()
    print("✓ SIMBAD TAP DISABLED - using reliable VizieR service!")
    print("✓ This provides maximum stability and reliability!")

def create_sample_csv_for_enhancement():
    """Create a small sample CSV to demonstrate the enhancement process."""
    
    print("=" * 60)
    print("Sample CSV Enhancement Demo")
    print("=" * 60)
    print()
    
    # Create sample data with missing catalog info
    sample_data = {
        'filename': ['test.fits'] * 10,
        'id': range(1, 11),
        'ra': [341.80877, 342.31682, 339.99572, 340.00722, 340.03693,
               340.05423, 340.22059, 340.35298, 340.54525, 340.72633],
        'dec': [57.181023, 57.143944, 57.295562, 57.295166, 57.293586,
                57.292723, 57.283001, 57.274634, 57.263864, 57.252521],
        'flux': np.random.exponential(1000, 10),
        'mag_calibrated': np.random.normal(16, 2, 10),
        
        # Missing Gaia data (empty)
        'gaia_source_id': [''] * 10,
        'gaia_ra': [np.nan] * 10,
        'gaia_dec': [np.nan] * 10,
        'phot_g_mean_mag': [np.nan] * 10,
        'phot_bp_mean_mag': [np.nan] * 10,
        'phot_rp_mean_mag': [np.nan] * 10,
        'bp_rp': [np.nan] * 10,
        'parallax': [np.nan] * 10,
        'pmra': [np.nan] * 10,
        'pmdec': [np.nan] * 10,
        
        # Missing SIMBAD data (mostly empty)
        'simbad_main_id': [''] * 10,
        'otype': [''] * 10,
        'sp_type': [''] * 10,
        'rv_value': [np.nan] * 10,
        'distance_result': [np.nan] * 10
    }
    
    # Add one existing SIMBAD entry to show it's preserved
    sample_data['simbad_main_id'][3] = 'ZTF J224001.61+571742.9'
    sample_data['otype'][3] = 'EB*'
    
    df = pd.DataFrame(sample_data)
    sample_csv = 'sample_catalog_test.csv'
    df.to_csv(sample_csv, index=False)
    
    print(f"Created sample CSV: {sample_csv}")
    print(f"  - 10 objects with coordinates")
    print(f"  - Missing Gaia data for all objects")
    print(f"  - Missing SIMBAD data for 9/10 objects")
    print()
    
    print("To enhance this file, run:")
    print(f"python catalog_enhancement.py {sample_csv} --batch-size 5 --delay 1.0")
    print()
    print("This will:")
    print("- Query SIMBAD for missing identifications")
    print("- Query Gaia DR3 for missing photometry")
    print("- Use multiple fallback services for reliability")
    print("- Create enhanced CSV with all available data")
    

if __name__ == '__main__':
    import sys
    
    if len(sys.argv) > 1 and sys.argv[1] == '--create-sample':
        create_sample_csv_for_enhancement()
    else:
        test_catalog_services()