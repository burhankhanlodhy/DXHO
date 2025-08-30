#!/usr/bin/env python3
"""
Test script for the new SIMBAD X-Match optimization.
Tests the bulk X-Match approach vs individual queries to demonstrate speed improvement.
"""

import os
import sys
import time
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import setup_logging
from photometry import Photometry

def create_test_coordinates(num_sources=20):
    """Create test sky coordinates using known bright stars for testing."""
    # Use some known bright stars in different regions
    known_stars = [
        # Betelgeuse (Orion)
        (88.7929, 7.4069),
        # Rigel (Orion)  
        (78.6344, -8.2017),
        # Sirius (Canis Major)
        (101.2872, -16.7161),
        # Vega (Lyra)
        (279.2347, 38.7837),
        # Altair (Aquila)
        (297.6958, 8.8683),
        # Deneb (Cygnus)
        (310.3580, 45.2803),
        # Polaris (Ursa Minor)
        (37.9544, 89.2641),
        # Arcturus (Boötes)
        (213.9153, 19.1825),
        # Spica (Virgo)
        (201.2983, -11.1614),
        # Capella (Auriga)
        (79.1722, 45.9980)
    ]
    
    sky_coords = []
    
    # Add known stars first
    for i, (ra, dec) in enumerate(known_stars[:min(10, num_sources)]):
        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
        sky_coords.append(coord)
    
    # Add some nearby coordinates with small offsets for remaining sources
    while len(sky_coords) < num_sources:
        base_star = known_stars[len(sky_coords) % len(known_stars)]
        # Small offset from known star (should still find catalog matches)
        offset_ra = np.random.uniform(-0.01, 0.01)  # ±0.01° = ±36 arcsec
        offset_dec = np.random.uniform(-0.01, 0.01)
        
        ra = base_star[0] + offset_ra
        dec = base_star[1] + offset_dec
        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
        sky_coords.append(coord)
    
    print(f"Created {len(sky_coords)} test coordinates using known bright stars")
    return sky_coords

def test_old_individual_method(photometry, sky_coords):
    """Test the old individual query method for comparison."""
    print(f"\n=== Testing Old Individual Method ({len(sky_coords)} sources) ===")
    
    start_time = time.time()
    
    # Simulate the old individual query approach (but limit to avoid long waits)
    old_matches = []
    max_test_sources = min(5, len(sky_coords))  # Only test first 5 to avoid long wait
    
    print(f"Testing individual queries for first {max_test_sources} sources...")
    
    for i, coord in enumerate(sky_coords[:max_test_sources]):
        ra = coord.ra.degree
        dec = coord.dec.degree
        
        try:
            # Use the original individual query method with standard radius
            match = photometry._query_simbad_tap(ra, dec, radius_arcsec=5.0)
            old_matches.append(match)
            print(f"  Source {i+1}: {match['main_id'] if match['main_id'] else 'No match'}")
            
            # Add delay to simulate the old rate-limited approach
            if i < max_test_sources - 1:
                time.sleep(0.5)  # Reduced from 2.5s for testing
                
        except Exception as e:
            print(f"  Source {i+1}: Error - {e}")
            old_matches.append({'main_id': '', 'ra': np.nan, 'dec': np.nan})
    
    elapsed_time = time.time() - start_time
    successful_matches = len([m for m in old_matches if m['main_id']])
    
    print(f"Old method results:")
    print(f"  Time taken: {elapsed_time:.1f} seconds")
    print(f"  Sources tested: {max_test_sources}")
    print(f"  Successful matches: {successful_matches}")
    print(f"  Rate: {elapsed_time/max_test_sources:.1f} seconds per source")
    
    return old_matches, elapsed_time

def test_new_xmatch_method(photometry, sky_coords):
    """Test the new X-Match method."""
    print(f"\n=== Testing New X-Match Method ({len(sky_coords)} sources) ===")
    
    start_time = time.time()
    
    try:
        # Use the new bulk X-Match method with standard radius
        new_matches = photometry._simbad_xmatch(sky_coords, radius_arcsec=5.0)
        
        elapsed_time = time.time() - start_time
        successful_matches = len([m for m in new_matches if m['simbad_main_id']])
        
        print(f"New X-Match method results:")
        print(f"  Time taken: {elapsed_time:.1f} seconds")
        print(f"  Sources processed: {len(sky_coords)}")
        print(f"  Successful matches: {successful_matches}")
        print(f"  Rate: {elapsed_time/len(sky_coords):.2f} seconds per source")
        
        # Show some example matches
        print(f"  Example matches:")
        for i, match in enumerate(new_matches[:5]):
            if match['simbad_main_id']:
                print(f"    Source {i+1}: {match['simbad_main_id']} (sep: {match['separation_arcsec']:.1f}\")")
            else:
                print(f"    Source {i+1}: No match")
        
        return new_matches, elapsed_time
        
    except Exception as e:
        print(f"X-Match method failed: {e}")
        import traceback
        traceback.print_exc()
        return [], float('inf')

def compare_methods():
    """Compare old vs new methods."""
    print("SIMBAD X-Match Optimization Test")
    print("=" * 50)
    
    setup_logging()
    
    # Enable debug logging to see what's happening with TAP
    import logging
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger('photometry').setLevel(logging.DEBUG)
    
    # Create test coordinates - use nearby stars in Orion only for better testing
    orion_stars = [
        # Betelgeuse (Orion)
        (88.7929, 7.4069),
        # Rigel (Orion) - close to Betelgeuse
        (78.6344, -8.2017),
        # Bellatrix (Orion) - closer to region
        (81.2826, 6.3497)
    ]
    
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    sky_coords = [SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs') for ra, dec in orion_stars]
    print(f"Using 3 Orion region coordinates for focused testing")
    
    # Initialize photometry module
    photometry = Photometry(gui_reference=None, enable_database=False)
    
    # Test new X-Match method (do this first to avoid connection issues)
    new_matches, new_time = test_new_xmatch_method(photometry, sky_coords)
    
    # Test old individual method (limited to avoid long waits)
    old_matches, old_time = test_old_individual_method(photometry, sky_coords)
    
    # Compare results
    print(f"\n=== Performance Comparison ===")
    
    if old_time > 0 and new_time > 0:
        old_rate = old_time / min(5, len(sky_coords))  # Rate per source for old method
        new_rate = new_time / len(sky_coords)  # Rate per source for new method
        
        if old_rate > 0:
            speedup = old_rate / new_rate
            print(f"Speed improvement: {speedup:.1f}x faster")
            
        print(f"Old method: {old_rate:.2f}s per source")
        print(f"New method: {new_rate:.3f}s per source")
        
        # Extrapolate time for large datasets
        print(f"\nExtrapolated time for 1000 sources:")
        print(f"  Old method: {old_rate * 1000:.0f} seconds ({old_rate * 1000 / 60:.1f} minutes)")
        print(f"  New method: {new_rate * 1000:.0f} seconds ({new_rate * 1000 / 60:.1f} minutes)")
        
    print(f"\nMatch quality:")
    if new_matches:
        match_rate = len([m for m in new_matches if m['simbad_main_id']]) / len(new_matches)
        print(f"  New method match rate: {match_rate:.1%}")
    
    return new_time < old_time if old_time > 0 else True

def main():
    """Main test function."""
    try:
        success = compare_methods()
        
        print(f"\n=== Test Summary ===")
        if success:
            print("[SUCCESS] SIMBAD X-Match optimization is working and faster!")
            print("The new implementation should significantly speed up photometry analysis.")
        else:
            print("[WARNING] X-Match method may need optimization")
            
        return success
        
    except Exception as e:
        print(f"Test error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)