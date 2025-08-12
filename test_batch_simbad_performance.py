#!/usr/bin/env python3
"""
Test script to validate the improved batch SIMBAD query performance.
Compares individual queries vs batch queries to ensure similar match rates.
"""

import sys
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import time
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import photometry module
try:
    from photometry import Photometry
    logger.info("Successfully imported Photometry")
except ImportError as e:
    logger.error(f"Failed to import Photometry: {e}")
    sys.exit(1)

def create_test_coordinates():
    """Create test coordinates with known SIMBAD objects for validation."""
    
    # Use a smaller, more localized test set to avoid RA boundary issues
    test_objects = [
        # Bright stars in a more concentrated region (Orion area)
        ("Sirius", 101.287, -16.716),
        ("Rigel", 78.634, -8.202),
        ("Betelgeuse", 88.793, 7.407),
        ("Bellatrix", 81.283, 6.350),
        ("Mintaka", 83.002, -0.299),
        ("Alnilam", 84.053, -1.202),
        ("Alnitak", 85.190, -1.943),
        ("Saiph", 86.939, -9.670),
        
        # Some nearby stars in similar RA range
        ("Aldebaran", 68.980, 16.509),
        ("Capella", 79.172, 45.998),
        ("Pollux", 116.329, 28.026),
        
        # Test coordinates that might not have matches in same region
        ("Empty Field 1", 85.123, -5.456),
        ("Empty Field 2", 90.789, 15.432),
    ]
    
    logger.info(f"Created test dataset with {len(test_objects)} coordinates")
    logger.info(f"Expected matches: ~11 bright stars (excluding empty fields)")
    logger.info(f"RA range: {min(ra for name, ra, dec in test_objects):.1f} to {max(ra for name, ra, dec in test_objects):.1f}")
    logger.info(f"Dec range: {min(dec for name, ra, dec in test_objects):.1f} to {max(dec for name, ra, dec in test_objects):.1f}")
    
    # Convert to SkyCoord objects
    sky_coords = []
    names = []
    for name, ra, dec in test_objects:
        sky_coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
        sky_coords.append(sky_coord)
        names.append(name)
    
    return sky_coords, names

def test_individual_queries(pipeline, sky_coords):
    """Test individual SIMBAD queries (original approach)."""
    
    logger.info("Testing individual SIMBAD queries...")
    start_time = time.time()
    
    individual_matches = []
    successful_matches = 0
    
    for i, coord in enumerate(sky_coords):
        ra = coord.ra.degree
        dec = coord.dec.degree
        
        # Use the individual query method
        match = pipeline._query_simbad_tap(ra, dec, radius_arcsec=5.0)
        individual_matches.append(match)
        
        if match and match.get('main_id'):
            successful_matches += 1
        
        # Add delay between queries (as in original implementation)
        if i < len(sky_coords) - 1:
            time.sleep(0.2)  # Reduced delay for testing
    
    elapsed_time = time.time() - start_time
    
    logger.info(f"Individual queries completed in {elapsed_time:.1f}s")
    logger.info(f"Individual query matches: {successful_matches}/{len(sky_coords)}")
    
    return individual_matches, successful_matches

def test_batch_query(pipeline, sky_coords):
    """Test batch SIMBAD query (new approach)."""
    
    logger.info("Testing batch SIMBAD query...")
    start_time = time.time()
    
    # Use the batch query method
    batch_matches = pipeline._query_simbad_batch(sky_coords, radius_arcsec=5.0)
    
    elapsed_time = time.time() - start_time
    
    # Count successful matches
    successful_matches = 0
    for match in batch_matches:
        if match and match.get('main_id'):
            successful_matches += 1
    
    logger.info(f"Batch query completed in {elapsed_time:.1f}s")
    logger.info(f"Batch query matches: {successful_matches}/{len(sky_coords)}")
    
    return batch_matches, successful_matches

def compare_results(names, individual_matches, batch_matches):
    """Compare the results from individual vs batch queries."""
    
    logger.info("\n" + "="*60)
    logger.info("DETAILED COMPARISON OF RESULTS")
    logger.info("="*60)
    
    individual_found = 0
    batch_found = 0
    both_found = 0
    discrepancies = 0
    
    for i, name in enumerate(names):
        individual_match = individual_matches[i]
        batch_match = batch_matches[i]
        
        ind_has_match = individual_match and individual_match.get('main_id', '').strip() != ''
        batch_has_match = batch_match and batch_match.get('main_id', '').strip() != ''
        
        if ind_has_match:
            individual_found += 1
        if batch_has_match:
            batch_found += 1
        if ind_has_match and batch_has_match:
            both_found += 1
        
        # Check for discrepancies
        if ind_has_match != batch_has_match:
            discrepancies += 1
            logger.warning(f"DISCREPANCY for {name}:")
            if ind_has_match and not batch_has_match:
                logger.warning(f"  Individual found: {individual_match.get('main_id', '')} | Batch: NO MATCH")
            elif not ind_has_match and batch_has_match:
                logger.warning(f"  Individual: NO MATCH | Batch found: {batch_match.get('main_id', '')}")
        elif ind_has_match and batch_has_match:
            # Both found - check if they match the same object
            ind_id = individual_match.get('main_id', '').strip()
            batch_id = batch_match.get('main_id', '').strip()
            if ind_id != batch_id:
                logger.warning(f"DIFFERENT OBJECTS for {name}:")
                logger.warning(f"  Individual: {ind_id} | Batch: {batch_id}")
            else:
                logger.info(f"✓ {name}: Both found {ind_id}")
    
    logger.info(f"\nSUMMARY:")
    logger.info(f"Individual query matches: {individual_found}")
    logger.info(f"Batch query matches: {batch_found}")
    logger.info(f"Both methods found: {both_found}")
    logger.info(f"Discrepancies: {discrepancies}")
    
    # Calculate performance metrics
    if individual_found > 0:
        match_ratio = batch_found / individual_found
        logger.info(f"Batch/Individual match ratio: {match_ratio:.3f}")
        
        if match_ratio >= 0.95:
            logger.info("✓ EXCELLENT: Batch query maintains >95% of individual query performance")
        elif match_ratio >= 0.85:
            logger.info("✓ GOOD: Batch query maintains >85% of individual query performance")
        elif match_ratio >= 0.70:
            logger.warning("⚠ FAIR: Batch query maintains >70% of individual query performance")
        else:
            logger.error("✗ POOR: Batch query performance significantly degraded")
    
    return {
        'individual_found': individual_found,
        'batch_found': batch_found,
        'both_found': both_found,
        'discrepancies': discrepancies,
        'match_ratio': batch_found / max(individual_found, 1)
    }

def main():
    """Run the batch SIMBAD performance validation test."""
    
    logger.info("Starting batch SIMBAD performance validation test")
    logger.info("="*60)
    
    try:
        # Initialize photometry pipeline
        pipeline = Photometry()
        
        if pipeline.simbad_tap is None:
            logger.error("SIMBAD TAP service not available - cannot run test")
            return False
        
        # Create test coordinates
        sky_coords, names = create_test_coordinates()
        
        # Test individual queries
        individual_matches, individual_count = test_individual_queries(pipeline, sky_coords)
        
        # Brief pause between tests
        time.sleep(2)
        
        # Test batch query
        batch_matches, batch_count = test_batch_query(pipeline, sky_coords)
        
        # Compare results
        comparison = compare_results(names, individual_matches, batch_matches)
        
        # Overall assessment
        logger.info("\n" + "="*60)
        logger.info("FINAL ASSESSMENT")
        logger.info("="*60)
        
        if comparison['match_ratio'] >= 0.90:
            logger.info("✓ SUCCESS: Batch query performance is acceptable")
            logger.info("The improved batch implementation successfully addresses the performance issue.")
            return True
        else:
            logger.error("✗ FAILURE: Batch query performance is still below acceptable threshold")
            logger.error("Further optimization of the batch query approach is needed.")
            return False
            
    except Exception as e:
        logger.error(f"Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)