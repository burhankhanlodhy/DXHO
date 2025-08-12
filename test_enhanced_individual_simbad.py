#!/usr/bin/env python3
"""
Test script to validate the enhanced individual SIMBAD queries with better connection error prevention.
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

from photometry import Photometry

def test_enhanced_individual_queries():
    """Test enhanced individual SIMBAD queries with better connection error prevention."""
    
    logger.info("Starting enhanced individual SIMBAD query test")
    logger.info("="*60)
    
    try:
        # Initialize photometry pipeline
        pipeline = Photometry()
        
        if pipeline.simbad_tap is None:
            logger.error("SIMBAD TAP service not available - cannot run test")
            return False
        
        # Create test coordinates - mix of bright stars that should be in SIMBAD
        test_objects = [
            ("Sirius", 101.287, -16.716),
            ("Rigel", 78.634, -8.202),
            ("Betelgeuse", 88.793, 7.407),
            ("Bellatrix", 81.283, 6.350),
            ("Aldebaran", 68.980, 16.509),
            ("Capella", 79.172, 45.998),
            ("Vega", 279.234, 38.784),
            ("Arcturus", 213.915, 19.182),
            ("Spica", 201.298, -11.161),
            ("Antares", 247.352, -26.432),
        ]
        
        logger.info(f"Testing enhanced individual queries with {len(test_objects)} bright stars")
        
        # Convert to SkyCoord objects
        sky_coords = []
        names = []
        for name, ra, dec in test_objects:
            sky_coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
            sky_coords.append(sky_coord)
            names.append(name)
        
        # Test the enhanced individual query method
        start_time = time.time()
        matches = []
        connection_errors = 0
        successful_matches = 0
        
        logger.info("Starting individual SIMBAD queries with enhanced error prevention...")
        
        for i, (coord, name) in enumerate(zip(sky_coords, names)):
            ra = coord.ra.degree
            dec = coord.dec.degree
            
            logger.info(f"Querying {name} ({i+1}/{len(test_objects)}): RA={ra:.3f}, Dec={dec:.3f}")
            
            try:
                # Use the enhanced query method
                match = pipeline._query_simbad_tap(ra, dec, radius_arcsec=5.0)
                matches.append(match)
                
                if match and match.get('main_id'):
                    successful_matches += 1
                    logger.info(f"  ✓ Found: {match['main_id']} ({match['otype']})")
                else:
                    logger.warning(f"  ✗ No match found for {name}")
                    
            except Exception as e:
                connection_errors += 1
                logger.error(f"  ✗ Query failed for {name}: {e}")
                matches.append({'main_id': '', 'otype': '', 'sp_type': ''})
            
            # Add delay between queries (this will be handled by the enhanced method)
            if i < len(test_objects) - 1:
                time.sleep(0.1)  # Small additional delay for testing
        
        elapsed_time = time.time() - start_time
        
        # Report results
        logger.info("\n" + "="*60)
        logger.info("ENHANCED INDIVIDUAL QUERY TEST RESULTS")
        logger.info("="*60)
        
        logger.info(f"Total queries: {len(test_objects)}")
        logger.info(f"Successful matches: {successful_matches}")
        logger.info(f"Connection errors: {connection_errors}")
        logger.info(f"Total time: {elapsed_time:.1f}s")
        logger.info(f"Average time per query: {elapsed_time/len(test_objects):.1f}s")
        logger.info(f"Success rate: {(successful_matches/len(test_objects)*100):.1f}%")
        logger.info(f"Error rate: {(connection_errors/len(test_objects)*100):.1f}%")
        
        # Check connection failure tracking
        if hasattr(pipeline, 'simbad_connection_failures'):
            logger.info(f"Final connection failure count: {pipeline.simbad_connection_failures}")
        
        # Detailed results
        logger.info("\nDETAILED RESULTS:")
        for i, (name, match) in enumerate(zip(names, matches)):
            if match and match.get('main_id'):
                logger.info(f"  {i+1:2d}. {name:12s}: {match['main_id']:20s} ({match['otype']:5s})")
            else:
                logger.info(f"  {i+1:2d}. {name:12s}: NO MATCH")
        
        # Assessment
        logger.info("\n" + "="*60)
        logger.info("ASSESSMENT")
        logger.info("="*60)
        
        if connection_errors == 0:
            logger.info("✓ EXCELLENT: No connection errors detected")
        elif connection_errors <= 2:
            logger.info("✓ GOOD: Minimal connection errors (≤2)")
        elif connection_errors <= 5:
            logger.warning("⚠ FAIR: Some connection errors detected (≤5)")
        else:
            logger.error("✗ POOR: Many connection errors detected (>5)")
        
        if successful_matches >= len(test_objects) * 0.9:
            logger.info("✓ EXCELLENT: >90% match success rate")
        elif successful_matches >= len(test_objects) * 0.8:
            logger.info("✓ GOOD: >80% match success rate")  
        elif successful_matches >= len(test_objects) * 0.7:
            logger.warning("⚠ FAIR: >70% match success rate")
        else:
            logger.error("✗ POOR: <70% match success rate")
        
        # Overall success criteria
        success = (connection_errors <= 2 and successful_matches >= len(test_objects) * 0.8)
        
        if success:
            logger.info("✓ OVERALL: Enhanced individual SIMBAD queries are working well")
            logger.info("  Connection error prevention measures appear to be effective")
            return True
        else:
            logger.error("✗ OVERALL: Enhanced individual SIMBAD queries need further improvement")
            logger.error("  Consider additional connection error prevention measures")
            return False
            
    except Exception as e:
        logger.error(f"Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_enhanced_individual_queries()
    sys.exit(0 if success else 1)