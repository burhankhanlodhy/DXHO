#!/usr/bin/env python3
"""
Debug script to investigate why batch SIMBAD query gets no matches.
"""

import sys
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import logging

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

from photometry import Photometry

def debug_batch_coordinate_matching():
    """Debug the coordinate matching in batch SIMBAD queries."""
    
    logger.info("Starting batch SIMBAD coordinate matching debug")
    
    # Initialize photometry pipeline
    pipeline = Photometry()
    
    # Test with a single bright star (Sirius)
    test_coord = SkyCoord(ra=101.287*u.degree, dec=-16.716*u.degree, frame='icrs')
    logger.info(f"Testing with Sirius: RA={test_coord.ra.degree:.3f}, Dec={test_coord.dec.degree:.3f}")
    
    # First test individual query
    individual_result = pipeline._query_simbad_tap(test_coord.ra.degree, test_coord.dec.degree)
    logger.info(f"Individual query result: {individual_result}")
    
    # Now test batch query with just this one coordinate
    batch_result = pipeline._query_simbad_batch([test_coord])
    logger.info(f"Batch query result: {batch_result[0] if batch_result else 'None'}")
    
    # Let's manually debug the field calculation
    ras = [test_coord.ra.degree]
    decs = [test_coord.dec.degree]
    radius_arcsec = 5.0
    radius_deg = radius_arcsec / 3600.0
    
    # Calculate field boundaries
    padding_factor = 1.5  # For single coordinate
    ra_min = min(ras) - radius_deg*padding_factor
    ra_max = max(ras) + radius_deg*padding_factor
    dec_min = min(decs) - radius_deg*padding_factor
    dec_max = max(decs) + radius_deg*padding_factor
    
    logger.info(f"Field boundaries:")
    logger.info(f"  RA: {ra_min:.6f} to {ra_max:.6f}")
    logger.info(f"  Dec: {dec_min:.6f} to {dec_max:.6f}")
    logger.info(f"  Source RA: {test_coord.ra.degree:.6f}, Dec: {test_coord.dec.degree:.6f}")
    
    # Test a manual SIMBAD query to see what we get in this region
    try:
        from astroquery.simbad import Simbad
        
        simbad = Simbad()
        simbad.ROW_LIMIT = 10
        
        # Query around Sirius
        result_table = simbad.query_region(test_coord, radius=10*u.arcsec)
        
        if result_table:
            logger.info(f"Manual SIMBAD region query found {len(result_table)} objects:")
            for row in result_table:
                logger.info(f"  {row['MAIN_ID']}: RA={row['RA']}, Dec={row['DEC']}")
        else:
            logger.warning("Manual SIMBAD region query found no objects")
            
    except Exception as e:
        logger.error(f"Manual SIMBAD query failed: {e}")
    
    return individual_result, batch_result

if __name__ == "__main__":
    debug_batch_coordinate_matching()