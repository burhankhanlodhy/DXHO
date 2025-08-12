#!/usr/bin/env python3
"""
Simple test script to verify SQLite database functionality.
"""

import os
import sys
import numpy as np
from datetime import datetime
from database import ObservatoryDatabase, DatabaseConfig

def test_database_functionality():
    """Test basic database operations."""
    print("Testing SQLite database functionality...")
    
    # Initialize database
    db = ObservatoryDatabase()
    
    if not db.initialize_connection():
        print("Failed to initialize database connection")
        return False
    
    print("+ Database connection initialized")
    
    # Test connection
    if not db.test_connection():
        print("Failed to test database connection")
        return False
    
    print("+ Database connection tested successfully")
    
    # Create a test session
    session_data = {
        'session_name': 'Test Session',
        'fits_file_path': '/test/path/test.fits',
        'image_center_ra': 180.0,
        'image_center_dec': 0.0,
        'field_of_view_deg': 1.0,
        'pixel_scale_arcsec': 1.0,
        'observation_date': datetime.now().isoformat(),
        'camera_info': 'Test Camera',
        'exposure_time': 30.0,
        'iso_value': 800,
        'aperture_fnum': 2.8,
        'focal_length_mm': 50.0,
        'fwhm_pixels': 3.5,
        'detection_threshold': 5.0,
        'total_sources_detected': 100,
        'astrometry_solved': True,
        'wcs_available': True,
        'notes': 'Test session for SQLite migration'
    }
    
    session_id = db.create_session(session_data)
    if not session_id:
        print("Failed to create test session")
        return False
    
    print(f"+ Created test session with ID: {session_id}")
    
    # Create test photometry sources
    sources_data = [
        {
            'pixel_x': 100.5,
            'pixel_y': 200.3,
            'ra': 180.1,
            'dec': 0.1,
            'flux': 1000.0,
            'flux_err': 50.0,
            'background_flux': 100.0,
            'aperture_radius': 5.0,
            'fwhm': 3.2,
            'ellipticity': 0.1,
            'theta': 45.0,
            'peak_counts': 2000.0,
            'signal_to_noise': 20.0,
            'saturated': False,
            'edge_source': False
        },
        {
            'pixel_x': 300.7,
            'pixel_y': 400.9,
            'ra': 180.2,
            'dec': 0.2,
            'flux': 1500.0,
            'flux_err': 60.0,
            'background_flux': 110.0,
            'aperture_radius': 5.0,
            'fwhm': 3.8,
            'ellipticity': 0.2,
            'theta': -30.0,
            'peak_counts': 3000.0,
            'signal_to_noise': 25.0,
            'saturated': False,
            'edge_source': False
        }
    ]
    
    source_ids = db.insert_photometry_sources(session_id, sources_data)
    if len(source_ids) != len(sources_data):
        print(f"Failed to insert photometry sources. Expected {len(sources_data)}, got {len(source_ids)}")
        return False
    
    print(f"+ Inserted {len(source_ids)} photometry sources")
    
    # Test session statistics
    stats = db.get_session_statistics()
    if len(stats) == 0:
        print("Failed to get session statistics")
        return False
    
    print(f"+ Retrieved session statistics: {len(stats)} sessions")
    print(f"  Session details: {stats.iloc[0]['session_name']}, Sources: {stats.iloc[0]['sources_stored']}")
    
    # Test regional queries (using rectangular approximation)
    region_sources = db.query_sources_in_region(180.0, 0.0, 1.0)
    if len(region_sources) != len(sources_data):
        print(f"Failed regional query. Expected {len(sources_data)} sources, got {len(region_sources)}")
        return False
    
    print(f"+ Regional query returned {len(region_sources)} sources")
    
    # Test maintenance status
    maintenance_status = db.get_maintenance_status()
    print(f"+ Retrieved maintenance status: {len(maintenance_status)} tables")
    
    # Test database vacuum
    vacuum_result = db.vacuum_database()
    print(f"+ Database vacuum: {vacuum_result}")
    
    # Close connection
    db.close()
    print("+ Database connection closed")
    
    print("\n" + "=" * 50)
    print("ALL TESTS PASSED! SQLite migration successful!")
    print("=" * 50)
    
    return True

if __name__ == "__main__":
    success = test_database_functionality()
    sys.exit(0 if success else 1)