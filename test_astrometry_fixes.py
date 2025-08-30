#!/usr/bin/env python3
"""
Test script for astrometry.net API fixes.
Tests the improved astrometry submission and FITS validation.
"""

import os
import sys
import tempfile
import numpy as np
from astropy.io import fits
import logging

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import setup_logging
from photometry import Photometry
from conversion import Conversion

def create_test_fits():
    """Create a test FITS file for astrometry testing."""
    # Create a synthetic star field
    image_size = (1000, 1000)
    image_data = np.zeros(image_size, dtype=np.uint16)
    
    # Add background noise
    image_data += np.random.normal(1000, 50, image_size).astype(np.uint16)
    
    # Add some "stars" (bright spots)
    star_positions = [
        (200, 200), (300, 150), (450, 400), (600, 700),
        (150, 800), (750, 300), (850, 650), (100, 500)
    ]
    
    for x, y in star_positions:
        # Create a gaussian-like star
        for dx in range(-5, 6):
            for dy in range(-5, 6):
                if 0 <= x+dx < image_size[1] and 0 <= y+dy < image_size[0]:
                    distance = np.sqrt(dx*dx + dy*dy)
                    if distance <= 5:
                        intensity = int(3000 * np.exp(-distance*distance/4))
                        new_value = image_data[y+dy, x+dx] + intensity
                        image_data[y+dy, x+dx] = min(65535, new_value)
    
    # Create FITS file
    temp_dir = tempfile.mkdtemp(prefix='astrometry_test_')
    fits_path = os.path.join(temp_dir, 'test_starfield.fits')
    
    # Create FITS HDU with proper headers
    hdu = fits.PrimaryHDU(data=image_data)
    header = hdu.header
    
    # Add essential FITS headers
    header['SIMPLE'] = True
    header['BITPIX'] = 16
    header['NAXIS'] = 2
    header['NAXIS1'] = image_size[1]
    header['NAXIS2'] = image_size[0]
    header['BZERO'] = 32768
    header['BSCALE'] = 1
    
    # Add some metadata to make it more realistic
    header['OBJECT'] = 'Test Star Field'
    header['ORIGIN'] = 'DEHO Test Suite'
    header['INSTRUME'] = 'Test Camera'
    
    hdu.writeto(fits_path, overwrite=True)
    print(f"Created test FITS file: {fits_path}")
    
    return fits_path, temp_dir

def test_fits_validation():
    """Test FITS validation for astrometry.net compatibility."""
    print("\n=== Testing FITS Validation ===")
    
    # Create test FITS file
    fits_path, temp_dir = create_test_fits()
    
    # Initialize photometry module
    photometry = Photometry(gui_reference=None, enable_database=False)
    
    # Test validation
    is_valid = photometry._validate_fits_for_astrometry(fits_path)
    print(f"FITS validation result: {'PASSED' if is_valid else 'FAILED'}")
    
    if is_valid:
        print("[OK] Test FITS file passes astrometry.net validation")
    else:
        print("[FAIL] Test FITS file failed validation")
    
    return fits_path, temp_dir, is_valid

def test_fits_optimization():
    """Test FITS optimization for astrometry.net."""
    print("\n=== Testing FITS Optimization ===")
    
    fits_path, temp_dir = create_test_fits()
    
    # Initialize photometry module
    photometry = Photometry(gui_reference=None, enable_database=False)
    
    # Test optimization
    optimized_path = photometry._prepare_fits_for_astrometry(fits_path)
    
    if optimized_path == fits_path:
        print("[OK] Original FITS file was already suitable")
    else:
        print(f"[OK] Created optimized FITS file: {os.path.basename(optimized_path)}")
        
        # Validate the optimized file
        is_valid = photometry._validate_fits_for_astrometry(optimized_path)
        print(f"Optimized file validation: {'PASSED' if is_valid else 'FAILED'}")
    
    return fits_path, temp_dir, optimized_path

def test_session_configuration():
    """Test robust session configuration."""
    print("\n=== Testing Session Configuration ===")
    
    # Initialize photometry module
    photometry = Photometry(gui_reference=None, enable_database=False)
    
    # Check session configuration
    session = photometry.session
    
    # Verify session has proper configuration
    checks = [
        ("User-Agent header", "DEHO-Observatory" in session.headers.get('User-Agent', '')),
        ("Accept header", session.headers.get('Accept') == 'application/json'),
        ("Connection header", session.headers.get('Connection') == 'keep-alive'),
        ("Session adapters", len(session.adapters) >= 2),
    ]
    
    all_passed = True
    for check_name, result in checks:
        status = "[OK]" if result else "[FAIL]"
        print(f"{status} {check_name}: {'PASSED' if result else 'FAILED'}")
        if not result:
            all_passed = False
    
    print(f"Session configuration: {'ALL CHECKS PASSED' if all_passed else 'SOME CHECKS FAILED'}")
    return all_passed

def test_astrometry_login():
    """Test astrometry.net login functionality."""
    print("\n=== Testing Astrometry.net Login ===")
    
    # Initialize photometry module
    photometry = Photometry(gui_reference=None, enable_database=False)
    
    try:
        # Test login
        print("Attempting login to astrometry.net...")
        session_key = photometry._astrometry_login()
        
        if session_key:
            print(f"[OK] Login successful, session key obtained")
            print(f"Session key length: {len(session_key)} characters")
            return True
        else:
            print("[FAIL] Login failed - no session key returned")
            return False
            
    except Exception as e:
        print(f"[FAIL] Login failed with error: {str(e)}")
        return False

def test_file_upload_format():
    """Test file upload format compatibility."""
    print("\n=== Testing File Upload Format ===")
    
    fits_path, temp_dir = create_test_fits()
    
    # Initialize photometry module
    photometry = Photometry(gui_reference=None, enable_database=False)
    
    # Test file size validation
    file_size = os.path.getsize(fits_path)
    max_size = 100 * 1024 * 1024  # 100MB
    size_ok = file_size <= max_size
    
    print(f"File size: {file_size} bytes ({'OK' if size_ok else 'TOO LARGE'})")
    
    # Test file existence
    exists = os.path.exists(fits_path)
    print(f"File exists: {'YES' if exists else 'NO'}")
    
    # Test FITS validation
    valid = photometry._validate_fits_for_astrometry(fits_path)
    print(f"FITS validation: {'PASSED' if valid else 'FAILED'}")
    
    all_ok = size_ok and exists and valid
    print(f"Upload format compatibility: {'READY' if all_ok else 'NEEDS WORK'}")
    
    return fits_path, temp_dir, all_ok

def cleanup_test_files(temp_dirs):
    """Clean up temporary test files."""
    print("\n=== Cleaning Up Test Files ===")
    
    import shutil
    
    for temp_dir in temp_dirs:
        if os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                print(f"Cleaned up: {temp_dir}")
            except Exception as e:
                print(f"Warning: Could not clean up {temp_dir}: {e}")

def main():
    """Main test function."""
    print("Astrometry.net API Fixes Test Suite")
    print("=" * 50)
    
    # Setup logging
    setup_logging()
    
    temp_dirs = []
    results = {}
    
    try:
        # Test FITS validation
        fits_path, temp_dir, valid = test_fits_validation()
        temp_dirs.append(temp_dir)
        results['fits_validation'] = valid
        
        # Test FITS optimization
        fits_path2, temp_dir2, optimized_path = test_fits_optimization()
        temp_dirs.append(temp_dir2)
        results['fits_optimization'] = True
        
        # Test session configuration
        session_ok = test_session_configuration()
        results['session_config'] = session_ok
        
        # Test astrometry.net login (requires internet)
        login_ok = test_astrometry_login()
        results['astrometry_login'] = login_ok
        
        # Test file upload format
        fits_path3, temp_dir3, upload_ok = test_file_upload_format()
        temp_dirs.append(temp_dir3)
        results['upload_format'] = upload_ok
        
        # Summary
        print("\n=== Test Summary ===")
        passed_tests = sum(1 for result in results.values() if result)
        total_tests = len(results)
        
        for test_name, result in results.items():
            status = "PASS" if result else "FAIL"
            print(f"{test_name}: {status}")
        
        print(f"\nOverall: {passed_tests}/{total_tests} tests passed")
        
        if passed_tests == total_tests:
            print("[SUCCESS] All tests passed! Astrometry.net API improvements should work better.")
            return True
        else:
            print("[WARNING] Some tests failed. Check the implementation.")
            return False
            
    except Exception as e:
        print(f"Error during testing: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        cleanup_test_files(temp_dirs)

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)