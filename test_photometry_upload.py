#!/usr/bin/env python3
"""
Test script to verify photometry.py can upload to astrometry.net correctly.
"""

import os
import sys
import tempfile
import numpy as np
from astropy.io import fits

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import setup_logging
from photometry import Photometry

def create_test_fits():
    """Create a test FITS file for upload testing."""
    # Create a small test image with some stars
    image_data = np.random.randint(1000, 2000, (200, 200), dtype=np.uint16)
    
    # Add a few bright "stars"
    star_positions = [(50, 50), (100, 150), (150, 100)]
    for x, y in star_positions:
        for dx in range(-3, 4):
            for dy in range(-3, 4):
                if 0 <= x+dx < 200 and 0 <= y+dy < 200:
                    distance = np.sqrt(dx*dx + dy*dy)
                    if distance <= 3:
                        intensity = int(2000 * np.exp(-distance*distance/2))
                        new_value = image_data[y+dy, x+dx] + intensity
                        image_data[y+dy, x+dx] = min(65535, new_value)
    
    # Create FITS file
    temp_dir = tempfile.mkdtemp(prefix='photometry_upload_test_')
    fits_path = os.path.join(temp_dir, 'test_upload.fits')
    
    # Create FITS HDU
    hdu = fits.PrimaryHDU(data=image_data)
    header = hdu.header
    
    # Add essential headers
    header['SIMPLE'] = True
    header['BITPIX'] = 16
    header['NAXIS'] = 2
    header['NAXIS1'] = 200
    header['NAXIS2'] = 200
    header['BZERO'] = 32768
    header['BSCALE'] = 1
    
    # Add some metadata
    header['OBJECT'] = 'Test Upload'
    header['ORIGIN'] = 'DEHO Test Suite'
    header['INSTRUME'] = 'Test Camera'
    
    hdu.writeto(fits_path, overwrite=True)
    print(f"Created test FITS file: {fits_path}")
    print(f"File size: {os.path.getsize(fits_path)} bytes")
    
    return fits_path, temp_dir

def test_astrometry_upload():
    """Test uploading to astrometry.net via photometry module."""
    print("\n=== Testing Photometry Module Upload ===")
    
    # Create test FITS file
    fits_path, temp_dir = create_test_fits()
    
    # Initialize photometry module
    photometry = Photometry(gui_reference=None, enable_database=False)
    
    try:
        # Test login
        print("Testing login...")
        session_key = photometry._astrometry_login()
        if not session_key:
            print("[FAIL] Login failed")
            return False, temp_dir
        
        print(f"[OK] Login successful: {session_key}")
        
        # Test upload
        print("Testing upload...")
        submission_id = photometry._submit_image_for_solving(fits_path, session_key)
        
        if submission_id:
            print(f"[OK] Upload successful: Submission ID {submission_id}")
            return True, temp_dir
        else:
            print("[FAIL] Upload failed")
            return False, temp_dir
            
    except Exception as e:
        print(f"[FAIL] Exception during test: {e}")
        import traceback
        traceback.print_exc()
        return False, temp_dir

def cleanup_test_files(temp_dir):
    """Clean up test files."""
    print(f"\n=== Cleaning up: {temp_dir} ===")
    import shutil
    try:
        shutil.rmtree(temp_dir)
        print("[OK] Cleanup successful")
    except Exception as e:
        print(f"Warning: Cleanup failed: {e}")

def main():
    """Main test function."""
    print("Photometry Module Upload Test")
    print("=" * 40)
    
    setup_logging()
    
    try:
        success, temp_dir = test_astrometry_upload()
        
        print(f"\n=== Test Result ===")
        if success:
            print("[SUCCESS] Photometry module can upload to astrometry.net!")
        else:
            print("[FAIL] Photometry module upload failed")
            
        return success
        
    except Exception as e:
        print(f"Test error: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        try:
            cleanup_test_files(temp_dir)
        except:
            pass

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)