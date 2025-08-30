#!/usr/bin/env python3
"""
Debug script for astrometry.net upload issues.
Tests the API call step by step to identify the problem.
"""

import os
import sys
import json
import requests
import tempfile
import numpy as np
from astropy.io import fits

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import setup_logging

def create_minimal_test_fits():
    """Create a minimal test FITS file."""
    # Create a small test image
    image_data = np.random.randint(1000, 5000, (100, 100), dtype=np.uint16)
    
    # Create FITS file
    temp_dir = tempfile.mkdtemp(prefix='debug_astrometry_')
    fits_path = os.path.join(temp_dir, 'debug_test.fits')
    
    # Create minimal FITS HDU
    hdu = fits.PrimaryHDU(data=image_data)
    header = hdu.header
    
    # Essential headers only
    header['SIMPLE'] = True
    header['BITPIX'] = 16
    header['NAXIS'] = 2
    header['NAXIS1'] = 100
    header['NAXIS2'] = 100
    header['BZERO'] = 32768
    header['BSCALE'] = 1
    
    hdu.writeto(fits_path, overwrite=True)
    print(f"Created minimal test FITS: {fits_path}")
    print(f"File size: {os.path.getsize(fits_path)} bytes")
    
    return fits_path, temp_dir

def test_api_login():
    """Test basic API login."""
    print("\n=== Testing API Login ===")
    
    api_key = "oecnptibpffoyzrl"
    base_url = "https://nova.astrometry.net/api/"
    
    session = requests.Session()
    
    try:
        login_url = f"{base_url}login"
        login_data = {'request-json': json.dumps({"apikey": api_key})}
        
        print(f"Login URL: {login_url}")
        print(f"Login data: {login_data}")
        
        response = session.post(login_url, data=login_data, timeout=30)
        print(f"Login response status: {response.status_code}")
        print(f"Login response headers: {dict(response.headers)}")
        print(f"Login response text: {response.text}")
        
        if response.status_code == 200:
            result = response.json()
            if result.get('status') == 'success':
                session_key = result.get('session')
                print(f"[OK] Login successful: {session_key}")
                return session, session_key
            else:
                print(f"[FAIL] Login failed: {result}")
                return None, None
        else:
            print(f"[FAIL] HTTP error: {response.status_code}")
            return None, None
            
    except Exception as e:
        print(f"[FAIL] Login exception: {e}")
        return None, None

def test_original_upload_format(session, session_key, fits_path):
    """Test upload with the original format that used to work."""
    print("\n=== Testing Original Upload Format ===")
    
    upload_url = "https://nova.astrometry.net/api/upload"
    
    # Original working format
    submission_data = {
        "session": session_key,
        "allow_commercial_use": "d",
        "allow_modifications": "d", 
        "publicly_visible": "y",
        "scale_units": "degwidth",
        "scale_type": "ul",
        "scale_lower": 0.1,
        "scale_upper": 10,
    }
    
    print(f"Upload URL: {upload_url}")
    print(f"Submission data: {submission_data}")
    print(f"File: {fits_path} ({os.path.getsize(fits_path)} bytes)")
    
    try:
        with open(fits_path, 'rb') as f:
            files = {
                'file': (os.path.basename(fits_path), f, 'application/fits'),
                'request-json': (None, json.dumps(submission_data), 'application/json')
            }
            
            print("Sending upload request...")
            response = session.post(upload_url, files=files, timeout=60)
            
        print(f"Upload response status: {response.status_code}")
        print(f"Upload response headers: {dict(response.headers)}")
        print(f"Upload response text: {response.text}")
        
        if response.status_code == 200:
            result = response.json()
            print(f"Upload result: {result}")
            if result.get('status') == 'success':
                print(f"[OK] Upload successful: {result.get('subid')}")
                return True, result.get('subid')
            else:
                print(f"[FAIL] Upload failed: {result}")
                return False, None
        else:
            print(f"[FAIL] HTTP error {response.status_code}")
            return False, None
            
    except Exception as e:
        print(f"[FAIL] Upload exception: {e}")
        return False, None

def test_web_compatible_format(session, session_key, fits_path):
    """Test upload with web-compatible format."""
    print("\n=== Testing Web-Compatible Format ===")
    
    upload_url = "https://nova.astrometry.net/api/upload"
    
    # Try to mimic web form exactly
    submission_data = {
        "session": session_key,
        "allow_commercial_use": "d",
        "allow_modifications": "d", 
        "publicly_visible": "y"
    }
    
    print(f"Minimal submission data: {submission_data}")
    
    try:
        with open(fits_path, 'rb') as f:
            files = {
                'file': (os.path.basename(fits_path), f, 'application/octet-stream')
            }
            data = {
                'request-json': json.dumps(submission_data)
            }
            
            print("Sending web-compatible request...")
            response = session.post(upload_url, files=files, data=data, timeout=60)
            
        print(f"Upload response status: {response.status_code}")
        print(f"Upload response text: {response.text[:500]}...")
        
        if response.status_code == 200:
            result = response.json()
            if result.get('status') == 'success':
                print(f"[OK] Web-compatible upload successful: {result.get('subid')}")
                return True, result.get('subid')
            else:
                print(f"[FAIL] Web-compatible upload failed: {result}")
                return False, None
        else:
            print(f"[FAIL] HTTP error {response.status_code}")
            return False, None
            
    except Exception as e:
        print(f"[FAIL] Upload exception: {e}")
        return False, None

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
    """Main debug function."""
    print("Astrometry.net Upload Debug Script")
    print("=" * 50)
    
    setup_logging()
    
    temp_dir = None
    
    try:
        # Create test FITS file
        fits_path, temp_dir = create_minimal_test_fits()
        
        # Test login
        session, session_key = test_api_login()
        if not session_key:
            print("Cannot proceed without valid session key")
            return False
        
        # Test original format
        success1, subid1 = test_original_upload_format(session, session_key, fits_path)
        
        # Test web-compatible format
        success2, subid2 = test_web_compatible_format(session, session_key, fits_path)
        
        # Results
        print(f"\n=== Results ===")
        print(f"Original format: {'SUCCESS' if success1 else 'FAILED'}")
        print(f"Web-compatible format: {'SUCCESS' if success2 else 'FAILED'}")
        
        if success1 or success2:
            print("[SUCCESS] At least one format works!")
            if success1:
                print(f"Original format submission ID: {subid1}")
            if success2:
                print(f"Web-compatible format submission ID: {subid2}")
            return True
        else:
            print("[FAIL] Both formats failed")
            return False
            
    except Exception as e:
        print(f"Debug script error: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        if temp_dir:
            cleanup_test_files(temp_dir)

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)