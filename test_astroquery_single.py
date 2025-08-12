"""
Test script for SIMBAD enhancement using astroquery on a single coordinate
"""
from simbad_enhancement_astroquery import SIMBADEnhancer

def test_single_coordinate():
    """Test SIMBAD query on a single coordinate that should have a match."""
    enhancer = SIMBADEnhancer()
    
    # Test with coordinates from your CSV that have known SIMBAD matches
    # Using row 4 which already has ZTF J224001.61+571742.9
    ra = 340.00722
    dec = 57.295166
    
    print(f"Testing SIMBAD query for RA={ra}, Dec={dec}")
    
    # Try 3 arcsec search first
    result = enhancer.query_simbad_cone(ra, dec, 3.0)
    if result:
        print("3 arcsec search SUCCESS:")
        for key, value in result.items():
            print(f"  {key}: {value}")
    else:
        print("3 arcsec search: No match")
        
        # Try 5 arcsec search
        result = enhancer.query_simbad_cone(ra, dec, 5.0)
        if result:
            print("5 arcsec search SUCCESS:")
            for key, value in result.items():
                print(f"  {key}: {value}")
        else:
            print("5 arcsec search: No match")
    
    # Test source name extraction
    if result:
        source_name = enhancer.get_source_name(result)
        print(f"Source name: '{source_name}'")

if __name__ == '__main__':
    test_single_coordinate()