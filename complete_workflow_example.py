"""
Complete Example: SIMBAD-First Catalog Query with KD-Tree and VizieR Enhancement
================================================================================

This example demonstrates the complete workflow:
1. Query SIMBAD for a field of objects
2. Build KD-tree spatial index on SIMBAD coordinates  
3. Populate object types from SIMBAD OTYPE/OTYPE_S
4. Optionally enrich with VizieR photometry and cross-identifications

Key Features:
- Rate-limited to prevent SIMBAD API violations
- KD-tree spatial indexing for efficient coordinate matching
- Proper handling of OTYPE_S (preferred) vs OTYPE fields
- VizieR table stacking with dtype conflict resolution
- String conversion for problematic ID columns (URAT1, PS1, etc.)
"""

import pandas as pd
import numpy as np
import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from enhanced_catalog_query import SIMBADKDTreeCatalog

def complete_workflow_example():
    """Complete workflow example showing all features."""
    
    print("Enhanced Catalog Query System - Complete Workflow")
    print("=" * 60)
    
    # Initialize the enhanced catalog system
    # Rate limit to 5 calls/sec to prevent SIMBAD blocking
    catalog = SIMBADKDTreeCatalog(rate_limit_calls_per_sec=5.0)
    
    # Define your field of interest
    # Example: Small field around a star-forming region
    field_ra = 83.633  # degrees (Orion Nebula region) 
    field_dec = -5.391  # degrees
    field_radius = 0.2  # degrees (~12 arcmin radius)
    
    print(f"\nTarget Field: RA={field_ra}°, Dec={field_dec}°, radius={field_radius}° (~{field_radius*60:.1f} arcmin)")
    
    # STEP 1: Query SIMBAD for the field
    print(f"\n{'='*60}")
    print("STEP 1: Querying SIMBAD for field objects")
    print("="*60)
    
    simbad_catalog = catalog.query_simbad_field(field_ra, field_dec, field_radius)
    
    if simbad_catalog.empty:
        print("No SIMBAD objects found - trying larger field...")
        simbad_catalog = catalog.query_simbad_field(field_ra, field_dec, field_radius * 2)
    
    if not simbad_catalog.empty:
        print(f"SUCCESS: Found {len(simbad_catalog)} SIMBAD objects")
        
        # Show sample of object types
        object_types = simbad_catalog['otype_final'].value_counts()
        print(f"Object types found: {dict(object_types.head())}")
        
        # Show some examples
        print("\nSample objects:")
        for _, obj in simbad_catalog.head(5).iterrows():
            print(f"  {obj['main_id']} ({obj['otype_final']}) at {obj['ra']:.4f}, {obj['dec']:.4f}")
    else:
        print("WARNING: No SIMBAD objects found in field")
        print("Creating mock data for demonstration...")
        
        # Create realistic mock data for demonstration
        simbad_catalog = pd.DataFrame({
            'ra': np.random.uniform(field_ra - field_radius, field_ra + field_radius, 20),
            'dec': np.random.uniform(field_dec - field_radius, field_dec + field_radius, 20),
            'main_id': [f'HD {12000 + i}' for i in range(20)],
            'otype_final': np.random.choice(['*', 'HII', 'Em*', 'YSO'], 20),
            'sp_type': [''] * 20,
            'rv_value': np.random.normal(20, 10, 20),
            'distance_result': np.random.uniform(400, 600, 20)  # pc
        })
        print(f"Created {len(simbad_catalog)} mock objects for testing")
    
    # STEP 2: Build KD-tree spatial index 
    print(f"\n{'='*60}")
    print("STEP 2: Building KD-tree spatial index")  
    print("="*60)
    
    if catalog.build_kdtree(simbad_catalog):
        print(f"SUCCESS: KD-tree built with {len(simbad_catalog)} objects")
        print("Spatial indexing active for fast coordinate matching")
    else:
        print("WARNING: KD-tree not available (install scikit-learn for optimal performance)")
        print("Will use slower linear search fallback")
    
    # STEP 3: Test object type population for target coordinates
    print(f"\n{'='*60}")
    print("STEP 3: Object type population for target coordinates")
    print("="*60)
    
    # Simulate some target coordinates (e.g., from photometry)
    target_sources = pd.DataFrame({
        'source_id': range(1, 11),
        'ra': np.random.uniform(field_ra - field_radius/2, field_ra + field_radius/2, 10),
        'dec': np.random.uniform(field_dec - field_radius/2, field_dec + field_radius/2, 10),
        'mag_v': np.random.uniform(10, 16, 10)  # V magnitude
    })
    
    print(f"Matching {len(target_sources)} target coordinates against SIMBAD catalog")
    
    # Populate object types using SIMBAD data
    matches = catalog.populate_object_types(target_sources, search_radius_arcsec=10.0)
    
    # Analyze results  
    successful_matches = matches[matches['simbad_main_id'] != '']
    print(f"\nMatching Results:")
    print(f"  Total targets: {len(matches)}")
    print(f"  Successful matches: {len(successful_matches)}")
    print(f"  Match rate: {len(successful_matches)/len(matches)*100:.1f}%")
    
    if not successful_matches.empty:
        print(f"\nSample matches:")
        for _, match in successful_matches.head(3).iterrows():
            print(f"  Source {match['target_index']}: {match['simbad_main_id']} "
                  f"({match['simbad_otype']}) - {match['match_distance_arcsec']:.2f}\"")
    
    # STEP 4: VizieR photometry enrichment
    print(f"\n{'='*60}")
    print("STEP 4: VizieR photometry and cross-identification enrichment")
    print("="*60)
    
    print("Querying VizieR for additional photometry and cross-IDs...")
    print("(This handles dtype conflicts and ID column conversion automatically)")
    
    # Use a subset for VizieR to minimize API calls
    vizier_targets = target_sources.head(3)  # Test with 3 objects
    
    vizier_data = catalog.query_vizier_photometry(vizier_targets, radius_arcsec=5.0)
    
    if vizier_data is not None and not vizier_data.empty:
        print(f"\nVizieR Enrichment Results:")
        print(f"  Total entries: {len(vizier_data)}")
        
        catalogs = vizier_data['vizier_catalog'].unique() if 'vizier_catalog' in vizier_data.columns else []
        print(f"  Catalogs accessed: {list(catalogs)}")
        
        # Show available photometric columns
        mag_cols = [col for col in vizier_data.columns if 'mag' in col.lower() or 'phot' in col.lower()]
        if mag_cols:
            print(f"  Photometric columns: {mag_cols}")
        
        # Show ID columns (demonstrating string conversion)
        id_cols = [col for col in vizier_data.columns if col in ['URAT1', 'PS1', 'Source', '_2MASS', 'AllWISE']]
        if id_cols:
            print(f"  Cross-ID columns (converted to str): {id_cols}")
    else:
        print("  No VizieR data retrieved (expected with limited API calls)")
    
    # STEP 5: Demonstrate combined results
    print(f"\n{'='*60}")
    print("STEP 5: Combined catalog results")
    print("="*60)
    
    # Merge SIMBAD matches with original targets
    enhanced_catalog = target_sources.copy()
    
    # Add SIMBAD information
    for _, match in matches.iterrows():
        idx = match['target_index']
        if idx < len(enhanced_catalog):
            enhanced_catalog.loc[idx, 'simbad_id'] = match['simbad_main_id']
            enhanced_catalog.loc[idx, 'object_type'] = match['simbad_otype'] 
            enhanced_catalog.loc[idx, 'spectral_type'] = match['simbad_sp_type']
            enhanced_catalog.loc[idx, 'match_distance'] = match['match_distance_arcsec']
    
    # Fill empty SIMBAD fields
    enhanced_catalog['simbad_id'] = enhanced_catalog['simbad_id'].fillna('')
    enhanced_catalog['object_type'] = enhanced_catalog['object_type'].fillna('')
    enhanced_catalog['spectral_type'] = enhanced_catalog['spectral_type'].fillna('')
    
    print("Enhanced catalog with SIMBAD object types:")
    print("=" * 80)
    cols_to_show = ['source_id', 'ra', 'dec', 'mag_v', 'simbad_id', 'object_type']
    available_cols = [col for col in cols_to_show if col in enhanced_catalog.columns]
    print(enhanced_catalog[available_cols].to_string(index=False, float_format='%.3f'))
    
    # Summary statistics
    print(f"\n{'='*60}")
    print("SUMMARY STATISTICS")
    print("="*60)
    
    match_stats = enhanced_catalog['object_type'].value_counts()
    print(f"Object type distribution:")
    for otype, count in match_stats.items():
        if otype != '':  # Skip empty types
            print(f"  {otype}: {count}")
    
    unmatched = len(enhanced_catalog[enhanced_catalog['simbad_id'] == ''])
    print(f"\nUnmatched objects: {unmatched} / {len(enhanced_catalog)}")
    
    print(f"\n{'='*60}")
    print("WORKFLOW COMPLETED SUCCESSFULLY")
    print("="*60)
    print("\nKey Features Demonstrated:")
    print("[OK] SIMBAD field querying with rate limiting")
    print("[OK] KD-tree spatial indexing (when available)")
    print("[OK] Object type population from OTYPE/OTYPE_S")
    print("[OK] VizieR photometry enrichment with dtype handling")
    print("[OK] Automatic string conversion for ID columns")
    print("[OK] Error-resilient table stacking")
    
    return enhanced_catalog

if __name__ == '__main__':
    enhanced_catalog = complete_workflow_example()