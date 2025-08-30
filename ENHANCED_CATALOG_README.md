# Enhanced Catalog Query System

## Overview

This module implements a comprehensive astronomical catalog query system with the following key features:

1. **SIMBAD-First Strategy**: Query SIMBAD first to build a reference catalog for your field
2. **KD-Tree Spatial Indexing**: Efficient spatial searching using sklearn KD-trees
3. **Object Type Population**: Proper handling of SIMBAD OTYPE/OTYPE_S fields
4. **VizieR Enrichment**: Optional photometry and cross-ID enrichment with proper error handling
5. **Rate Limiting**: Built-in protection against API violations
6. **Dtype Conflict Resolution**: Automatic handling of problematic ID columns and table stacking

## Key Features

### SIMBAD Rate Limiting Protection
- Built-in rate limiter (default: 5 calls/second)
- Prevents SIMBAD blocking due to excessive requests
- Configurable rate limits

### KD-Tree Spatial Indexing
- Uses sklearn KD-tree for O(log N) coordinate matching
- Haversine distance metric for accurate sky coordinates
- Falls back gracefully when sklearn unavailable

### Object Type Handling
- Prioritizes SIMBAD OTYPE_S over OTYPE when available
- Extracts spectral types, radial velocities, distances
- Handles missing or invalid data gracefully

### VizieR Integration
- Queries multiple catalogs (2MASS, Gaia, AllWISE, SDSS, Hipparcos)
- Converts problematic ID columns (URAT1, PS1, Source) to strings
- Handles dtype conflicts during table stacking
- Renames coordinate columns consistently
- Error-resilient: one bad catalog doesn't kill the whole query

## Installation

```bash
pip install numpy pandas astropy astroquery scikit-learn requests
```

## Basic Usage

```python
from enhanced_catalog_query import SIMBADKDTreeCatalog
import pandas as pd

# Initialize with rate limiting
catalog = SIMBADKDTreeCatalog(rate_limit_calls_per_sec=5.0)

# Step 1: Query SIMBAD for field
field_ra, field_dec = 83.6, -5.4  # degrees
simbad_catalog = catalog.query_simbad_field(field_ra, field_dec, 0.2)

# Step 2: Build KD-tree spatial index
catalog.build_kdtree(simbad_catalog)

# Step 3: Populate object types for your targets
targets = pd.DataFrame({
    'ra': [83.65, 83.62, 83.68],
    'dec': [-5.39, -5.42, -5.35]
})

matches = catalog.populate_object_types(targets, search_radius_arcsec=10.0)

# Step 4: Optional VizieR enrichment
vizier_data = catalog.query_vizier_photometry(targets, radius_arcsec=5.0)
```

## Complete Workflow Example

See `complete_workflow_example.py` for a full demonstration showing:
- Field querying with fallback strategies
- KD-tree construction and usage
- Object type population and statistics
- VizieR photometry enrichment
- Combined catalog creation

## API Reference

### SIMBADKDTreeCatalog

#### Constructor
```python
catalog = SIMBADKDTreeCatalog(rate_limit_calls_per_sec=5.0)
```

#### Methods

##### query_simbad_field(ra_center, dec_center, radius_deg)
Query SIMBAD for all objects in a field.

**Parameters:**
- `ra_center`: Field center RA in degrees
- `dec_center`: Field center Dec in degrees  
- `radius_deg`: Search radius in degrees

**Returns:** DataFrame with SIMBAD objects

##### build_kdtree(simbad_df)
Build KD-tree spatial index from SIMBAD catalog.

**Parameters:**
- `simbad_df`: DataFrame with SIMBAD objects (must have 'ra', 'dec' columns)

**Returns:** True if successful, False otherwise

##### populate_object_types(target_coords, search_radius_arcsec=5.0)
Populate object types for target coordinates using SIMBAD data.

**Parameters:**
- `target_coords`: DataFrame with 'ra', 'dec' columns
- `search_radius_arcsec`: Search radius in arcseconds

**Returns:** DataFrame with matches and object types

##### query_vizier_photometry(coords_df, radius_arcsec=3.0)
Query VizieR for additional photometry and cross-IDs.

**Parameters:**
- `coords_df`: DataFrame with coordinates  
- `radius_arcsec`: Search radius in arcseconds

**Returns:** Combined photometric catalog DataFrame or None

## Data Handling

### SIMBAD Object Types
The system properly handles SIMBAD object types:
- Uses OTYPE_S (specific type) when available
- Falls back to OTYPE (general type) 
- Provides 'otype_final' column with the best available type

### VizieR Dtype Conflicts
Automatically handles common VizieR issues:
- Converts ID columns (URAT1, PS1, Source) to strings
- Renames coordinate columns consistently (RA_ICRS → ra_viz)
- Uses error-resilient table stacking with metadata_conflicts='silent'
- Individual catalog failures don't break the entire query

### Rate Limiting
Built-in protection against API violations:
- Configurable calls per second limit
- Automatic delays between requests
- Prevents SIMBAD blocking

## Error Handling

The system is designed to be robust:
- Graceful fallbacks when services unavailable
- Individual catalog failures don't break queries  
- Missing data handled with NaN/empty strings
- Comprehensive logging for debugging

## Performance Considerations

### KD-Tree Benefits
- O(log N) coordinate matching vs O(N) linear search
- Significant speedup for large catalogs (>1000 objects)
- Haversine distance for accurate sky coordinates

### Rate Limiting
- Default 5 calls/second prevents blocking
- Adjustable based on your usage patterns
- Batch processing with delays between batches

### Memory Usage
- Caches query results to avoid repeated API calls
- KD-tree memory scales linearly with catalog size
- VizieR results are processed incrementally

## Dependencies

### Required
- numpy >= 1.21.0
- pandas >= 1.3.0  
- astropy >= 5.0.0
- requests >= 2.28.0

### Optional but Recommended
- astroquery >= 0.4.6 (for direct catalog access)
- scikit-learn >= 1.0.0 (for KD-tree spatial indexing)

### Fallback Behavior
- Without astroquery: Uses HTTP API fallbacks
- Without sklearn: Uses slower linear coordinate matching
- Without requests: Limited to astroquery-only functionality

## Testing

Run the test suite:
```bash
python test_enhanced_catalog.py
python complete_workflow_example.py
```

## Troubleshooting

### Common Issues

**No SIMBAD objects found:**
- Try larger search radius
- Check coordinate validity (0 ≤ RA ≤ 360, -90 ≤ Dec ≤ 90)
- Some sky regions have sparse SIMBAD coverage

**VizieR stacking errors:**
- Usually handled automatically with fallback to largest table
- Check logs for specific catalog failures
- Reduce number of catalogs if persistent issues

**Rate limiting:**
- Reduce calls_per_second parameter
- Add delays between large batch operations
- Monitor SIMBAD response times

**KD-tree not available:**
- Install scikit-learn for optimal performance
- System falls back to linear search automatically

## Examples

See the included example files:
- `test_enhanced_catalog.py` - Basic functionality test
- `complete_workflow_example.py` - Full workflow demonstration

## License

This code is part of the DEHO Observatory project and follows the project's licensing terms.