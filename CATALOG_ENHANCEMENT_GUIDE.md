# Comprehensive Catalog Enhancement System

## Overview

This enhanced system addresses your specific requirements:

1. **✅ VizieR Fallback**: Uses VizieR service when SIMBAD TAP is down
2. **✅ Gaia Data Enhancement**: Automatically fills missing Gaia DR3 data
3. **✅ Multiple Fallback Strategies**: Ensures maximum reliability
4. **✅ Source Name Generation**: Creates clean, readable object names

## Key Improvements

### Service Reliability
- **Primary**: astroquery SIMBAD & Gaia services
- **Fallback 1**: VizieR SIMBAD catalog queries  
- **Fallback 2**: Direct HTTP requests to SIMBAD
- **Result**: System works even when individual services are down

### Gaia Data Enhancement
- Automatically detects missing Gaia data (empty `gaia_source_id`)
- Queries Gaia DR3 for: `source_id`, `ra`, `dec`, `phot_g_mean_mag`, `phot_bp_mean_mag`, `phot_rp_mean_mag`, `bp_rp`, `parallax`, `pmra`, `pmdec`
- Uses precise 2" radius matching for better accuracy

### SIMBAD Data Enhancement  
- Detects missing SIMBAD data (empty `simbad_main_id`)
- Queries for: `main_id`, `otype`, `sp_type`, `rv_value`, `distance_result`
- Generates clean `source_name` from SIMBAD identifiers
- Uses 3" then 5" radius search strategy

## Files Created

### 1. `catalog_enhancement.py` - Main Enhancement System
**Advanced multi-service catalog enhancement with comprehensive fallbacks**

**Key Features:**
- **Unified Processing**: Handles both SIMBAD and Gaia data in one pass
- **Service Detection**: Automatically detects available services
- **Smart Caching**: Prevents duplicate API calls
- **Batch Processing**: Configurable batch sizes and delays
- **Miss Logging**: Separate logs for SIMBAD and Gaia misses

**Service Hierarchy:**
```
SIMBAD Data:
├── astroquery.simbad (primary)
├── VizieR SIMBAD catalog (fallback 1) 
└── Direct SIMBAD HTTP (fallback 2)

Gaia Data:
└── astroquery.gaia (Gaia DR3)
```

### 2. `test_catalog_services.py` - Service Testing & Demo
**Test individual services and create sample data**

## Quick Start

### Basic Enhancement
```bash
# Enhance your CSV with both SIMBAD and Gaia data
python catalog_enhancement.py v1.csv

# With custom settings
python catalog_enhancement.py v1.csv --output v1_full_enhanced.csv --batch-size 25 --delay 2.5
```

### Test Services First
```bash
# Test all catalog services
python test_catalog_services.py

# Create sample data for testing
python test_catalog_services.py --create-sample
```

## Expected Results

### Input CSV Analysis
Your `v1.csv` has:
- **58,971 total rows**
- **Missing SIMBAD data**: ~56,957 rows (95.8%)
- **Missing Gaia data**: Most rows (based on empty `gaia_source_id` columns)

### Output Files
After enhancement:
- **`v1_enhanced.csv`**: Complete enhanced dataset
- **`v1_enhanced.simbad_misses.csv`**: Objects without SIMBAD matches  
- **`v1_enhanced.gaia_misses.csv`**: Objects without Gaia matches

### Enhanced Columns

**New/Updated SIMBAD Columns:**
| Column | Description | Example Values |
|--------|-------------|----------------|
| `simbad_main_id` | SIMBAD identifier | `HD 216763`, `TYC 3956-1234-1` |
| `otype` | Object type | `Star`, `EB*`, `Galaxy` |
| `sp_type` | Spectral type | `G5V`, `K2III`, `M3V` |
| `rv_value` | Radial velocity (km/s) | `-12.5`, `8.3` |
| `distance_result` | Distance (pc) | `45.2`, `156.3` |
| `source_name` | Clean readable name | `HD 216763`, `NGC 7380` |

**New/Updated Gaia Columns:**
| Column | Description | Example Values |
|--------|-------------|----------------|
| `gaia_source_id` | Gaia DR3 source ID | `2007132539214702976` |
| `gaia_ra` | Gaia RA (degrees) | `341.80952` |
| `gaia_dec` | Gaia Dec (degrees) | `57.180333` |
| `phot_g_mean_mag` | Gaia G magnitude | `18.206`, `16.631` |
| `phot_bp_mean_mag` | Gaia BP magnitude | `18.767`, `17.355` |
| `phot_rp_mean_mag` | Gaia RP magnitude | `17.333`, `15.812` |
| `bp_rp` | BP-RP color index | `1.434`, `1.544` |
| `parallax` | Parallax (mas) | `0.429`, `0.281` |
| `pmra` | Proper motion RA (mas/yr) | `-1.720`, `-1.383` |
| `pmdec` | Proper motion Dec (mas/yr) | `0.097`, `-3.190` |

## Processing Strategy

### Intelligent Row Detection
The system analyzes your CSV and determines which rows need enhancement:

```python
# SIMBAD needed: empty simbad_main_id AND valid coordinates
simbad_needed = (df['simbad_main_id'].isna() | (df['simbad_main_id'] == '')) & valid_coords

# Gaia needed: empty gaia_source_id AND valid coordinates  
gaia_needed = (df['gaia_source_id'].isna() | (df['gaia_source_id'] == '')) & valid_coords
```

### Batch Processing Logic
```python
# Example processing for your data:
# - 56,957 rows need SIMBAD data
# - ~50,000+ rows need Gaia data
# - Processes both simultaneously to minimize API calls
# - Uses caching to avoid duplicate queries for same coordinates
```

### Service Fallback Flow
```
For each object coordinate:
┌─ Try astroquery SIMBAD
├─ If fails → Try VizieR SIMBAD catalog  
├─ If fails → Try direct SIMBAD HTTP
└─ If all fail → Log to simbad_misses.csv

┌─ Try astroquery Gaia DR3
└─ If fails → Log to gaia_misses.csv
```

## Configuration Options

### Batch Sizes by Dataset Size
```bash
# Small dataset (< 1,000 objects)
--batch-size 10 --delay 1.0

# Medium dataset (1,000 - 10,000 objects) 
--batch-size 25 --delay 2.0

# Large dataset (10,000+ objects) - Your case
--batch-size 50 --delay 3.0

# Very large/slow network
--batch-size 25 --delay 5.0
```

### Search Radii Strategy
- **SIMBAD**: 3" radius first, then 5" if no match
- **Gaia**: 2" radius first (more precise), then 5" if no match
- **Rationale**: Gaia positions are very precise, SIMBAD has more positional uncertainty

## Programmatic Usage

```python
from catalog_enhancement import CatalogEnhancer

# Initialize enhancer
enhancer = CatalogEnhancer()
enhancer.batch_size = 50
enhancer.delay_between_batches = 2.0

# Process your CSV
enhanced_csv = enhancer.enhance_csv_with_catalogs('v1.csv', 'v1_enhanced.csv')

# Check results
print(f"SIMBAD misses: {len(enhancer.simbad_misses)}")
print(f"Gaia misses: {len(enhancer.gaia_misses)}")
```

## Expected Processing Time

For your dataset (~56,957 rows):
- **Batch size 50, delay 2.0s**: ~38 hours
- **Batch size 25, delay 2.0s**: ~76 hours  
- **Batch size 100, delay 1.0s**: ~19 hours (if network/services stable)

**Recommendation**: Start with smaller batch (25) and shorter delay (1.0s), increase if stable.

## Service Status Monitoring

The system provides real-time service status:
```
2025-01-XX XX:XX:XX - INFO - astroquery available - using advanced catalog services
2025-01-XX XX:XX:XX - INFO - Processing batch 1 (rows 1 to 50)
2025-01-XX XX:XX:XX - INFO - Processing row 42 (42/50) [SIMBAD: True, Gaia: True]
2025-01-XX XX:XX:XX - WARNING - SIMBAD astroquery failed for RA=341.809, Dec=57.181: Connection timeout
2025-01-XX XX:XX:XX - DEBUG - Trying VizieR fallback for RA=341.809, Dec=57.181
2025-01-XX XX:XX:XX - INFO - VizieR query successful: HD 216763 (Star)
```

## Error Recovery

### Network Issues
- Automatic service fallbacks
- Configurable timeouts (10s default)
- Retry strategies built into astroquery
- Graceful degradation to HTTP requests

### Service Outages
- **SIMBAD TAP down**: Automatic VizieR fallback
- **VizieR down**: Direct SIMBAD HTTP fallback
- **All SIMBAD down**: Process continues with Gaia-only
- **Gaia down**: Process continues with SIMBAD-only

### Data Validation
- Coordinate range validation
- Magnitude range checking
- Proper motion bounds verification
- NaN handling for missing values

## Success Estimation

Based on your coordinate region (NGC 7380 area):
- **SIMBAD matches**: ~15-30% success rate (typical for photometric surveys)
- **Gaia matches**: ~70-90% success rate (Gaia has excellent sky coverage)
- **Combined enhancement**: Most objects will have at least partial catalog data

## Running the Enhancement

### Recommended Command for Your Dataset
```bash
# Conservative approach for large dataset
python catalog_enhancement.py v1.csv --batch-size 25 --delay 2.0 --output v1_full_enhanced.csv
```

### Monitor Progress
```bash
# Check log output for progress
tail -f catalog_enhancement.log

# Check intermediate results (system saves progress)
wc -l v1_full_enhanced.csv
```

This system provides maximum reliability for catalog enhancement by using multiple services and fallback strategies, ensuring you get the most complete dataset possible even when individual services are unavailable.