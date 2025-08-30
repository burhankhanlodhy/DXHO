"""
Comprehensive Catalog Enhancement Module
Automatically fetch missing SIMBAD and Gaia data using multiple services and fallback strategies.
"""

import pandas as pd
import numpy as np
import time
import logging
from datetime import datetime
from typing import Optional, Dict, List, Tuple, Any
import sqlite3
from contextlib import contextmanager
import os

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False
    logger.warning("requests not available")

try:
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False
    logger.warning("astropy not available for coordinate calculations")

try:
    from astroquery.gaia import Gaia
    from astroquery.simbad import Simbad
    from astroquery.vizier import Vizier
    ASTROQUERY_AVAILABLE = True
    logger.info("astroquery available - using advanced catalog services")
except ImportError:
    ASTROQUERY_AVAILABLE = False
    logger.warning("astroquery not available - using HTTP fallbacks only")


class CatalogCache:
    """Cache for catalog queries to avoid repeated API calls."""
    
    def __init__(self):
        self.simbad_cache = {}
        self.gaia_cache = {}
        self.vizier_cache = {}
    
    def _make_key(self, ra: float, dec: float, radius: float, service: str = "") -> str:
        """Create cache key from coordinates and radius."""
        ra_round = round(ra, 6)
        dec_round = round(dec, 6)
        radius_round = round(radius, 1)
        return f"{service}_{ra_round}_{dec_round}_{radius_round}"
    
    def get_simbad(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """Get cached SIMBAD result."""
        key = self._make_key(ra, dec, radius, "simbad")
        return self.simbad_cache.get(key)
    
    def set_simbad(self, ra: float, dec: float, radius: float, result: Optional[Dict]):
        """Store SIMBAD result in cache."""
        key = self._make_key(ra, dec, radius, "simbad")
        self.simbad_cache[key] = result
    
    def get_gaia(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """Get cached Gaia result."""
        key = self._make_key(ra, dec, radius, "gaia")
        return self.gaia_cache.get(key)
    
    def set_gaia(self, ra: float, dec: float, radius: float, result: Optional[Dict]):
        """Store Gaia result in cache."""
        key = self._make_key(ra, dec, radius, "gaia")
        self.gaia_cache[key] = result


class CatalogEnhancer:
    """Main class for enhancing astronomical data with multiple catalog services."""
    
    def __init__(self, db_path: Optional[str] = None):
        """Initialize catalog enhancer."""
        self.db_path = db_path or 'deho_observatory.db'
        self.cache = CatalogCache()
        self.simbad_misses = []
        self.gaia_misses = []
        self.batch_size = 50
        self.delay_between_batches = 2.0
        
        # Configure services if available
        if ASTROQUERY_AVAILABLE:
            self._setup_astroquery_services()
    
    def _setup_astroquery_services(self):
        """Setup astroquery services with proper configurations."""
        try:
            # Configure Gaia service
            Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
            Gaia.ROW_LIMIT = 50
            
            # Configure VizieR service as PRIMARY - use VizieR B catalogs directly
            self.vizier = Vizier(row_limit=10)
            # Use VizieR's own catalog columns, not SIMBAD
            self.vizier.columns = ['*']  # Get all available columns from VizieR catalogs
            
            # DO NOT use SIMBAD at all - use VizieR B catalogs instead
            self.simbad = None
            
            logger.info("astroquery services configured (VizieR B catalogs PRIMARY, NO SIMBAD)")
            
        except Exception as e:
            logger.warning(f"Failed to configure astroquery services: {e}")
            self.vizier = None
    
    @contextmanager
    def get_db_connection(self):
        """Get database connection with proper cleanup."""
        conn = sqlite3.connect(self.db_path)
        try:
            yield conn
        finally:
            conn.close()
    
    def find_missing_data_rows(self, csv_path: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Find rows missing SIMBAD and/or Gaia data.
        
        Returns:
            (full_df, missing_simbad_rows, missing_gaia_rows)
        """
        logger.info(f"Loading CSV file: {csv_path}")
        df = pd.read_csv(csv_path)
        
        # Find rows missing SIMBAD data
        simbad_missing_mask = (
            (df['simbad_main_id'].isna() | (df['simbad_main_id'] == '')) &
            df['ra'].notna() & df['dec'].notna() &
            (df['ra'] >= 0) & (df['ra'] <= 360) &
            (df['dec'] >= -90) & (df['dec'] <= 90)
        )
        
        # Find rows missing Gaia data
        gaia_missing_mask = (
            (df['gaia_source_id'].isna() | (df['gaia_source_id'] == '')) &
            df['ra'].notna() & df['dec'].notna() &
            (df['ra'] >= 0) & (df['ra'] <= 360) &
            (df['dec'] >= -90) & (df['dec'] <= 90)
        )
        
        missing_simbad = df[simbad_missing_mask].copy()
        missing_gaia = df[gaia_missing_mask].copy()
        
        logger.info(f"Found {len(missing_simbad)} rows needing SIMBAD data")
        logger.info(f"Found {len(missing_gaia)} rows needing Gaia data")
        logger.info(f"Total rows in CSV: {len(df)}")
        
        return df, missing_simbad, missing_gaia
    
    def query_gaia_astroquery(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """Query Gaia DR3 using optimized ADQL with RUWE and magnitude filtering."""
        try:
            coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            
            # Optimized ADQL query instead of simple cone search
            query = f"""
            SELECT TOP 10
                source_id, ra, dec, 
                phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
                bp_rp, parallax, pmra, pmdec, ruwe,
                DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec})) AS match_distance
            FROM gaiadr3.gaia_source 
            WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius/3600.0})) = 1
                AND phot_g_mean_mag < 19.0
                AND phot_g_mean_mag IS NOT NULL
                AND (ruwe IS NULL OR ruwe < 1.4)
                AND (phot_bp_mean_mag IS NOT NULL AND phot_rp_mean_mag IS NOT NULL)
            ORDER BY match_distance ASC
            """
            
            job = Gaia.launch_job_async(query)
            results = job.get_results()
            
            if results is not None and len(results) > 0:
                # Get the nearest source
                row = results[0]
                
                # Calculate match distance
                gaia_coord = SkyCoord(ra=row['ra']*u.deg, dec=row['dec']*u.deg)
                match_distance = coord.separation(gaia_coord).arcsec
                
                return {
                    'gaia_source_id': str(row['source_id']) if 'source_id' in row.colnames else '',
                    'gaia_ra': float(row['ra']) if 'ra' in row.colnames else np.nan,
                    'gaia_dec': float(row['dec']) if 'dec' in row.colnames else np.nan,
                    'phot_g_mean_mag': float(row['phot_g_mean_mag']) if 'phot_g_mean_mag' in row.colnames and row['phot_g_mean_mag'] is not np.ma.masked else np.nan,
                    'phot_bp_mean_mag': float(row['phot_bp_mean_mag']) if 'phot_bp_mean_mag' in row.colnames and row['phot_bp_mean_mag'] is not np.ma.masked else np.nan,
                    'phot_rp_mean_mag': float(row['phot_rp_mean_mag']) if 'phot_rp_mean_mag' in row.colnames and row['phot_rp_mean_mag'] is not np.ma.masked else np.nan,
                    'bp_rp': float(row['bp_rp']) if 'bp_rp' in row.colnames and row['bp_rp'] is not np.ma.masked else np.nan,
                    'parallax': float(row['parallax']) if 'parallax' in row.colnames and row['parallax'] is not np.ma.masked else np.nan,
                    'pmra': float(row['pmra']) if 'pmra' in row.colnames and row['pmra'] is not np.ma.masked else np.nan,
                    'pmdec': float(row['pmdec']) if 'pmdec' in row.colnames and row['pmdec'] is not np.ma.masked else np.nan,
                    'match_distance_arcsec': float(match_distance)
                }
            
            return None
            
        except Exception as e:
            logger.warning(f"Gaia astroquery failed for RA={ra}, Dec={dec}: {e}")
            return None
    
    def query_vizier_primary(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """Query VizieR B catalogs directly as PRIMARY service."""
        try:
            if not ASTROQUERY_AVAILABLE or self.vizier is None:
                return self.query_vizier_http_fallback(ra, dec, radius)
            
            coord = SkyCoord(ra*u.deg, dec*u.deg)
            
            # Try multiple VizieR B catalogs in order of preference
            catalogs_to_try = [
                'I/350/gaiaedr3',  # Gaia EDR3 catalog
                'I/355/gaiadr3',   # Gaia DR3 catalog  
                'B/mk/mktypes',    # MK spectral types
                'B/gcvs/gcvs_cat', # Variable stars
                'I/239/hip_main',  # Hipparcos main catalog
                'I/311/hip2',      # Hipparcos new reduction
                'I/280B/ascc',     # ASCC-2.5 catalog
                'B/pastel/pastel'  # PASTEL stellar parameters
            ]
            
            for catalog in catalogs_to_try:
                try:
                    tables = self.vizier.query_region(coord, radius=radius*u.arcsec, catalog=catalog)
                    
                    if tables and len(tables) > 0:
                        # Get first catalog table
                        table = tables[0]
                        if len(table) > 0:
                            row = table[0]
                            
                            # Calculate match distance
                            cat_ra = None
                            cat_dec = None
                            
                            # Try different coordinate column names
                            for ra_col in ['_RA', 'RA_ICRS', 'RAJ2000', 'RA', 'ra']:
                                if ra_col in row.colnames:
                                    cat_ra = float(row[ra_col]) if not np.ma.is_masked(row[ra_col]) else None
                                    break
                            
                            for dec_col in ['_DE', 'DE_ICRS', 'DEJ2000', 'DEC', 'Dec', 'dec']:
                                if dec_col in row.colnames:
                                    cat_dec = float(row[dec_col]) if not np.ma.is_masked(row[dec_col]) else None
                                    break
                            
                            if cat_ra and cat_dec:
                                cat_coord = SkyCoord(cat_ra*u.deg, cat_dec*u.deg)
                                match_distance = coord.separation(cat_coord).arcsec
                            else:
                                match_distance = np.nan
                            
                            # Extract data based on available columns
                            main_id = ''
                            otype = ''
                            sp_type = ''
                            
                            # Try different identifier columns
                            for id_col in ['Name', 'MAIN_ID', 'HIP', 'HD', 'HR', 'Source']:
                                if id_col in row.colnames and row[id_col]:
                                    main_id = str(row[id_col]).strip()
                                    break
                            
                            # Try different object type columns
                            for type_col in ['VarType', 'SpType', 'OType', 'Class']:
                                if type_col in row.colnames and row[type_col]:
                                    if not otype:  # Don't overwrite if already found
                                        otype = str(row[type_col]).strip()
                            
                            # Try different spectral type columns  
                            for sp_col in ['SpType', 'SpT', 'MK', 'Spectrum']:
                                if sp_col in row.colnames and row[sp_col]:
                                    sp_type = str(row[sp_col]).strip()
                                    break
                            
                            # Get magnitude/distance info
                            rv_value = np.nan
                            distance_result = np.nan
                            
                            # Try radial velocity columns
                            for rv_col in ['RV', 'Vrad', 'RadVel']:
                                if rv_col in row.colnames and not np.ma.is_masked(row[rv_col]):
                                    try:
                                        rv_value = float(row[rv_col])
                                        break
                                    except (ValueError, TypeError):
                                        pass
                            
                            # Try distance/parallax columns
                            for dist_col in ['Dist', 'Distance', 'Plx']:
                                if dist_col in row.colnames and not np.ma.is_masked(row[dist_col]):
                                    try:
                                        if dist_col == 'Plx':  # Convert parallax to distance
                                            plx = float(row[dist_col])
                                            if plx > 0:
                                                distance_result = 1000.0 / plx  # pc
                                        else:
                                            distance_result = float(row[dist_col])
                                        break
                                    except (ValueError, TypeError):
                                        pass
                            
                            if main_id or otype or sp_type:  # Found something useful
                                return {
                                    'main_id': main_id,
                                    'otype': otype,
                                    'sp_type': sp_type,
                                    'rv_value': rv_value,
                                    'distance_result': distance_result,
                                    'simbad_ra': cat_ra or np.nan,
                                    'simbad_dec': cat_dec or np.nan,
                                    'match_distance_arcsec': float(match_distance),
                                    'catalog_source': catalog
                                }
                
                except Exception as cat_error:
                    logger.debug(f"Catalog {catalog} failed: {cat_error}")
                    continue
            
            return None
            
        except Exception as e:
            logger.warning(f"VizieR primary query failed for RA={ra}, Dec={dec}: {e}")
            return None
    
    # VizieR is now PRIMARY - this method removed as it's redundant with query_vizier_primary
    
    def query_vizier_http_fallback(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """HTTP fallback for VizieR B catalog queries."""
        try:
            if not REQUESTS_AVAILABLE:
                return None
            
            # VizieR cone search URL
            vizier_url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
            
            # Try multiple VizieR B catalogs via HTTP
            catalogs_to_try = [
                ('I/350/gaiaedr3', ['Source', 'RA_ICRS', 'DE_ICRS', 'Gmag']),  # Gaia EDR3
                ('I/239/hip_main', ['HIP', 'RAhms', 'DEdms', 'SpType', 'Vmag']),  # Hipparcos
                ('I/280B/ascc', ['recno', 'RAJ2000', 'DEJ2000', 'Bmag', 'Vmag']),  # ASCC
                ('B/mk/mktypes', ['Name', 'RAJ2000', 'DEJ2000', 'SpT']),  # MK types
            ]
            
            for catalog, columns in catalogs_to_try:
                try:
                    params = {
                        '-source': catalog,
                        '-c': f"{ra} {dec}",
                        '-c.rs': radius,  # radius in arcsec
                        '-out.max': '5',
                        '-out': ','.join(columns)
                    }
                    
                    response = requests.get(vizier_url, params=params, timeout=10)
                    response.raise_for_status()
                    
                    lines = response.text.strip().split('\n')
                    
                    # Skip header lines (starting with #)
                    data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
                    
                    if data_lines:
                        # Parse the first data line
                        parts = data_lines[0].split('\t')
                        if len(parts) >= 3:  # At least have name and coordinates
                            try:
                                # Extract basic info
                                main_id = parts[0].strip() if parts[0] and parts[0] != '' else ''
                                
                                # Try to parse coordinates (different formats for different catalogs)
                                cat_ra = np.nan
                                cat_dec = np.nan
                                
                                if catalog == 'I/350/gaiaedr3' or catalog == 'I/355/gaiadr3':
                                    # Gaia uses decimal degrees
                                    if len(parts) >= 3:
                                        cat_ra = float(parts[1]) if parts[1] and parts[1] != '' else np.nan
                                        cat_dec = float(parts[2]) if parts[2] and parts[2] != '' else np.nan
                                
                                elif catalog == 'I/239/hip_main':
                                    # Hipparcos might use HMS/DMS or decimal
                                    if len(parts) >= 3:
                                        try:
                                            cat_ra = float(parts[1]) if parts[1] and parts[1] != '' else np.nan
                                            cat_dec = float(parts[2]) if parts[2] and parts[2] != '' else np.nan
                                        except ValueError:
                                            # Handle HMS/DMS format if needed
                                            pass
                                
                                else:
                                    # Default decimal degree parsing
                                    if len(parts) >= 3:
                                        cat_ra = float(parts[1]) if parts[1] and parts[1] != '' else np.nan
                                        cat_dec = float(parts[2]) if parts[2] and parts[2] != '' else np.nan
                                
                                # Calculate match distance
                                if ASTROPY_AVAILABLE and not (np.isnan(cat_ra) or np.isnan(cat_dec)):
                                    coord1 = SkyCoord(ra*u.deg, dec*u.deg)
                                    coord2 = SkyCoord(cat_ra*u.deg, cat_dec*u.deg)
                                    match_distance = coord1.separation(coord2).arcsec
                                else:
                                    match_distance = np.nan
                                
                                # Extract additional data based on catalog
                                otype = ''
                                sp_type = ''
                                
                                if catalog == 'B/mk/mktypes' and len(parts) >= 4:
                                    sp_type = parts[3].strip() if parts[3] and parts[3] != '' else ''
                                elif catalog == 'I/239/hip_main' and len(parts) >= 4:
                                    sp_type = parts[3].strip() if parts[3] and parts[3] != '' else ''
                                
                                if main_id:  # Found something useful
                                    return {
                                        'main_id': main_id,
                                        'otype': otype,
                                        'sp_type': sp_type,
                                        'rv_value': np.nan,  # Not available in most HTTP queries
                                        'distance_result': np.nan,  # Not available in most HTTP queries
                                        'simbad_ra': cat_ra,
                                        'simbad_dec': cat_dec,
                                        'match_distance_arcsec': match_distance,
                                        'catalog_source': catalog
                                    }
                                
                            except (ValueError, IndexError) as e:
                                logger.debug(f"Error parsing {catalog} response: {e}")
                                continue
                
                except Exception as cat_error:
                    logger.debug(f"HTTP catalog {catalog} failed: {cat_error}")
                    continue
            
            return None
            
        except Exception as e:
            logger.warning(f"VizieR HTTP fallback failed for RA={ra}, Dec={dec}: {e}")
            return None
    
    def query_simbad_http_fallback(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """Direct HTTP fallback for SIMBAD."""
        try:
            if not REQUESTS_AVAILABLE:
                return None
                
            script_url = "https://simbad.cds.unistra.fr/simbad/sim-script"
            
            script_query = f"""
format object "%MAIN_ID|%OTYPE|%SP_TYPE|%RVZ_RADVEL|%DISTANCE_RESULT|%RA(d)|%DEC(d)"
query coo {ra} {dec} radius={radius}s
"""
            
            response = requests.post(
                script_url,
                data={'script': script_query},
                timeout=10,
                headers={'User-Agent': 'DEHO-Observatory/1.0'}
            )
            response.raise_for_status()
            
            text = response.text.strip()
            lines = text.split('\n')
            data_lines = [line for line in lines if line.strip() and not line.startswith('::')]
            
            if data_lines:
                first_line = data_lines[0].strip()
                parts = first_line.split('|')
                
                if len(parts) >= 7:
                    try:
                        simbad_ra = float(parts[5]) if parts[5] and parts[5] != '' else np.nan
                        simbad_dec = float(parts[6]) if parts[6] and parts[6] != '' else np.nan
                        
                        if ASTROPY_AVAILABLE and not (np.isnan(simbad_ra) or np.isnan(simbad_dec)):
                            coord1 = SkyCoord(ra*u.deg, dec*u.deg)
                            coord2 = SkyCoord(simbad_ra*u.deg, simbad_dec*u.deg)
                            match_distance = coord1.separation(coord2).arcsec
                        else:
                            match_distance = np.nan
                        
                        return {
                            'main_id': parts[0] if parts[0] and parts[0] != '' else '',
                            'otype': parts[1] if parts[1] and parts[1] != '' else '',
                            'sp_type': parts[2] if parts[2] and parts[2] != '' else '',
                            'rv_value': float(parts[3]) if parts[3] and parts[3] != '' and parts[3] != '~' else np.nan,
                            'distance_result': float(parts[4]) if parts[4] and parts[4] != '' and parts[4] != '~' else np.nan,
                            'simbad_ra': simbad_ra,
                            'simbad_dec': simbad_dec,
                            'match_distance_arcsec': match_distance
                        }
                    except (ValueError, IndexError):
                        pass
            
            return None
            
        except Exception as e:
            logger.warning(f"SIMBAD HTTP fallback failed for RA={ra}, Dec={dec}: {e}")
            return None
    
    def query_gaia_data(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """Query Gaia data with caching and fallbacks."""
        # Check cache first
        cached_result = self.cache.get_gaia(ra, dec, radius)
        if cached_result is not None:
            logger.debug(f"Gaia cache hit for RA={ra}, Dec={dec}")
            return cached_result
        
        result = None
        
        # Try astroquery first if available
        if ASTROQUERY_AVAILABLE:
            result = self.query_gaia_astroquery(ra, dec, radius)
        
        # Cache result (even if None)
        self.cache.set_gaia(ra, dec, radius, result)
        return result
    
    def query_catalog_data(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """Query catalog data using VizieR B catalogs ONLY (NO SIMBAD)."""
        # Check cache first
        cached_result = self.cache.get_simbad(ra, dec, radius)  # Reuse cache structure
        if cached_result is not None:
            logger.debug(f"VizieR catalog cache hit for RA={ra}, Dec={dec}")
            return cached_result
        
        result = None
        
        # Strategy 1: VizieR B catalogs PRIMARY (astroquery VizieR)
        if ASTROQUERY_AVAILABLE and hasattr(self, 'vizier') and self.vizier is not None:
            result = self.query_vizier_primary(ra, dec, radius)
        
        # Strategy 2: VizieR HTTP fallback (still using VizieR B catalogs)
        if result is None:
            logger.debug(f"Trying VizieR B catalog HTTP fallback for RA={ra}, Dec={dec}")
            result = self.query_vizier_http_fallback(ra, dec, radius)
        
        # NO SIMBAD fallback - only VizieR B catalogs
        
        # Cache result (even if None)
        self.cache.set_simbad(ra, dec, radius, result)  # Reuse cache structure
        return result
    
    def get_source_name(self, simbad_data: Dict) -> str:
        """Extract clean source name from SIMBAD data."""
        if not simbad_data:
            return ''
        
        main_id = simbad_data.get('main_id', '')
        
        if main_id:
            source_name = main_id.strip()
            
            # Clean common prefixes for friendlier names
            if source_name.startswith('NAME '):
                source_name = source_name[5:]
            
            return source_name
        
        return ''
    
    def enhance_single_row(self, row: pd.Series, row_index: int, 
                          need_simbad: bool = True, need_gaia: bool = True) -> Dict:
        """Enhance a single row with catalog data."""
        ra = row['ra']
        dec = row['dec']
        
        result = {
            'index': row_index,
            'simbad_data': None,
            'gaia_data': None
        }
        
        # Query catalog data if needed (VizieR B catalogs, NO SIMBAD)
        if need_simbad:
            # Try 3 arcsecond search first, then 5 arcsecond
            simbad_data = self.query_catalog_data(ra, dec, 3.0)
            if not simbad_data:
                simbad_data = self.query_catalog_data(ra, dec, 5.0)
            
            if simbad_data:
                result['simbad_data'] = {
                    'simbad_main_id': simbad_data['main_id'],
                    'otype': simbad_data['otype'],
                    'sp_type': simbad_data['sp_type'],
                    'rv_value': simbad_data['rv_value'],
                    'distance_result': simbad_data['distance_result'],
                    'source_name': self.get_source_name(simbad_data)
                }
            else:
                self.simbad_misses.append({'index': row_index, 'ra': ra, 'dec': dec})
        
        # Query Gaia data if needed
        if need_gaia:
            # Try 2 arcsecond search for Gaia (more precise)
            gaia_data = self.query_gaia_data(ra, dec, 2.0)
            if not gaia_data:
                gaia_data = self.query_gaia_data(ra, dec, 5.0)
            
            if gaia_data:
                result['gaia_data'] = gaia_data
            else:
                self.gaia_misses.append({'index': row_index, 'ra': ra, 'dec': dec})
        
        return result
    
    def process_batch(self, batch_rows: pd.DataFrame, need_simbad: bool = True, need_gaia: bool = True) -> List[Dict]:
        """Process a batch of rows."""
        results = []
        
        for idx, (original_idx, row) in enumerate(batch_rows.iterrows()):
            logger.info(f"Processing row {original_idx} ({idx+1}/{len(batch_rows)})")
            result = self.enhance_single_row(row, original_idx, need_simbad, need_gaia)
            results.append(result)
            
            # Small delay between queries
            if idx < len(batch_rows) - 1:
                time.sleep(0.2)
        
        return results
    
    def enhance_csv_with_catalogs(self, csv_path: str, output_path: str = None) -> str:
        """Enhance CSV file with catalog data."""
        if output_path is None:
            base_name = os.path.splitext(csv_path)[0]
            output_path = f"{base_name}_enhanced.csv"
        
        # Find missing data
        df, missing_simbad, missing_gaia = self.find_missing_data_rows(csv_path)
        
        # Add source_name column if missing
        if 'source_name' not in df.columns:
            df['source_name'] = ''
        
        # Process all missing data in unified batches
        # Combine indices of rows needing any enhancement
        all_missing_indices = set()
        if len(missing_simbad) > 0:
            all_missing_indices.update(missing_simbad.index)
        if len(missing_gaia) > 0:
            all_missing_indices.update(missing_gaia.index)
        
        if not all_missing_indices:
            logger.info("No rows need catalog enhancement")
            return csv_path
        
        # Convert to sorted list and get rows
        all_missing_indices = sorted(list(all_missing_indices))
        rows_to_process = df.loc[all_missing_indices]
        
        logger.info(f"Processing {len(rows_to_process)} rows needing catalog enhancement")
        
        # Process in batches
        all_results = []
        
        for batch_start in range(0, len(rows_to_process), self.batch_size):
            batch_end = min(batch_start + self.batch_size, len(rows_to_process))
            batch = rows_to_process.iloc[batch_start:batch_end]
            
            logger.info(f"Processing batch {batch_start//self.batch_size + 1} "
                       f"(rows {batch_start+1} to {batch_end})")
            
            # Determine what each row in this batch needs
            simbad_indices = set(missing_simbad.index) if len(missing_simbad) > 0 else set()
            gaia_indices = set(missing_gaia.index) if len(missing_gaia) > 0 else set()
            
            batch_results = []
            for idx, (original_idx, row) in enumerate(batch.iterrows()):
                need_simbad = original_idx in simbad_indices
                need_gaia = original_idx in gaia_indices
                
                logger.info(f"Processing row {original_idx} ({idx+1}/{len(batch)}) "
                           f"[SIMBAD: {need_simbad}, Gaia: {need_gaia}]")
                
                result = self.enhance_single_row(row, original_idx, need_simbad, need_gaia)
                batch_results.append(result)
                
                if idx < len(batch) - 1:
                    time.sleep(0.2)
            
            all_results.extend(batch_results)
            
            # Delay between batches
            if batch_end < len(rows_to_process):
                logger.info(f"Waiting {self.delay_between_batches}s before next batch...")
                time.sleep(self.delay_between_batches)
        
        # Apply results to dataframe
        for result in all_results:
            idx = result['index']
            
            # Apply SIMBAD data
            if result['simbad_data']:
                for col, value in result['simbad_data'].items():
                    df.at[idx, col] = value
            
            # Apply Gaia data
            if result['gaia_data']:
                for col, value in result['gaia_data'].items():
                    if col != 'match_distance_arcsec':  # Don't add this to CSV
                        df.at[idx, col] = value
        
        # Save enhanced CSV
        df.to_csv(output_path, index=False)
        logger.info(f"Enhanced CSV saved to: {output_path}")
        
        # Save miss logs
        if self.simbad_misses:
            simbad_misses_path = f"{os.path.splitext(output_path)[0]}.simbad_misses.csv"
            simbad_misses_df = pd.DataFrame(self.simbad_misses)
            simbad_misses_df.to_csv(simbad_misses_path, index=False)
            logger.info(f"Logged {len(self.simbad_misses)} SIMBAD misses to: {simbad_misses_path}")
        
        if self.gaia_misses:
            gaia_misses_path = f"{os.path.splitext(output_path)[0]}.gaia_misses.csv"
            gaia_misses_df = pd.DataFrame(self.gaia_misses)
            gaia_misses_df.to_csv(gaia_misses_path, index=False)
            logger.info(f"Logged {len(self.gaia_misses)} Gaia misses to: {gaia_misses_path}")
        
        return output_path


def main():
    """Main function to run catalog enhancement."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Enhance astronomical data with SIMBAD and Gaia catalogs')
    parser.add_argument('csv_file', help='Path to CSV file to enhance')
    parser.add_argument('--output', '-o', help='Output CSV file path')
    parser.add_argument('--batch-size', type=int, default=50, help='Batch size for API calls')
    parser.add_argument('--delay', type=float, default=2.0, help='Delay between batches (seconds)')
    
    args = parser.parse_args()
    
    # Initialize enhancer
    enhancer = CatalogEnhancer()
    enhancer.batch_size = args.batch_size
    enhancer.delay_between_batches = args.delay
    
    # Enhance CSV
    output_path = enhancer.enhance_csv_with_catalogs(args.csv_file, args.output)
    
    logger.info("Catalog enhancement completed successfully")


if __name__ == '__main__':
    main()