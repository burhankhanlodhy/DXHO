"""
Enhanced Catalog Query Module
Query SIMBAD first, build KD-tree on coordinates, and populate object types.
Optionally enrich with VizieR photometry while avoiding dtype conflicts.
"""

import pandas as pd
import numpy as np
import time
import logging
from datetime import datetime
from typing import Optional, Dict, List, Tuple, Any, Union
from contextlib import contextmanager
import os

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

try:
    from sklearn.neighbors import KDTree
    KDTREE_AVAILABLE = True
except ImportError:
    KDTREE_AVAILABLE = False
    logger.warning("sklearn not available - KD-tree spatial indexing disabled")

try:
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.table import vstack
    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False
    logger.warning("astropy not available")

try:
    from astroquery.simbad import Simbad
    from astroquery.vizier import Vizier
    ASTROQUERY_AVAILABLE = True
except ImportError:
    ASTROQUERY_AVAILABLE = False
    logger.warning("astroquery not available")

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False
    logger.warning("requests not available")


class RateLimiter:
    """Simple rate limiter to prevent API violations."""
    
    def __init__(self, calls_per_second: float = 5.0):
        self.calls_per_second = calls_per_second
        self.min_interval = 1.0 / calls_per_second
        self.last_call = 0.0
    
    def wait(self):
        """Wait if necessary to maintain rate limit."""
        now = time.time()
        elapsed = now - self.last_call
        if elapsed < self.min_interval:
            sleep_time = self.min_interval - elapsed
            time.sleep(sleep_time)
        self.last_call = time.time()


class SIMBADKDTreeCatalog:
    """Enhanced catalog service with SIMBAD first strategy and KD-tree spatial indexing."""
    
    def __init__(self, rate_limit_calls_per_sec: float = 5.0):
        self.rate_limiter = RateLimiter(rate_limit_calls_per_sec)
        self.simbad_cache = {}
        self.vizier_cache = {}
        
        # KD-tree for spatial indexing
        self.simbad_kdtree = None
        self.simbad_catalog_data = None
        self.simbad_coordinates = None
        
        # Configure services
        self._setup_services()
    
    def _setup_services(self):
        """Setup astroquery services with optimal configurations."""
        if ASTROQUERY_AVAILABLE:
            # Configure SIMBAD
            self.simbad = Simbad()
            # Add available fields - check with list_votable_fields() if needed
            try:
                self.simbad.add_votable_fields('otype', 'sp', 'rv')
                # Try to add distance field with correct name
                try:
                    self.simbad.add_votable_fields('parallax')
                except:
                    try:
                        self.simbad.add_votable_fields('distance')
                    except:
                        pass  # Distance field not available
            except Exception as e:
                logger.warning(f"Some SIMBAD fields not available: {e}")
            
            self.simbad.ROW_LIMIT = 50
            
            # Configure VizieR  
            self.vizier = Vizier(row_limit=20)
            self.vizier.columns = ['*']  # Get all columns but handle carefully
            
            logger.info("astroquery services configured with SIMBAD primary strategy")
        else:
            self.simbad = None
            self.vizier = None
            logger.warning("astroquery not available - using HTTP fallbacks only")
    
    def query_simbad_field(self, ra_center: float, dec_center: float, radius_deg: float) -> pd.DataFrame:
        """
        Query SIMBAD for a field and build local catalog with KD-tree indexing.
        
        Args:
            ra_center: Field center RA in degrees
            dec_center: Field center Dec in degrees
            radius_deg: Search radius in degrees
            
        Returns:
            DataFrame with SIMBAD objects in field
        """
        logger.info(f"Querying SIMBAD field: RA={ra_center:.4f}, Dec={dec_center:.4f}, radius={radius_deg:.4f}°")
        
        # Rate limiting
        self.rate_limiter.wait()
        
        simbad_objects = []
        
        if ASTROQUERY_AVAILABLE and self.simbad is not None:
            try:
                # Use astroquery for field query
                coord = SkyCoord(ra_center*u.deg, dec_center*u.deg)
                results = self.simbad.query_region(coord, radius=radius_deg*u.deg)
                
                if results is not None and len(results) > 0:
                    for row in results:
                        try:
                            # Extract coordinates - handle different column names
                            ra = np.nan
                            dec = np.nan
                            
                            # Try different coordinate column names
                            for ra_col in ['RA', 'RA_d', 'RAJ2000', 'RA_ICRS']:
                                if ra_col in row.colnames and row[ra_col] is not np.ma.masked:
                                    ra = float(row[ra_col])
                                    break
                            
                            for dec_col in ['DEC', 'DE_d', 'DEJ2000', 'DEC_ICRS']:
                                if dec_col in row.colnames and row[dec_col] is not np.ma.masked:
                                    dec = float(row[dec_col])
                                    break
                            
                            # Skip if coordinates are invalid
                            if np.isnan(ra) or np.isnan(dec):
                                continue
                            
                            # Extract object data with flexible field names
                            obj_data = {
                                'ra': ra,
                                'dec': dec, 
                                'main_id': str(row['MAIN_ID']) if row['MAIN_ID'] is not np.ma.masked else '',
                                'otype': str(row['OTYPE']) if 'OTYPE' in row.colnames and row['OTYPE'] is not np.ma.masked else '',
                                'otype_s': '',  # OTYPE_S not always available via astroquery
                                'sp_type': str(row['SP_TYPE']) if 'SP_TYPE' in row.colnames and row['SP_TYPE'] is not np.ma.masked else '',
                                'rv_value': float(row['RV_V']) if 'RV_V' in row.colnames and row['RV_V'] is not np.ma.masked else np.nan,
                                'distance_result': np.nan  # Will be calculated from parallax if available
                            }
                            
                            # Handle parallax to distance conversion
                            if 'PLX_VALUE' in row.colnames and row['PLX_VALUE'] is not np.ma.masked:
                                try:
                                    plx_mas = float(row['PLX_VALUE'])
                                    if plx_mas > 0:
                                        obj_data['distance_result'] = 1000.0 / plx_mas  # Convert to pc
                                except:
                                    pass
                            
                            # Use OTYPE_S if available, fallback to OTYPE
                            if obj_data['otype_s']:
                                obj_data['otype_final'] = obj_data['otype_s']
                            else:
                                obj_data['otype_final'] = obj_data['otype']
                            
                            simbad_objects.append(obj_data)
                            
                        except Exception as e:
                            logger.warning(f"Error processing SIMBAD row: {e}")
                            continue
                
                logger.info(f"Retrieved {len(simbad_objects)} objects from SIMBAD astroquery")
                
            except Exception as e:
                logger.warning(f"SIMBAD astroquery failed: {e}")
                # Fall back to HTTP if astroquery fails
                simbad_objects = self._query_simbad_http_field(ra_center, dec_center, radius_deg)
        else:
            # HTTP fallback
            simbad_objects = self._query_simbad_http_field(ra_center, dec_center, radius_deg)
        
        # Convert to DataFrame
        if simbad_objects:
            df = pd.DataFrame(simbad_objects)
            logger.info(f"SIMBAD field query returned {len(df)} objects")
            return df
        else:
            logger.warning("No SIMBAD objects found in field")
            return pd.DataFrame()
    
    def _query_simbad_http_field(self, ra_center: float, dec_center: float, radius_deg: float) -> List[Dict]:
        """HTTP fallback for SIMBAD field queries."""
        if not REQUESTS_AVAILABLE:
            logger.warning("No HTTP capability available")
            return []
        
        try:
            script_url = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
            
            # TAP query for field
            adql_query = f"""
            SELECT MAIN_ID, OTYPE, OTYPE_S, SP_TYPE, RV_VALUE, DISTANCE_RESULT, RA, DEC
            FROM basic JOIN ident ON oidref = oid
            WHERE CONTAINS(POINT('ICRS', RA, DEC), 
                          CIRCLE('ICRS', {ra_center}, {dec_center}, {radius_deg})) = 1
            """
            
            response = requests.post(
                script_url,
                data={
                    'REQUEST': 'doQuery',
                    'LANG': 'ADQL', 
                    'FORMAT': 'TSV',
                    'QUERY': adql_query
                },
                timeout=30,
                headers={'User-Agent': 'DEHO-Observatory/1.0'}
            )
            response.raise_for_status()
            
            lines = response.text.strip().split('\n')
            if len(lines) < 2:  # No data
                return []
            
            headers = lines[0].split('\t')
            objects = []
            
            for line in lines[1:]:
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    try:
                        obj = {}
                        for i, header in enumerate(headers):
                            value = parts[i].strip() if i < len(parts) else ''
                            
                            if header in ['RA', 'DEC', 'RV_VALUE', 'DISTANCE_RESULT']:
                                try:
                                    obj[header.lower()] = float(value) if value and value != '' else np.nan
                                except ValueError:
                                    obj[header.lower()] = np.nan
                            else:
                                obj[header.lower()] = value
                        
                        # Use OTYPE_S if available, fallback to OTYPE
                        otype_s = obj.get('otype_s', '')
                        otype = obj.get('otype', '')
                        obj['otype_final'] = otype_s if otype_s else otype
                        
                        if not (np.isnan(obj.get('ra', np.nan)) or np.isnan(obj.get('dec', np.nan))):
                            objects.append(obj)
                    except Exception as e:
                        logger.debug(f"Error parsing SIMBAD HTTP row: {e}")
                        continue
            
            logger.info(f"Retrieved {len(objects)} objects from SIMBAD HTTP")
            return objects
            
        except Exception as e:
            logger.warning(f"SIMBAD HTTP field query failed: {e}")
            return []
    
    def build_kdtree(self, simbad_df: pd.DataFrame) -> bool:
        """
        Build KD-tree from SIMBAD coordinates for efficient spatial searches.
        
        Args:
            simbad_df: DataFrame with SIMBAD objects containing 'ra' and 'dec' columns
            
        Returns:
            True if KD-tree built successfully
        """
        if not KDTREE_AVAILABLE:
            logger.warning("KD-tree not available - spatial indexing disabled")
            return False
        
        if simbad_df.empty:
            logger.warning("No SIMBAD data to build KD-tree")
            return False
        
        try:
            # Extract coordinates
            coords = simbad_df[['ra', 'dec']].values
            valid_mask = ~(np.isnan(coords).any(axis=1))
            valid_coords = coords[valid_mask]
            valid_data = simbad_df[valid_mask].reset_index(drop=True)
            
            if len(valid_coords) == 0:
                logger.warning("No valid coordinates for KD-tree")
                return False
            
            # Convert to radians for haversine distance
            coords_rad = np.radians(valid_coords)
            
            # Build KD-tree
            self.simbad_kdtree = KDTree(coords_rad, metric='haversine')
            self.simbad_catalog_data = valid_data
            self.simbad_coordinates = valid_coords
            
            logger.info(f"Built KD-tree with {len(valid_coords)} SIMBAD objects")
            return True
            
        except Exception as e:
            logger.error(f"Failed to build KD-tree: {e}")
            return False
    
    def find_nearest_simbad_objects(self, ra: float, dec: float, radius_arcsec: float = 5.0, k: int = 5) -> pd.DataFrame:
        """
        Find nearest SIMBAD objects using KD-tree spatial index.
        
        Args:
            ra: Right ascension in degrees
            dec: Declination in degrees  
            radius_arcsec: Search radius in arcseconds
            k: Maximum number of objects to return
            
        Returns:
            DataFrame with nearest objects and distances
        """
        if self.simbad_kdtree is None or self.simbad_catalog_data is None:
            logger.warning("KD-tree not initialized")
            return pd.DataFrame()
        
        try:
            # Convert search position to radians
            search_coord_rad = np.radians([[ra, dec]])
            
            # Convert radius to radians (for haversine distance on unit sphere)
            radius_rad = np.radians(radius_arcsec / 3600.0)
            
            # Find neighbors within radius
            indices = self.simbad_kdtree.query_radius(search_coord_rad, r=radius_rad, return_distance=False)[0]
            
            if len(indices) == 0:
                return pd.DataFrame()
            
            # Get distances to all found objects
            distances = self.simbad_kdtree.query(search_coord_rad, k=min(k, len(self.simbad_catalog_data)), return_distance=True)
            dist_rad, idx_array = distances
            
            # Filter to only objects within radius
            within_radius = dist_rad[0] <= radius_rad
            valid_indices = idx_array[0][within_radius]
            valid_distances_rad = dist_rad[0][within_radius]
            
            if len(valid_indices) == 0:
                return pd.DataFrame()
            
            # Get matching objects
            matches = self.simbad_catalog_data.iloc[valid_indices].copy()
            
            # Add distance information (convert back to arcseconds)
            matches['match_distance_arcsec'] = np.degrees(valid_distances_rad) * 3600.0
            
            # Sort by distance
            matches = matches.sort_values('match_distance_arcsec').reset_index(drop=True)
            
            return matches
            
        except Exception as e:
            logger.error(f"KD-tree search failed: {e}")
            return pd.DataFrame()
    
    def populate_object_types(self, target_coords: pd.DataFrame, search_radius_arcsec: float = 5.0) -> pd.DataFrame:
        """
        Populate object types for target coordinates using SIMBAD data.
        
        Args:
            target_coords: DataFrame with 'ra' and 'dec' columns
            search_radius_arcsec: Search radius in arcseconds
            
        Returns:
            DataFrame with SIMBAD matches and populated object types
        """
        if target_coords.empty:
            return pd.DataFrame()
        
        results = []
        
        for idx, row in target_coords.iterrows():
            ra, dec = row['ra'], row['dec']
            
            # Find nearest SIMBAD objects using KD-tree
            matches = self.find_nearest_simbad_objects(ra, dec, search_radius_arcsec, k=1)
            
            if not matches.empty:
                best_match = matches.iloc[0]
                
                result = {
                    'target_index': idx,
                    'target_ra': ra,
                    'target_dec': dec,
                    'simbad_main_id': best_match['main_id'],
                    'simbad_otype': best_match['otype_final'],  # Use OTYPE_S preferentially 
                    'simbad_sp_type': best_match.get('sp_type', ''),
                    'simbad_rv_value': best_match.get('rv_value', np.nan),
                    'simbad_distance': best_match.get('distance_result', np.nan),
                    'simbad_ra': best_match['ra'],
                    'simbad_dec': best_match['dec'],
                    'match_distance_arcsec': best_match['match_distance_arcsec']
                }
            else:
                result = {
                    'target_index': idx,
                    'target_ra': ra,
                    'target_dec': dec,
                    'simbad_main_id': '',
                    'simbad_otype': '',
                    'simbad_sp_type': '',
                    'simbad_rv_value': np.nan,
                    'simbad_distance': np.nan,
                    'simbad_ra': np.nan,
                    'simbad_dec': np.nan,
                    'match_distance_arcsec': np.nan
                }
            
            results.append(result)
        
        return pd.DataFrame(results)
    
    def query_vizier_photometry(self, coords_df: pd.DataFrame, radius_arcsec: float = 3.0) -> Optional[pd.DataFrame]:
        """
        Query VizieR for additional photometry, handling dtype conflicts properly.
        
        Args:
            coords_df: DataFrame with coordinates
            radius_arcsec: Search radius
            
        Returns:
            Combined photometric catalog or None
        """
        if not ASTROQUERY_AVAILABLE or self.vizier is None:
            logger.warning("VizieR not available")
            return None
        
        # Catalogs for photometry (avoiding SIMBAD dependency) 
        photometry_catalogs = [
            'II/246/out',      # 2MASS Point Sources
            'I/350/gaiaedr3',  # Gaia EDR3
            'II/328/allwise',  # AllWISE
            'V/147/sdss12',    # SDSS DR12
            'I/239/hip_main'   # Hipparcos
        ]
        
        all_tables = []
        
        for catalog_id in photometry_catalogs:
            try:
                logger.info(f"Querying VizieR catalog: {catalog_id}")
                
                # Rate limiting
                self.rate_limiter.wait()
                
                # Query each coordinate separately for better control
                catalog_results = []
                
                for idx, row in coords_df.iterrows():
                    coord = SkyCoord(row['ra']*u.deg, row['dec']*u.deg)
                    
                    try:
                        tables = self.vizier.query_region(coord, radius=radius_arcsec*u.arcsec, catalog=catalog_id)
                        
                        if tables and len(tables) > 0:
                            table = tables[0]
                            if len(table) > 0:
                                # Convert table to pandas DataFrame
                                df = table.to_pandas()
                                
                                # Convert common problematic ID columns to string
                                id_columns = ['URAT1', 'PS1', 'Source', '_2MASS', 'AllWISE', 'SDSS12', 'HIP']
                                for col in id_columns:
                                    if col in df.columns:
                                        df[col] = df[col].astype(str)
                                
                                # Standardize coordinate column names
                                coord_mapping = {
                                    'RA_ICRS': 'ra_viz',
                                    'DE_ICRS': 'dec_viz', 
                                    'RAJ2000': 'ra_viz',
                                    'DEJ2000': 'dec_viz',
                                    'RA': 'ra_viz',
                                    'DEC': 'dec_viz'
                                }
                                
                                df = df.rename(columns=coord_mapping)
                                
                                # Add source information
                                df['vizier_catalog'] = catalog_id
                                df['target_index'] = idx
                                
                                # Select only useful columns to avoid conflicts
                                useful_cols = ['target_index', 'vizier_catalog']
                                if 'ra_viz' in df.columns:
                                    useful_cols.append('ra_viz')
                                if 'dec_viz' in df.columns:
                                    useful_cols.append('dec_viz')
                                
                                # Add magnitude columns (catalog-specific)
                                mag_cols = [c for c in df.columns if any(x in c.lower() for x in ['mag', 'flux', 'phot'])]
                                useful_cols.extend(mag_cols)
                                
                                # Add ID columns
                                id_cols = [c for c in df.columns if c in id_columns]
                                useful_cols.extend(id_cols)
                                
                                # Select only available columns
                                available_cols = [c for c in useful_cols if c in df.columns]
                                df_subset = df[available_cols]
                                
                                catalog_results.append(df_subset)
                    
                    except Exception as coord_error:
                        logger.debug(f"Coordinate query failed for {catalog_id}: {coord_error}")
                        continue
                
                # Combine results for this catalog
                if catalog_results:
                    try:
                        catalog_combined = pd.concat(catalog_results, ignore_index=True)
                        all_tables.append(catalog_combined)
                        logger.info(f"Collected {len(catalog_combined)} entries from {catalog_id}")
                    except Exception as combine_error:
                        logger.warning(f"Failed to combine {catalog_id} results: {combine_error}")
                        # Add individual tables if combining fails
                        all_tables.extend(catalog_results)
                
            except Exception as catalog_error:
                logger.warning(f"VizieR catalog {catalog_id} failed completely: {catalog_error}")
                continue
        
        # Stack all catalog results with proper error handling
        if all_tables:
            try:
                # Ensure all tables have consistent index before stacking
                for i, table in enumerate(all_tables):
                    all_tables[i] = table.reset_index(drop=True)
                
                # Use pandas concat with careful handling of mixed dtypes
                final_table = pd.concat(all_tables, ignore_index=True, sort=False)
                
                logger.info(f"VizieR photometry: combined {len(final_table)} total entries from {len(all_tables)} catalog queries")
                return final_table
                
            except Exception as stack_error:
                logger.error(f"Failed to stack VizieR tables: {stack_error}")
                # Return largest table if stacking fails
                largest_table = max(all_tables, key=len) if all_tables else None
                if largest_table is not None:
                    logger.info(f"Returning largest table with {len(largest_table)} entries")
                return largest_table
        
        logger.warning("No VizieR photometry data retrieved")
        return None


def main():
    """Example usage of the enhanced catalog system."""
    # Initialize the catalog service
    catalog = SIMBADKDTreeCatalog(rate_limit_calls_per_sec=5.0)
    
    # Example field center (adjust to your field)
    field_ra = 15.5  # degrees
    field_dec = 45.0  # degrees  
    field_radius = 0.5  # degrees
    
    print(f"Querying SIMBAD field: RA={field_ra}, Dec={field_dec}, radius={field_radius}°")
    
    # Step 1: Query SIMBAD for field
    simbad_df = catalog.query_simbad_field(field_ra, field_dec, field_radius)
    
    if simbad_df.empty:
        print("No SIMBAD objects found in field")
        return
    
    print(f"Found {len(simbad_df)} SIMBAD objects")
    
    # Step 2: Build KD-tree
    if catalog.build_kdtree(simbad_df):
        print("KD-tree built successfully")
    else:
        print("Failed to build KD-tree")
        return
    
    # Step 3: Example target coordinates to match
    target_coords = pd.DataFrame({
        'ra': [field_ra + 0.01, field_ra - 0.02, field_ra + 0.05],
        'dec': [field_dec + 0.01, field_dec - 0.02, field_dec - 0.03]
    })
    
    print(f"Matching {len(target_coords)} target coordinates")
    
    # Step 4: Populate object types
    matches = catalog.populate_object_types(target_coords, search_radius_arcsec=10.0)
    
    print(f"Found {len(matches[matches['simbad_main_id'] != ''])} matches with SIMBAD objects")
    print("\nMatches:")
    for _, row in matches.iterrows():
        if row['simbad_main_id']:
            print(f"  {row['simbad_main_id']} ({row['simbad_otype']}) - {row['match_distance_arcsec']:.1f}\"")
    
    # Step 5: Optional VizieR photometry
    print("\nQuerying VizieR for additional photometry...")
    viz_data = catalog.query_vizier_photometry(target_coords, radius_arcsec=5.0)
    
    if viz_data is not None:
        print(f"Retrieved {len(viz_data)} VizieR entries")
        print(f"Catalogs: {viz_data['vizier_catalog'].unique()}")
    else:
        print("No VizieR data retrieved")


if __name__ == '__main__':
    main()