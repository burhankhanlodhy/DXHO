"""
SIMBAD Data Enhancement Module
Automatically fetch missing SIMBAD data and source names for astronomical objects.
"""

import pandas as pd
import numpy as np
import requests
import time
import json
import logging
from datetime import datetime
from typing import Optional, Dict, List, Tuple, Any
import sqlite3
from contextlib import contextmanager
import os
from astropy.coordinates import SkyCoord
import astropy.units as u

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SIMLADCache:
    """Simple cache for SIMBAD queries to avoid repeated API calls."""
    
    def __init__(self):
        self.cache = {}
    
    def _make_key(self, ra: float, dec: float, radius: float) -> str:
        """Create cache key from rounded coordinates and radius."""
        # Round to avoid floating point precision issues
        ra_round = round(ra, 6)
        dec_round = round(dec, 6)
        radius_round = round(radius, 1)
        return f"{ra_round}_{dec_round}_{radius_round}"
    
    def get(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """Get cached result if available."""
        key = self._make_key(ra, dec, radius)
        return self.cache.get(key)
    
    def set(self, ra: float, dec: float, radius: float, result: Dict):
        """Store result in cache."""
        key = self._make_key(ra, dec, radius)
        self.cache[key] = result


class SIMBADEnhancer:
    """Main class for enhancing astronomical data with SIMBAD information."""
    
    def __init__(self, db_path: Optional[str] = None):
        """
        Initialize SIMBAD enhancer.
        
        Args:
            db_path: Path to SQLite database file
        """
        self.db_path = db_path or 'deho_observatory.db'
        self.cache = SIMLADCache()
        self.misses = []
        self.batch_size = 200
        self.delay_between_batches = 1.0  # seconds
        
        # SIMBAD TAP service URL
        self.simbad_url = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
    
    @contextmanager
    def get_db_connection(self):
        """Get database connection with proper cleanup."""
        conn = sqlite3.connect(self.db_path)
        try:
            yield conn
        finally:
            conn.close()
    
    def find_missing_simbad_rows(self, csv_path: str) -> pd.DataFrame:
        """
        Find rows in CSV that need SIMBAD data.
        
        Args:
            csv_path: Path to CSV file
            
        Returns:
            DataFrame with rows needing SIMBAD data
        """
        logger.info(f"Loading CSV file: {csv_path}")
        df = pd.read_csv(csv_path)
        
        # Find rows where simbad_main_id is blank/NaN and ra/dec are valid
        missing_mask = (
            (df['simbad_main_id'].isna() | (df['simbad_main_id'] == '')) &
            df['ra'].notna() & 
            df['dec'].notna() &
            (df['ra'] >= 0) & (df['ra'] <= 360) &
            (df['dec'] >= -90) & (df['dec'] <= 90)
        )
        
        missing_rows = df[missing_mask].copy()
        logger.info(f"Found {len(missing_rows)} rows needing SIMBAD data out of {len(df)} total rows")
        
        return df, missing_rows
    
    def query_simbad_cone(self, ra: float, dec: float, radius: float) -> Optional[Dict]:
        """
        Query SIMBAD for objects within radius of coordinates.
        
        Args:
            ra: Right ascension in degrees
            dec: Declination in degrees  
            radius: Search radius in arcseconds
            
        Returns:
            Dictionary with SIMBAD data or None if no match
        """
        # Check cache first
        cached_result = self.cache.get(ra, dec, radius)
        if cached_result is not None:
            logger.debug(f"Cache hit for RA={ra}, Dec={dec}, radius={radius}")
            return cached_result
        
        try:
            # Use SIMBAD script interface instead of TAP for better reliability
            script_url = "https://simbad.cds.unistra.fr/simbad/sim-script"
            
            # Build script query
            script_query = f"""
format object "%MAIN_ID|%OTYPE|%SP_TYPE|%RVZ_RADVEL|%DISTANCE_RESULT|%RA(d)|%DEC(d)"
query coo {ra} {dec} radius={radius}s
"""
            
            response = requests.post(
                script_url, 
                data={'script': script_query},
                timeout=30,
                headers={'User-Agent': 'DEHO-Observatory/1.0'}
            )
            response.raise_for_status()
            
            text = response.text.strip()
            
            # Parse response
            lines = text.split('\n')
            data_lines = [line for line in lines if line.strip() and not line.startswith('::')]
            
            if data_lines:
                # Get the first result (nearest object)
                first_line = data_lines[0].strip()
                parts = first_line.split('|')
                
                if len(parts) >= 7:
                    # Parse coordinates and calculate distance
                    try:
                        simbad_ra = float(parts[5]) if parts[5] and parts[5] != '' else np.nan
                        simbad_dec = float(parts[6]) if parts[6] and parts[6] != '' else np.nan
                        
                        # Calculate match distance using astropy
                        if not (np.isnan(simbad_ra) or np.isnan(simbad_dec)):
                            coord1 = SkyCoord(ra*u.deg, dec*u.deg)
                            coord2 = SkyCoord(simbad_ra*u.deg, simbad_dec*u.deg) 
                            match_distance = coord1.separation(coord2).arcsec
                        else:
                            match_distance = np.nan
                            
                    except (ValueError, IndexError):
                        simbad_ra = np.nan
                        simbad_dec = np.nan
                        match_distance = np.nan
                    
                    result = {
                        'main_id': parts[0] if parts[0] and parts[0] != '' else '',
                        'otype': parts[1] if parts[1] and parts[1] != '' else '',
                        'sp_type': parts[2] if parts[2] and parts[2] != '' else '',
                        'rv_value': float(parts[3]) if parts[3] and parts[3] != '' else np.nan,
                        'distance_result': float(parts[4]) if parts[4] and parts[4] != '' else np.nan,
                        'simbad_ra': simbad_ra,
                        'simbad_dec': simbad_dec,
                        'match_distance_arcsec': match_distance
                    }
                    
                    # Cache the result
                    self.cache.set(ra, dec, radius, result)
                    return result
            
            # No match found
            result = None
            self.cache.set(ra, dec, radius, result)
            return result
                
        except Exception as e:
            logger.warning(f"SIMBAD query failed for RA={ra}, Dec={dec}, radius={radius}: {e}")
            return None
    
    def get_source_name(self, simbad_data: Dict) -> str:
        """
        Extract source name from SIMBAD data.
        
        Args:
            simbad_data: Dictionary containing SIMBAD query result
            
        Returns:
            Source name string
        """
        if not simbad_data:
            return ''
        
        main_id = simbad_data.get('main_id', '')
        
        # Clean up the main_id to create a readable source name
        if main_id:
            # Remove extra spaces and common prefixes for cleaner names
            source_name = main_id.strip()
            
            # If it starts with common catalog prefixes, keep as is
            # Otherwise, use the full main_id
            return source_name
        
        return ''
    
    def enhance_single_row(self, row: pd.Series, row_index: int) -> Dict:
        """
        Enhance a single row with SIMBAD data.
        
        Args:
            row: Pandas Series representing the row
            row_index: Original index of the row
            
        Returns:
            Dictionary with enhanced data
        """
        ra = row['ra']
        dec = row['dec']
        
        # Try 3 arcsecond search first
        simbad_data = self.query_simbad_cone(ra, dec, 3.0)
        
        # If no match, try 5 arcsecond search
        if not simbad_data:
            simbad_data = self.query_simbad_cone(ra, dec, 5.0)
            search_radius = 5.0
        else:
            search_radius = 3.0
        
        if simbad_data:
            # Extract source name
            source_name = self.get_source_name(simbad_data)
            
            return {
                'index': row_index,
                'simbad_main_id': simbad_data['main_id'],
                'otype': simbad_data['otype'],
                'sp_type': simbad_data['sp_type'],
                'rv_value': simbad_data['rv_value'],
                'distance_result': simbad_data['distance_result'],
                'source_name': source_name,
                'search_radius': search_radius,
                'match_distance_arcsec': simbad_data['match_distance_arcsec']
            }
        else:
            # Log miss for later analysis
            self.misses.append({
                'index': row_index,
                'ra': ra,
                'dec': dec,
                'radius_tried': search_radius
            })
            
            return {
                'index': row_index,
                'simbad_main_id': '',
                'otype': '',
                'sp_type': '',
                'rv_value': np.nan,
                'distance_result': np.nan,
                'source_name': '',
                'search_radius': search_radius,
                'match_distance_arcsec': np.nan
            }
    
    def process_batch(self, batch_rows: pd.DataFrame) -> List[Dict]:
        """
        Process a batch of rows with SIMBAD queries.
        
        Args:
            batch_rows: DataFrame containing rows to process
            
        Returns:
            List of enhancement dictionaries
        """
        results = []
        
        for idx, (original_idx, row) in enumerate(batch_rows.iterrows()):
            logger.info(f"Processing row {original_idx} ({idx+1}/{len(batch_rows)})")
            result = self.enhance_single_row(row, original_idx)
            results.append(result)
            
            # Small delay between individual queries to be gentle
            if idx < len(batch_rows) - 1:  # Don't delay after last item
                time.sleep(0.1)
        
        return results
    
    def enhance_csv_with_simbad(self, csv_path: str, output_path: str = None) -> str:
        """
        Enhance CSV file with SIMBAD data.
        
        Args:
            csv_path: Input CSV file path
            output_path: Output CSV file path (auto-generated if None)
            
        Returns:
            Path to output CSV file
        """
        if output_path is None:
            base_name = os.path.splitext(csv_path)[0]
            output_path = f"{base_name}_with_simbad.csv"
        
        # Find missing rows
        df, missing_rows = self.find_missing_simbad_rows(csv_path)
        
        if len(missing_rows) == 0:
            logger.info("No rows need SIMBAD enhancement")
            return csv_path
        
        # Add source_name column if it doesn't exist
        if 'source_name' not in df.columns:
            df['source_name'] = ''
        
        # Process in batches
        all_results = []
        
        for batch_start in range(0, len(missing_rows), self.batch_size):
            batch_end = min(batch_start + self.batch_size, len(missing_rows))
            batch = missing_rows.iloc[batch_start:batch_end]
            
            logger.info(f"Processing batch {batch_start//self.batch_size + 1} "
                       f"(rows {batch_start+1} to {batch_end})")
            
            batch_results = self.process_batch(batch)
            all_results.extend(batch_results)
            
            # Delay between batches
            if batch_end < len(missing_rows):
                logger.info(f"Waiting {self.delay_between_batches}s before next batch...")
                time.sleep(self.delay_between_batches)
        
        # Apply results to dataframe
        for result in all_results:
            idx = result['index']
            df.at[idx, 'simbad_main_id'] = result['simbad_main_id']
            df.at[idx, 'otype'] = result['otype'] 
            df.at[idx, 'sp_type'] = result['sp_type']
            df.at[idx, 'rv_value'] = result['rv_value']
            df.at[idx, 'distance_result'] = result['distance_result']
            df.at[idx, 'source_name'] = result['source_name']
        
        # Save enhanced CSV
        df.to_csv(output_path, index=False)
        logger.info(f"Enhanced CSV saved to: {output_path}")
        
        # Save misses log
        if self.misses:
            misses_path = f"{os.path.splitext(output_path)[0]}.misses.csv"
            misses_df = pd.DataFrame(self.misses)
            misses_df.to_csv(misses_path, index=False)
            logger.info(f"Logged {len(self.misses)} misses to: {misses_path}")
        
        return output_path
    
    def update_database_schema(self):
        """Add source_name column to database if it doesn't exist."""
        with self.get_db_connection() as conn:
            cursor = conn.cursor()
            
            # Check if source_name column exists in photometry_sources
            cursor.execute("PRAGMA table_info(photometry_sources)")
            columns = [row[1] for row in cursor.fetchall()]
            
            if 'source_name' not in columns:
                logger.info("Adding source_name column to photometry_sources table")
                cursor.execute(
                    "ALTER TABLE photometry_sources ADD COLUMN source_name TEXT"
                )
                conn.commit()
                
            # Check if source_name column exists in simbad_matches  
            cursor.execute("PRAGMA table_info(simbad_matches)")
            columns = [row[1] for row in cursor.fetchall()]
            
            if 'source_name' not in columns:
                logger.info("Adding source_name column to simbad_matches table")
                cursor.execute(
                    "ALTER TABLE simbad_matches ADD COLUMN source_name TEXT"
                )
                conn.commit()
    
    def sync_csv_to_database(self, csv_path: str):
        """
        Sync enhanced CSV data back to database.
        
        Args:
            csv_path: Path to enhanced CSV file
        """
        logger.info(f"Syncing CSV data to database: {csv_path}")
        
        # Ensure schema is updated
        self.update_database_schema()
        
        df = pd.read_csv(csv_path)
        
        with self.get_db_connection() as conn:
            cursor = conn.cursor()
            
            # Update existing records with new SIMBAD data
            for _, row in df.iterrows():
                if pd.notna(row.get('simbad_main_id')) and row['simbad_main_id']:
                    # Update photometry_sources with source_name
                    if pd.notna(row.get('source_name')):
                        cursor.execute("""
                            UPDATE photometry_sources 
                            SET source_name = ?
                            WHERE ra = ? AND dec = ?
                        """, (row['source_name'], row['ra'], row['dec']))
                    
                    # Insert/update SIMBAD match data
                    cursor.execute("""
                        INSERT OR REPLACE INTO simbad_matches 
                        (source_id, main_id, ra, dec, object_type, spectral_type, 
                         distance_pc, radial_velocity_km_s, source_name, match_timestamp)
                        SELECT 
                            ps.source_id, ?, ?, ?, ?, ?, ?, ?, ?, datetime('now')
                        FROM photometry_sources ps 
                        WHERE ps.ra = ? AND ps.dec = ?
                    """, (
                        row['simbad_main_id'],
                        row['ra'], row['dec'],
                        row.get('otype', ''),
                        row.get('sp_type', ''),
                        row.get('distance_result') if pd.notna(row.get('distance_result')) else None,
                        row.get('rv_value') if pd.notna(row.get('rv_value')) else None,
                        row.get('source_name', ''),
                        row['ra'], row['dec']
                    ))
            
            conn.commit()
            logger.info("Database sync completed")


def main():
    """Main function to run SIMBAD enhancement."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Enhance astronomical data with SIMBAD')
    parser.add_argument('csv_file', help='Path to CSV file to enhance')
    parser.add_argument('--output', '-o', help='Output CSV file path')
    parser.add_argument('--db', help='Database file path')
    parser.add_argument('--batch-size', type=int, default=200, help='Batch size for API calls')
    parser.add_argument('--delay', type=float, default=1.0, help='Delay between batches (seconds)')
    parser.add_argument('--sync-db', action='store_true', help='Sync results to database')
    
    args = parser.parse_args()
    
    # Initialize enhancer
    enhancer = SIMBADEnhancer(db_path=args.db)
    enhancer.batch_size = args.batch_size
    enhancer.delay_between_batches = args.delay
    
    # Enhance CSV
    output_path = enhancer.enhance_csv_with_simbad(args.csv_file, args.output)
    
    # Sync to database if requested
    if args.sync_db:
        enhancer.sync_csv_to_database(output_path)
    
    logger.info("SIMBAD enhancement completed successfully")


if __name__ == '__main__':
    main()