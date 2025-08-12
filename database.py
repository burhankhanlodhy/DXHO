"""
Database module for the Digital Exoplanet Hunting Observatory.
SQLite integration for astronomical data storage.
"""

import sqlite3
import math
import numpy as np
import pandas as pd
from datetime import datetime
import logging
import os
from contextlib import contextmanager
import json
from typing import Optional, Dict, List, Any, Tuple

logger = logging.getLogger(__name__)


class DatabaseConfig:
    """Database configuration management."""
    
    def __init__(self):
        """Initialize database configuration with defaults."""
        self.config = {
            'database': os.getenv('DEHO_DB_PATH', os.path.join(os.path.dirname(__file__), 'deho_observatory.db'))
        }
    
    def from_file(self, config_file: str) -> None:
        """Load configuration from JSON file."""
        try:
            with open(config_file, 'r') as f:
                file_config = json.load(f)
                self.config.update(file_config)
        except Exception as e:
            logger.warning(f"Could not load database config from {config_file}: {e}")
    
    def get_database_path(self) -> str:
        """Get SQLite database path."""
        return self.config['database']


class ObservatoryDatabase:
    """Main database interface for astronomical data storage."""
    
    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize database connection.
        
        Args:
            config: Database configuration object
        """
        self.config = config or DatabaseConfig()
        self.connection = None
        self.db_path = self.config.get_database_path()
        
        # Try to load config from file if it exists
        config_file = os.path.join(os.path.dirname(__file__), 'database_config.json')
        if os.path.exists(config_file):
            self.config.from_file(config_file)
            self.db_path = self.config.get_database_path()
    
    def initialize_connection(self) -> bool:
        """
        Initialize database connection and create database if it doesn't exist.
        
        Returns:
            bool: True if successful
        """
        try:
            # Create directory if it doesn't exist
            db_dir = os.path.dirname(self.db_path)
            if db_dir and not os.path.exists(db_dir):
                os.makedirs(db_dir)
                
            self.connection = sqlite3.connect(self.db_path, check_same_thread=False)
            self.connection.row_factory = sqlite3.Row  # Enable column access by name
            
            # Enable foreign key constraints
            self.connection.execute("PRAGMA foreign_keys = ON")
            
            logger.info(f"SQLite database connection initialized: {self.db_path}")
            return True
        except Exception as e:
            logger.error(f"Failed to initialize database connection: {e}")
            return False
    
    @contextmanager
    def get_connection(self):
        """
        Context manager for database connections.
        
        Yields:
            sqlite3.Connection: Database connection
        """
        if not self.connection:
            if not self.initialize_connection():
                raise RuntimeError("Failed to initialize database connection")
        
        try:
            yield self.connection
        except Exception as e:
            self.connection.rollback()
            logger.error(f"Database error: {e}")
            raise
        finally:
            self.connection.commit()
    
    def test_connection(self) -> bool:
        """
        Test database connection and basic functionality.
        
        Returns:
            bool: True if connection is working
        """
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                
                # Test basic connection
                cursor.execute("SELECT sqlite_version();")
                version = cursor.fetchone()[0]
                logger.info(f"Connected to SQLite version: {version}")
                
                # List existing tables
                cursor.execute("""
                    SELECT name FROM sqlite_master 
                    WHERE type='table' AND name NOT LIKE 'sqlite_%';
                """)
                
                tables = [row[0] for row in cursor.fetchall()]
                logger.info(f"Available tables: {tables}")
                
                return True
                    
        except Exception as e:
            logger.error(f"Database connection test failed: {e}")
            return False
    
    def create_schema(self, schema_file: Optional[str] = None) -> bool:
        """
        Create database schema from SQL file.
        
        Args:
            schema_file: Path to SQL schema file
            
        Returns:
            bool: True if successful
        """
        if not schema_file:
            schema_file = os.path.join(os.path.dirname(__file__), 'database_schema_sqlite.sql')
        
        try:
            with open(schema_file, 'r') as f:
                schema_sql = f.read()
            
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.executescript(schema_sql)
                logger.info("SQLite database schema created successfully")
                return True
                    
        except Exception as e:
            logger.error(f"Failed to create database schema: {e}")
            return False
    
    def _convert_nan_to_none(self, value: Any) -> Any:
        """Convert numpy NaN/inf to None for SQLite NULL."""
        if isinstance(value, (int, float)) and (np.isnan(value) or np.isinf(value)):
            return None
        return value
    
    def _calculate_angular_distance(self, ra1: float, dec1: float, ra2: float, dec2: float) -> float:
        """
        Calculate angular distance between two points using haversine formula.
        
        Args:
            ra1, dec1: First point coordinates in degrees
            ra2, dec2: Second point coordinates in degrees
            
        Returns:
            Angular distance in degrees
        """
        # Convert to radians
        ra1_rad = math.radians(ra1)
        dec1_rad = math.radians(dec1)
        ra2_rad = math.radians(ra2)
        dec2_rad = math.radians(dec2)
        
        # Haversine formula
        delta_ra = ra2_rad - ra1_rad
        delta_dec = dec2_rad - dec1_rad
        
        a = (math.sin(delta_dec / 2) ** 2 + 
             math.cos(dec1_rad) * math.cos(dec2_rad) * math.sin(delta_ra / 2) ** 2)
        c = 2 * math.asin(math.sqrt(a))
        
        # Convert back to degrees
        return math.degrees(c)
    
    def create_session(self, session_data: Dict[str, Any]) -> Optional[int]:
        """
        Create a new observatory session record.
        
        Args:
            session_data: Dictionary with session information
            
        Returns:
            int: Session ID if successful, None otherwise
        """
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                
                # Convert NaN values to None
                clean_data = {k: self._convert_nan_to_none(v) for k, v in session_data.items()}
                
                # Convert boolean values to integers for SQLite
                clean_data['astrometry_solved'] = int(clean_data.get('astrometry_solved', False))
                clean_data['wcs_available'] = int(clean_data.get('wcs_available', False))
                
                insert_query = """
                    INSERT INTO observatory_sessions (
                        session_name, fits_file_path, image_center_ra, image_center_dec,
                        field_of_view_deg, pixel_scale_arcsec, observation_date,
                        camera_info, exposure_time, iso_value, aperture_fnum,
                        focal_length_mm, fwhm_pixels, detection_threshold,
                        total_sources_detected, astrometry_solved, wcs_available, notes
                    ) VALUES (
                        :session_name, :fits_file_path, :image_center_ra, :image_center_dec,
                        :field_of_view_deg, :pixel_scale_arcsec, :observation_date,
                        :camera_info, :exposure_time, :iso_value, :aperture_fnum,
                        :focal_length_mm, :fwhm_pixels, :detection_threshold,
                        :total_sources_detected, :astrometry_solved, :wcs_available, :notes
                    );
                """
                
                cursor.execute(insert_query, clean_data)
                session_id = cursor.lastrowid
                
                logger.info(f"Created observatory session {session_id}: {clean_data.get('session_name', 'Unknown')}")
                return session_id
                    
        except Exception as e:
            logger.error(f"Failed to create observatory session: {e}")
            return None
    
    def insert_photometry_sources(self, session_id: int, sources_data: List[Dict[str, Any]]) -> List[int]:
        """
        Insert photometry sources for a session.
        
        Args:
            session_id: Observatory session ID
            sources_data: List of source data dictionaries
            
        Returns:
            List[int]: List of inserted source IDs
        """
        source_ids = []
        
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                
                insert_query = """
                    INSERT INTO photometry_sources (
                        session_id, pixel_x, pixel_y, ra, dec, flux, flux_err,
                        background_flux, aperture_radius, fwhm, ellipticity, theta,
                        peak_counts, signal_to_noise, saturated, edge_source
                    ) VALUES (
                        :session_id, :pixel_x, :pixel_y, :ra, :dec,
                        :flux, :flux_err, :background_flux, :aperture_radius,
                        :fwhm, :ellipticity, :theta, :peak_counts,
                        :signal_to_noise, :saturated, :edge_source
                    );
                """
                
                for source_data in sources_data:
                    # Add session_id and clean NaN values
                    source_data['session_id'] = session_id
                    clean_data = {k: self._convert_nan_to_none(v) for k, v in source_data.items()}
                    
                    # Convert boolean values to integers for SQLite
                    clean_data['saturated'] = int(clean_data.get('saturated', False))
                    clean_data['edge_source'] = int(clean_data.get('edge_source', False))
                    
                    cursor.execute(insert_query, clean_data)
                    source_id = cursor.lastrowid
                    source_ids.append(source_id)
                
                logger.info(f"Inserted {len(source_ids)} photometry sources for session {session_id}")
                
        except Exception as e:
            logger.error(f"Failed to insert photometry sources: {e}")
        
        return source_ids
    
    def insert_gaia_matches(self, matches_data: List[Dict[str, Any]]) -> int:
        """
        Insert Gaia catalog matches.
        
        Args:
            matches_data: List of Gaia match data dictionaries
            
        Returns:
            int: Number of matches inserted
        """
        inserted_count = 0
        
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                
                insert_query = """
                    INSERT OR IGNORE INTO gaia_matches (
                        source_id, gaia_source_id, ra, dec, parallax, parallax_error,
                        pmra, pmra_error, pmdec, pmdec_error, phot_g_mean_mag,
                        phot_g_mean_flux, phot_g_mean_flux_error, phot_bp_mean_mag,
                        phot_bp_mean_flux, phot_bp_mean_flux_error, phot_rp_mean_mag,
                        phot_rp_mean_flux, phot_rp_mean_flux_error, bp_rp,
                        phot_g_n_obs, phot_bp_n_obs, phot_rp_n_obs,
                        astrometric_excess_noise, ruwe, match_distance_arcsec
                    ) VALUES (
                        :source_id, :gaia_source_id, :ra, :dec,
                        :parallax, :parallax_error, :pmra, :pmra_error,
                        :pmdec, :pmdec_error, :phot_g_mean_mag,
                        :phot_g_mean_flux, :phot_g_mean_flux_error,
                        :phot_bp_mean_mag, :phot_bp_mean_flux,
                        :phot_bp_mean_flux_error, :phot_rp_mean_mag,
                        :phot_rp_mean_flux, :phot_rp_mean_flux_error,
                        :bp_rp, :phot_g_n_obs, :phot_bp_n_obs,
                        :phot_rp_n_obs, :astrometric_excess_noise,
                        :ruwe, :match_distance_arcsec
                    );
                """
                
                for match_data in matches_data:
                    clean_data = {k: self._convert_nan_to_none(v) for k, v in match_data.items()}
                    cursor.execute(insert_query, clean_data)
                    if cursor.rowcount > 0:
                        inserted_count += 1
                
                logger.info(f"Inserted {inserted_count} Gaia matches")
                
        except Exception as e:
            logger.error(f"Failed to insert Gaia matches: {e}")
            inserted_count = 0
        
        return inserted_count
    
    def insert_simbad_matches(self, matches_data: List[Dict[str, Any]]) -> int:
        """
        Insert SIMBAD catalog matches.
        
        Args:
            matches_data: List of SIMBAD match data dictionaries
            
        Returns:
            int: Number of matches inserted
        """
        inserted_count = 0
        
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                
                insert_query = """
                    INSERT OR IGNORE INTO simbad_matches (
                        source_id, main_id, ra, dec, object_type, spectral_type,
                        distance_pc, radial_velocity_km_s, match_distance_arcsec
                    ) VALUES (
                        :source_id, :main_id, :ra, :dec,
                        :object_type, :spectral_type, :distance_pc,
                        :radial_velocity_km_s, :match_distance_arcsec
                    );
                """
                
                for match_data in matches_data:
                    clean_data = {k: self._convert_nan_to_none(v) for k, v in match_data.items()}
                    cursor.execute(insert_query, clean_data)
                    if cursor.rowcount > 0:
                        inserted_count += 1
                
                logger.info(f"Inserted {inserted_count} SIMBAD matches")
                
        except Exception as e:
            logger.error(f"Failed to insert SIMBAD matches: {e}")
            inserted_count = 0
        
        return inserted_count
    
    def query_sources_in_region(self, ra: float, dec: float, radius_deg: float) -> pd.DataFrame:
        """
        Query sources within a rectangular region (SQLite approximation).
        For circular searches, use _calculate_angular_distance in Python.
        
        Args:
            ra: Right ascension in degrees
            dec: Declination in degrees
            radius_deg: Search radius in degrees
            
        Returns:
            pd.DataFrame: Query results
        """
        try:
            with self.get_connection() as conn:
                # Rectangular approximation for initial filtering
                query = """
                    SELECT * FROM complete_source_catalog
                    WHERE ra IS NOT NULL AND dec IS NOT NULL
                    AND ra BETWEEN ? AND ?
                    AND dec BETWEEN ? AND ?;
                """
                
                ra_min = ra - radius_deg
                ra_max = ra + radius_deg
                dec_min = max(-90, dec - radius_deg)
                dec_max = min(90, dec + radius_deg)
                
                df = pd.read_sql_query(
                    query, conn,
                    params=[ra_min, ra_max, dec_min, dec_max]
                )
                
                # Filter by actual circular distance in Python
                if len(df) > 0:
                    distances = df.apply(
                        lambda row: self._calculate_angular_distance(
                            ra, dec, row['ra'], row['dec']
                        ), axis=1
                    )
                    df = df[distances <= radius_deg].copy()
                    df['distance_deg'] = distances[distances <= radius_deg]
                    df = df.sort_values('distance_deg')
                
                logger.info(f"Found {len(df)} sources within {radius_deg}° of RA={ra}, Dec={dec}")
                return df
                
        except Exception as e:
            logger.error(f"Failed to query sources in region: {e}")
            return pd.DataFrame()
    
    def get_session_statistics(self) -> pd.DataFrame:
        """
        Get statistics for all observatory sessions.
        
        Returns:
            pd.DataFrame: Session statistics
        """
        try:
            with self.get_connection() as conn:
                query = """
                    SELECT * FROM session_statistics 
                    ORDER BY processing_timestamp DESC;
                """
                
                df = pd.read_sql_query(query, conn)
                return df
                
        except Exception as e:
            logger.error(f"Failed to get session statistics: {e}")
            return pd.DataFrame()
    
    def get_maintenance_status(self) -> pd.DataFrame:
        """
        Get database maintenance status.
        
        Returns:
            pd.DataFrame: Table sizes and statistics
        """
        try:
            with self.get_connection() as conn:
                query = """
                    SELECT * FROM maintenance_status 
                    ORDER BY table_name;
                """
                
                df = pd.read_sql_query(query, conn)
                return df
                
        except Exception as e:
            logger.error(f"Failed to get maintenance status: {e}")
            return pd.DataFrame()
    
    def vacuum_database(self) -> str:
        """
        Perform database maintenance (VACUUM).
        
        Returns:
            str: Maintenance log message
        """
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("VACUUM;")
                cursor.execute("ANALYZE;")
                
                logger.info("Database maintenance completed")
                return "Database VACUUM and ANALYZE completed successfully"
                
        except Exception as e:
            logger.error(f"Database maintenance failed: {e}")
            return f"Maintenance failed: {e}"
    
    def close(self):
        """Close database connection."""
        if self.connection:
            self.connection.close()
            self.connection = None
            logger.info("Database connection closed")