-- Digital Exoplanet Hunting Observatory Database Schema
-- SQLite version for astronomical data storage
-- 
-- Units Convention:
-- - Coordinates: degrees (RA 0-360, Dec -90 to +90)
-- - Parallax: milliarcseconds (mas)
-- - Distances: parsec (pc)
-- - Magnitudes: mag
-- - Store NULL for missing/NaN values

-- Enable foreign key constraints
PRAGMA foreign_keys = ON;

-- Observatory sessions table
CREATE TABLE IF NOT EXISTS observatory_sessions (
    session_id INTEGER PRIMARY KEY AUTOINCREMENT,
    session_name TEXT NOT NULL,
    fits_file_path TEXT NOT NULL,
    image_center_ra REAL CHECK (image_center_ra >= 0 AND image_center_ra <= 360), -- degrees
    image_center_dec REAL CHECK (image_center_dec >= -90 AND image_center_dec <= 90), -- degrees  
    field_of_view_deg REAL CHECK (field_of_view_deg > 0 AND field_of_view_deg <= 180), -- degrees
    pixel_scale_arcsec REAL CHECK (pixel_scale_arcsec > 0), -- arcsec/pixel
    observation_date TEXT, -- ISO 8601 format
    camera_info TEXT,
    exposure_time REAL CHECK (exposure_time > 0), -- seconds
    iso_value INTEGER CHECK (iso_value > 0),
    aperture_fnum REAL CHECK (aperture_fnum > 0),
    focal_length_mm REAL CHECK (focal_length_mm > 0),
    fwhm_pixels REAL CHECK (fwhm_pixels > 0),
    detection_threshold REAL CHECK (detection_threshold > 0),
    total_sources_detected INTEGER CHECK (total_sources_detected >= 0),
    astrometry_solved INTEGER DEFAULT 0, -- 0 = False, 1 = True
    wcs_available INTEGER DEFAULT 0, -- 0 = False, 1 = True
    processing_timestamp TEXT DEFAULT (datetime('now')), -- ISO 8601 format
    notes TEXT,
    
    -- Ensure unique session names per FITS file
    UNIQUE (fits_file_path, session_name)
);

-- Index for coordinate searches
CREATE INDEX IF NOT EXISTS idx_sessions_ra ON observatory_sessions (image_center_ra);
CREATE INDEX IF NOT EXISTS idx_sessions_dec ON observatory_sessions (image_center_dec);
CREATE INDEX IF NOT EXISTS idx_sessions_coords ON observatory_sessions (image_center_ra, image_center_dec);

-- Photometry results table (detected sources)
CREATE TABLE IF NOT EXISTS photometry_sources (
    source_id INTEGER PRIMARY KEY AUTOINCREMENT,
    session_id INTEGER REFERENCES observatory_sessions(session_id) ON DELETE CASCADE,
    
    -- Source position with range validation
    pixel_x REAL NOT NULL CHECK (pixel_x >= 0), -- pixel coordinates
    pixel_y REAL NOT NULL CHECK (pixel_y >= 0),
    ra REAL CHECK (ra >= 0 AND ra <= 360), -- degrees, NULL if no WCS
    dec REAL CHECK (dec >= -90 AND dec <= 90), -- degrees, NULL if no WCS
    
    -- Photometry measurements with sanity checks
    flux REAL NOT NULL, -- counts (can be negative for background subtraction)
    flux_err REAL CHECK (flux_err >= 0), -- counts error (must be positive)
    background_flux REAL, -- counts
    aperture_radius REAL CHECK (aperture_radius > 0), -- pixels
    
    -- Source properties with validation
    fwhm REAL CHECK (fwhm > 0), -- pixels
    ellipticity REAL CHECK (ellipticity >= 0 AND ellipticity <= 1),
    theta REAL CHECK (theta >= -180 AND theta <= 180), -- position angle, degrees
    peak_counts REAL CHECK (peak_counts > 0),
    signal_to_noise REAL,
    
    -- Flags
    saturated INTEGER DEFAULT 0, -- 0 = False, 1 = True
    edge_source INTEGER DEFAULT 0, -- 0 = False, 1 = True
    
    measurement_timestamp TEXT DEFAULT (datetime('now')),
    
    -- Ensure unique pixel coordinates per session (prevent duplicate detections)
    UNIQUE (session_id, pixel_x, pixel_y)
);

-- Spatial indexes for coordinate searches
CREATE INDEX IF NOT EXISTS idx_photometry_ra ON photometry_sources (ra) WHERE ra IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_photometry_dec ON photometry_sources (dec) WHERE dec IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_photometry_coords ON photometry_sources (ra, dec) WHERE ra IS NOT NULL AND dec IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_photometry_session ON photometry_sources(session_id);

-- Gaia DR3 catalog matches
CREATE TABLE IF NOT EXISTS gaia_matches (
    gaia_id INTEGER PRIMARY KEY AUTOINCREMENT,
    source_id INTEGER REFERENCES photometry_sources(source_id) ON DELETE CASCADE,
    
    -- Gaia identifiers
    gaia_source_id INTEGER NOT NULL, -- Gaia DR3 source_id
    
    -- Astrometry (Gaia coordinates are authoritative) with range checks
    ra REAL NOT NULL CHECK (ra >= 0 AND ra <= 360), -- degrees
    dec REAL NOT NULL CHECK (dec >= -90 AND dec <= 90), -- degrees  
    parallax REAL, -- milliarcseconds (mas), can be negative
    parallax_error REAL CHECK (parallax_error >= 0), -- mas
    
    -- Proper motion with reasonable limits
    pmra REAL CHECK (pmra >= -10000 AND pmra <= 10000), -- mas/year
    pmra_error REAL CHECK (pmra_error >= 0), -- mas/year
    pmdec REAL CHECK (pmdec >= -10000 AND pmdec <= 10000), -- mas/year  
    pmdec_error REAL CHECK (pmdec_error >= 0), -- mas/year
    
    -- Photometry with astronomical magnitude limits
    phot_g_mean_mag REAL CHECK (phot_g_mean_mag >= -5 AND phot_g_mean_mag <= 30), -- Gaia G magnitude
    phot_g_mean_flux REAL CHECK (phot_g_mean_flux >= 0), -- e-/s
    phot_g_mean_flux_error REAL CHECK (phot_g_mean_flux_error >= 0), -- e-/s
    phot_bp_mean_mag REAL CHECK (phot_bp_mean_mag >= -5 AND phot_bp_mean_mag <= 30), -- Gaia BP magnitude
    phot_bp_mean_flux REAL CHECK (phot_bp_mean_flux >= 0), -- e-/s
    phot_bp_mean_flux_error REAL CHECK (phot_bp_mean_flux_error >= 0), -- e-/s
    phot_rp_mean_mag REAL CHECK (phot_rp_mean_mag >= -5 AND phot_rp_mean_mag <= 30), -- Gaia RP magnitude
    phot_rp_mean_flux REAL CHECK (phot_rp_mean_flux >= 0), -- e-/s
    phot_rp_mean_flux_error REAL CHECK (phot_rp_mean_flux_error >= 0), -- e-/s
    bp_rp REAL CHECK (bp_rp >= -5 AND bp_rp <= 10), -- BP-RP color
    
    -- Quality indicators with validation
    phot_g_n_obs INTEGER CHECK (phot_g_n_obs >= 0),
    phot_bp_n_obs INTEGER CHECK (phot_bp_n_obs >= 0), 
    phot_rp_n_obs INTEGER CHECK (phot_rp_n_obs >= 0),
    astrometric_excess_noise REAL CHECK (astrometric_excess_noise >= 0),
    ruwe REAL CHECK (ruwe >= 0), -- Renormalized unit weight error
    
    -- Cross-match information
    match_distance_arcsec REAL CHECK (match_distance_arcsec >= 0 AND match_distance_arcsec <= 3600), -- arcseconds
    match_timestamp TEXT DEFAULT (datetime('now')),
    
    -- Prevent duplicate Gaia matches for the same source
    UNIQUE (source_id),
    -- Prevent duplicate source matches to same Gaia source
    UNIQUE (gaia_source_id)
);

-- Spatial indexes for Gaia sources
CREATE INDEX IF NOT EXISTS idx_gaia_ra ON gaia_matches (ra);
CREATE INDEX IF NOT EXISTS idx_gaia_dec ON gaia_matches (dec);
CREATE INDEX IF NOT EXISTS idx_gaia_coords ON gaia_matches (ra, dec);
CREATE INDEX IF NOT EXISTS idx_gaia_source_id ON gaia_matches(gaia_source_id);
CREATE INDEX IF NOT EXISTS idx_gaia_match ON gaia_matches(source_id);

-- SIMBAD catalog matches
CREATE TABLE IF NOT EXISTS simbad_matches (
    simbad_id INTEGER PRIMARY KEY AUTOINCREMENT,
    source_id INTEGER REFERENCES photometry_sources(source_id) ON DELETE CASCADE,
    
    -- SIMBAD identifiers  
    main_id TEXT, -- Main SIMBAD identifier
    
    -- Coordinates (from SIMBAD) with range validation
    ra REAL CHECK (ra >= 0 AND ra <= 360), -- degrees
    dec REAL CHECK (dec >= -90 AND dec <= 90), -- degrees
    
    -- Object classification
    object_type TEXT, -- SIMBAD object type
    spectral_type TEXT, -- Spectral classification
    
    -- Physical properties with validation
    distance_pc REAL CHECK (distance_pc > 0 AND distance_pc < 1e12), -- distance in parsecs
    radial_velocity_km_s REAL CHECK (radial_velocity_km_s >= -100000 AND radial_velocity_km_s <= 100000), -- km/s
    
    -- Cross-match information
    match_distance_arcsec REAL CHECK (match_distance_arcsec >= 0 AND match_distance_arcsec <= 3600), -- arcseconds
    match_timestamp TEXT DEFAULT (datetime('now')),
    
    -- Prevent duplicate SIMBAD matches for the same source
    UNIQUE (source_id),
    -- Prevent duplicate matches to the same SIMBAD object (if main_id provided)
    UNIQUE (main_id)
);

-- Spatial indexes for SIMBAD sources
CREATE INDEX IF NOT EXISTS idx_simbad_ra ON simbad_matches (ra) WHERE ra IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_simbad_dec ON simbad_matches (dec) WHERE dec IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_simbad_coords ON simbad_matches (ra, dec) WHERE ra IS NOT NULL AND dec IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_simbad_match ON simbad_matches(source_id);
CREATE INDEX IF NOT EXISTS idx_simbad_main_id ON simbad_matches(main_id);

-- ============================================================================
-- STAGING TABLES FOR DATA VALIDATION AND QUALITY CONTROL
-- ============================================================================

-- Staging table for photometry sources (relaxed constraints for initial validation)
CREATE TABLE IF NOT EXISTS staging_photometry_sources (
    staging_id INTEGER PRIMARY KEY AUTOINCREMENT,
    session_id INTEGER, -- Not enforced at staging level
    
    -- Source position (no constraints for initial load)
    pixel_x REAL,
    pixel_y REAL,
    ra REAL,
    dec REAL,
    
    -- Photometry measurements
    flux REAL,
    flux_err REAL,
    background_flux REAL,
    aperture_radius REAL,
    
    -- Source properties
    fwhm REAL,
    ellipticity REAL,
    theta REAL,
    peak_counts REAL,
    signal_to_noise REAL,
    
    -- Flags
    saturated INTEGER DEFAULT 0,
    edge_source INTEGER DEFAULT 0,
    
    -- Validation metadata
    validation_status TEXT DEFAULT 'pending', -- pending, valid, invalid, processed
    validation_errors TEXT, -- JSON array of validation error messages
    staging_timestamp TEXT DEFAULT (datetime('now'))
);

-- Create indexes on staging tables for performance
CREATE INDEX IF NOT EXISTS idx_staging_photometry_status ON staging_photometry_sources(validation_status);
CREATE INDEX IF NOT EXISTS idx_staging_photometry_session ON staging_photometry_sources(session_id);

-- Create useful views for common queries

-- Complete source catalog view with all available data
CREATE VIEW IF NOT EXISTS complete_source_catalog AS
SELECT 
    ps.source_id,
    ps.session_id,
    s.session_name,
    s.fits_file_path,
    
    -- Coordinates (prefer Gaia if available, otherwise photometry)
    COALESCE(gm.ra, ps.ra) as ra,
    COALESCE(gm.dec, ps.dec) as dec,
    ps.pixel_x,
    ps.pixel_y,
    
    -- Photometry 
    ps.flux,
    ps.flux_err,
    ps.background_flux,
    ps.signal_to_noise,
    
    -- Gaia data
    gm.gaia_source_id,
    gm.parallax,
    gm.phot_g_mean_mag,
    gm.phot_bp_mean_mag,
    gm.phot_rp_mean_mag,
    gm.bp_rp,
    gm.pmra,
    gm.pmdec,
    
    -- SIMBAD data
    sm.main_id as simbad_name,
    sm.object_type,
    sm.spectral_type,
    sm.distance_pc,
    sm.radial_velocity_km_s,
    
    -- Match quality
    gm.match_distance_arcsec as gaia_match_arcsec,
    sm.match_distance_arcsec as simbad_match_arcsec
    
FROM photometry_sources ps
LEFT JOIN observatory_sessions s ON ps.session_id = s.session_id
LEFT JOIN gaia_matches gm ON ps.source_id = gm.source_id  
LEFT JOIN simbad_matches sm ON ps.source_id = sm.source_id;

-- Statistics view
CREATE VIEW IF NOT EXISTS session_statistics AS
SELECT 
    s.session_id,
    s.session_name,
    s.total_sources_detected,
    COUNT(ps.source_id) as sources_stored,
    COUNT(gm.gaia_id) as gaia_matches,
    COUNT(sm.simbad_id) as simbad_matches,
    s.processing_timestamp,
    s.astrometry_solved,
    s.wcs_available
FROM observatory_sessions s
LEFT JOIN photometry_sources ps ON s.session_id = ps.session_id
LEFT JOIN gaia_matches gm ON ps.source_id = gm.source_id
LEFT JOIN simbad_matches sm ON ps.source_id = sm.source_id  
GROUP BY s.session_id, s.session_name, s.total_sources_detected, s.processing_timestamp, s.astrometry_solved, s.wcs_available
ORDER BY s.processing_timestamp DESC;

-- ============================================================================
-- SQLITE-SPECIFIC FUNCTIONS FOR ASTRONOMICAL CALCULATIONS
-- ============================================================================

-- Note: SQLite doesn't have built-in astronomical functions like PostgreSQL's Q3C
-- We'll implement basic distance calculations using SQLite's math functions

-- Haversine formula for angular distance (not as efficient as Q3C but functional)
-- Distance calculation will be done in Python code for better performance

-- Create maintenance status view
CREATE VIEW IF NOT EXISTS maintenance_status AS
SELECT 
    'observatory_sessions' as table_name,
    COUNT(*) as row_count,
    MAX(processing_timestamp) as last_activity
FROM observatory_sessions
UNION ALL
SELECT 
    'photometry_sources',
    COUNT(*),
    MAX(measurement_timestamp)
FROM photometry_sources
UNION ALL
SELECT 
    'staging_photometry_sources',
    COUNT(*),
    MAX(staging_timestamp)
FROM staging_photometry_sources;

-- Example queries for reference (without Q3C):
-- 
-- -- Find sources within approximate region (rectangular search)
-- SELECT * FROM complete_source_catalog 
-- WHERE ra BETWEEN 179.0 AND 181.0 
--   AND dec BETWEEN -1.0 AND 1.0;
--
-- -- For circular searches, distance calculations should be done in Python
-- -- using the haversine formula or astropy for better precision