-- Digital Exoplanet Hunting Observatory Database Schema
-- PostgreSQL with Q3C extension for astronomical coordinates
-- 
-- Units Convention:
-- - Coordinates: degrees (RA 0-360, Dec -90 to +90)
-- - Parallax: milliarcseconds (mas)
-- - Distances: parsec (pc)
-- - Magnitudes: mag
-- - Store NULL for missing/NaN values

-- Enable Q3C extension for astronomical coordinate queries
CREATE EXTENSION IF NOT EXISTS q3c;

-- Create database schema
CREATE SCHEMA IF NOT EXISTS deho;
SET search_path TO deho, public;

-- Observatory sessions table
CREATE TABLE IF NOT EXISTS observatory_sessions (
    session_id SERIAL PRIMARY KEY,
    session_name VARCHAR(255) NOT NULL,
    fits_file_path TEXT NOT NULL,
    image_center_ra DOUBLE PRECISION CHECK (image_center_ra >= 0 AND image_center_ra <= 360), -- degrees
    image_center_dec DOUBLE PRECISION CHECK (image_center_dec >= -90 AND image_center_dec <= 90), -- degrees  
    field_of_view_deg DOUBLE PRECISION CHECK (field_of_view_deg > 0 AND field_of_view_deg <= 180), -- degrees
    pixel_scale_arcsec DOUBLE PRECISION CHECK (pixel_scale_arcsec > 0), -- arcsec/pixel
    observation_date TIMESTAMP,
    camera_info TEXT,
    exposure_time DOUBLE PRECISION CHECK (exposure_time > 0), -- seconds
    iso_value INTEGER CHECK (iso_value > 0),
    aperture_fnum DOUBLE PRECISION CHECK (aperture_fnum > 0),
    focal_length_mm DOUBLE PRECISION CHECK (focal_length_mm > 0),
    fwhm_pixels DOUBLE PRECISION CHECK (fwhm_pixels > 0),
    detection_threshold DOUBLE PRECISION CHECK (detection_threshold > 0),
    total_sources_detected INTEGER CHECK (total_sources_detected >= 0),
    astrometry_solved BOOLEAN DEFAULT FALSE,
    wcs_available BOOLEAN DEFAULT FALSE,
    processing_timestamp TIMESTAMP DEFAULT NOW(),
    notes TEXT,
    
    -- Ensure unique session names per FITS file
    CONSTRAINT unique_session_per_fits UNIQUE (fits_file_path, session_name)
);

-- Index for coordinate searches
CREATE INDEX IF NOT EXISTS idx_sessions_coords ON observatory_sessions USING btree (q3c_ang2ipix(image_center_ra, image_center_dec));

-- Photometry results table (detected sources)
CREATE TABLE IF NOT EXISTS photometry_sources (
    source_id SERIAL PRIMARY KEY,
    session_id INTEGER REFERENCES observatory_sessions(session_id) ON DELETE CASCADE,
    
    -- Source position with range validation
    pixel_x DOUBLE PRECISION NOT NULL CHECK (pixel_x >= 0), -- pixel coordinates
    pixel_y DOUBLE PRECISION NOT NULL CHECK (pixel_y >= 0),
    ra DOUBLE PRECISION CHECK (ra >= 0 AND ra <= 360), -- degrees, NULL if no WCS
    dec DOUBLE PRECISION CHECK (dec >= -90 AND dec <= 90), -- degrees, NULL if no WCS
    
    -- Photometry measurements with sanity checks
    flux DOUBLE PRECISION NOT NULL, -- counts (can be negative for background subtraction)
    flux_err DOUBLE PRECISION CHECK (flux_err >= 0), -- counts error (must be positive)
    background_flux DOUBLE PRECISION, -- counts
    aperture_radius DOUBLE PRECISION CHECK (aperture_radius > 0), -- pixels
    
    -- Source properties with validation
    fwhm DOUBLE PRECISION CHECK (fwhm > 0), -- pixels
    ellipticity DOUBLE PRECISION CHECK (ellipticity >= 0 AND ellipticity <= 1),
    theta DOUBLE PRECISION CHECK (theta >= -180 AND theta <= 180), -- position angle, degrees
    peak_counts DOUBLE PRECISION CHECK (peak_counts > 0),
    signal_to_noise DOUBLE PRECISION,
    
    -- Flags
    saturated BOOLEAN DEFAULT FALSE,
    edge_source BOOLEAN DEFAULT FALSE,
    
    measurement_timestamp TIMESTAMP DEFAULT NOW(),
    
    -- Ensure unique pixel coordinates per session (prevent duplicate detections)
    CONSTRAINT unique_pixel_per_session UNIQUE (session_id, pixel_x, pixel_y)
);

-- Spatial index for coordinate searches (only when coordinates available)
CREATE INDEX IF NOT EXISTS idx_photometry_coords ON photometry_sources USING btree (q3c_ang2ipix(ra, dec)) WHERE ra IS NOT NULL AND dec IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_photometry_session ON photometry_sources(session_id);

-- Gaia DR3 catalog matches
CREATE TABLE IF NOT EXISTS gaia_matches (
    gaia_id SERIAL PRIMARY KEY,
    source_id INTEGER REFERENCES photometry_sources(source_id) ON DELETE CASCADE,
    
    -- Gaia identifiers
    gaia_source_id BIGINT NOT NULL, -- Gaia DR3 source_id
    
    -- Astrometry (Gaia coordinates are authoritative) with range checks
    ra DOUBLE PRECISION NOT NULL CHECK (ra >= 0 AND ra <= 360), -- degrees
    dec DOUBLE PRECISION NOT NULL CHECK (dec >= -90 AND dec <= 90), -- degrees  
    parallax DOUBLE PRECISION, -- milliarcseconds (mas), can be negative
    parallax_error DOUBLE PRECISION CHECK (parallax_error >= 0), -- mas
    
    -- Proper motion with reasonable limits
    pmra DOUBLE PRECISION CHECK (pmra >= -10000 AND pmra <= 10000), -- mas/year
    pmra_error DOUBLE PRECISION CHECK (pmra_error >= 0), -- mas/year
    pmdec DOUBLE PRECISION CHECK (pmdec >= -10000 AND pmdec <= 10000), -- mas/year  
    pmdec_error DOUBLE PRECISION CHECK (pmdec_error >= 0), -- mas/year
    
    -- Photometry with astronomical magnitude limits
    phot_g_mean_mag DOUBLE PRECISION CHECK (phot_g_mean_mag >= -5 AND phot_g_mean_mag <= 30), -- Gaia G magnitude
    phot_g_mean_flux DOUBLE PRECISION CHECK (phot_g_mean_flux >= 0), -- e-/s
    phot_g_mean_flux_error DOUBLE PRECISION CHECK (phot_g_mean_flux_error >= 0), -- e-/s
    phot_bp_mean_mag DOUBLE PRECISION CHECK (phot_bp_mean_mag >= -5 AND phot_bp_mean_mag <= 30), -- Gaia BP magnitude
    phot_bp_mean_flux DOUBLE PRECISION CHECK (phot_bp_mean_flux >= 0), -- e-/s
    phot_bp_mean_flux_error DOUBLE PRECISION CHECK (phot_bp_mean_flux_error >= 0), -- e-/s
    phot_rp_mean_mag DOUBLE PRECISION CHECK (phot_rp_mean_mag >= -5 AND phot_rp_mean_mag <= 30), -- Gaia RP magnitude
    phot_rp_mean_flux DOUBLE PRECISION CHECK (phot_rp_mean_flux >= 0), -- e-/s
    phot_rp_mean_flux_error DOUBLE PRECISION CHECK (phot_rp_mean_flux_error >= 0), -- e-/s
    bp_rp DOUBLE PRECISION CHECK (bp_rp >= -5 AND bp_rp <= 10), -- BP-RP color
    
    -- Quality indicators with validation
    phot_g_n_obs INTEGER CHECK (phot_g_n_obs >= 0),
    phot_bp_n_obs INTEGER CHECK (phot_bp_n_obs >= 0), 
    phot_rp_n_obs INTEGER CHECK (phot_rp_n_obs >= 0),
    astrometric_excess_noise DOUBLE PRECISION CHECK (astrometric_excess_noise >= 0),
    ruwe DOUBLE PRECISION CHECK (ruwe >= 0), -- Renormalized unit weight error
    
    -- Cross-match information
    match_distance_arcsec DOUBLE PRECISION CHECK (match_distance_arcsec >= 0 AND match_distance_arcsec <= 3600), -- arcseconds
    match_timestamp TIMESTAMP DEFAULT NOW(),
    
    -- Prevent duplicate Gaia matches for the same source
    CONSTRAINT unique_gaia_per_source UNIQUE (source_id),
    -- Prevent duplicate source matches to same Gaia source
    CONSTRAINT unique_source_per_gaia UNIQUE (gaia_source_id)
);

-- Spatial index for Gaia sources
CREATE INDEX IF NOT EXISTS idx_gaia_coords ON gaia_matches USING btree (q3c_ang2ipix(ra, dec));
CREATE INDEX IF NOT EXISTS idx_gaia_source_id ON gaia_matches(gaia_source_id);
CREATE INDEX IF NOT EXISTS idx_gaia_match ON gaia_matches(source_id);

-- SIMBAD catalog matches
CREATE TABLE IF NOT EXISTS simbad_matches (
    simbad_id SERIAL PRIMARY KEY,
    source_id INTEGER REFERENCES photometry_sources(source_id) ON DELETE CASCADE,
    
    -- SIMBAD identifiers  
    main_id VARCHAR(255), -- Main SIMBAD identifier
    
    -- Coordinates (from SIMBAD) with range validation
    ra DOUBLE PRECISION CHECK (ra >= 0 AND ra <= 360), -- degrees
    dec DOUBLE PRECISION CHECK (dec >= -90 AND dec <= 90), -- degrees
    
    -- Object classification
    object_type VARCHAR(50), -- SIMBAD object type
    spectral_type VARCHAR(100), -- Spectral classification
    
    -- Physical properties with validation
    distance_pc DOUBLE PRECISION CHECK (distance_pc > 0 AND distance_pc < 1e12), -- distance in parsecs
    radial_velocity_km_s DOUBLE PRECISION CHECK (radial_velocity_km_s >= -100000 AND radial_velocity_km_s <= 100000), -- km/s
    
    -- Cross-match information
    match_distance_arcsec DOUBLE PRECISION CHECK (match_distance_arcsec >= 0 AND match_distance_arcsec <= 3600), -- arcseconds
    match_timestamp TIMESTAMP DEFAULT NOW(),
    
    -- Prevent duplicate SIMBAD matches for the same source
    CONSTRAINT unique_simbad_per_source UNIQUE (source_id),
    -- Prevent duplicate matches to the same SIMBAD object (if main_id provided)
    CONSTRAINT unique_source_per_simbad UNIQUE (main_id) DEFERRABLE INITIALLY DEFERRED
);

-- Spatial index for SIMBAD sources
CREATE INDEX IF NOT EXISTS idx_simbad_coords ON simbad_matches USING btree (q3c_ang2ipix(ra, dec)) WHERE ra IS NOT NULL AND dec IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_simbad_match ON simbad_matches(source_id);
CREATE INDEX IF NOT EXISTS idx_simbad_main_id ON simbad_matches(main_id);

-- ============================================================================
-- STAGING TABLES FOR DATA VALIDATION AND QUALITY CONTROL
-- ============================================================================

-- Staging table for photometry sources (relaxed constraints for initial validation)
CREATE TABLE IF NOT EXISTS staging_photometry_sources (
    staging_id SERIAL PRIMARY KEY,
    session_id INTEGER, -- Not enforced at staging level
    
    -- Source position (no constraints for initial load)
    pixel_x DOUBLE PRECISION,
    pixel_y DOUBLE PRECISION,
    ra DOUBLE PRECISION,
    dec DOUBLE PRECISION,
    
    -- Photometry measurements
    flux DOUBLE PRECISION,
    flux_err DOUBLE PRECISION,
    background_flux DOUBLE PRECISION,
    aperture_radius DOUBLE PRECISION,
    
    -- Source properties
    fwhm DOUBLE PRECISION,
    ellipticity DOUBLE PRECISION,
    theta DOUBLE PRECISION,
    peak_counts DOUBLE PRECISION,
    signal_to_noise DOUBLE PRECISION,
    
    -- Flags
    saturated BOOLEAN DEFAULT FALSE,
    edge_source BOOLEAN DEFAULT FALSE,
    
    -- Validation metadata
    validation_status VARCHAR(20) DEFAULT 'pending', -- pending, valid, invalid, processed
    validation_errors TEXT[], -- Array of validation error messages
    staging_timestamp TIMESTAMP DEFAULT NOW()
);

-- Staging table for Gaia matches (relaxed constraints)
CREATE TABLE IF NOT EXISTS staging_gaia_matches (
    staging_id SERIAL PRIMARY KEY,
    source_staging_id INTEGER, -- References staging table, not final table
    
    gaia_source_id BIGINT,
    ra DOUBLE PRECISION,
    dec DOUBLE PRECISION,
    parallax DOUBLE PRECISION,
    parallax_error DOUBLE PRECISION,
    pmra DOUBLE PRECISION,
    pmra_error DOUBLE PRECISION,
    pmdec DOUBLE PRECISION,
    pmdec_error DOUBLE PRECISION,
    phot_g_mean_mag DOUBLE PRECISION,
    phot_g_mean_flux DOUBLE PRECISION,
    phot_g_mean_flux_error DOUBLE PRECISION,
    phot_bp_mean_mag DOUBLE PRECISION,
    phot_bp_mean_flux DOUBLE PRECISION,
    phot_bp_mean_flux_error DOUBLE PRECISION,
    phot_rp_mean_mag DOUBLE PRECISION,
    phot_rp_mean_flux DOUBLE PRECISION,
    phot_rp_mean_flux_error DOUBLE PRECISION,
    bp_rp DOUBLE PRECISION,
    phot_g_n_obs INTEGER,
    phot_bp_n_obs INTEGER,
    phot_rp_n_obs INTEGER,
    astrometric_excess_noise DOUBLE PRECISION,
    ruwe DOUBLE PRECISION,
    match_distance_arcsec DOUBLE PRECISION,
    
    -- Validation metadata
    validation_status VARCHAR(20) DEFAULT 'pending',
    validation_errors TEXT[],
    staging_timestamp TIMESTAMP DEFAULT NOW()
);

-- Staging table for SIMBAD matches (relaxed constraints)
CREATE TABLE IF NOT EXISTS staging_simbad_matches (
    staging_id SERIAL PRIMARY KEY,
    source_staging_id INTEGER, -- References staging table
    
    main_id VARCHAR(255),
    ra DOUBLE PRECISION,
    dec DOUBLE PRECISION,
    object_type VARCHAR(50),
    spectral_type VARCHAR(100),
    distance_pc DOUBLE PRECISION,
    radial_velocity_km_s DOUBLE PRECISION,
    match_distance_arcsec DOUBLE PRECISION,
    
    -- Validation metadata
    validation_status VARCHAR(20) DEFAULT 'pending',
    validation_errors TEXT[],
    staging_timestamp TIMESTAMP DEFAULT NOW()
);

-- Create indexes on staging tables for performance
CREATE INDEX IF NOT EXISTS idx_staging_photometry_status ON staging_photometry_sources(validation_status);
CREATE INDEX IF NOT EXISTS idx_staging_photometry_session ON staging_photometry_sources(session_id);
CREATE INDEX IF NOT EXISTS idx_staging_gaia_status ON staging_gaia_matches(validation_status);
CREATE INDEX IF NOT EXISTS idx_staging_simbad_status ON staging_simbad_matches(validation_status);

-- Create useful views for common queries

-- Complete source catalog view with all available data
CREATE OR REPLACE VIEW complete_source_catalog AS
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
CREATE OR REPLACE VIEW session_statistics AS
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

-- Grant permissions (adjust as needed)
-- GRANT ALL PRIVILEGES ON SCHEMA deho TO deho_user;
-- GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA deho TO deho_user;
-- GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA deho TO deho_user;

-- ============================================================================
-- DATABASE MAINTENANCE FUNCTIONS AND PROCEDURES
-- ============================================================================

-- Function to validate and move data from staging to production tables
CREATE OR REPLACE FUNCTION validate_and_promote_staging()
RETURNS TABLE(
    promoted_sources INTEGER,
    promoted_gaia INTEGER,
    promoted_simbad INTEGER,
    validation_errors TEXT
) AS $$
DECLARE
    sources_count INTEGER := 0;
    gaia_count INTEGER := 0;
    simbad_count INTEGER := 0;
    error_messages TEXT := '';
BEGIN
    -- Validate and promote photometry sources
    UPDATE staging_photometry_sources 
    SET validation_status = 'valid'
    WHERE validation_status = 'pending'
      AND pixel_x >= 0 AND pixel_y >= 0
      AND (ra IS NULL OR (ra >= 0 AND ra <= 360))
      AND (dec IS NULL OR (dec >= -90 AND dec <= 90))
      AND flux_err >= 0
      AND aperture_radius > 0
      AND fwhm > 0;
    
    -- Insert validated sources into production table
    INSERT INTO photometry_sources (
        session_id, pixel_x, pixel_y, ra, dec, flux, flux_err,
        background_flux, aperture_radius, fwhm, ellipticity, theta,
        peak_counts, signal_to_noise, saturated, edge_source
    )
    SELECT session_id, pixel_x, pixel_y, ra, dec, flux, flux_err,
           background_flux, aperture_radius, fwhm, ellipticity, theta,
           peak_counts, signal_to_noise, saturated, edge_source
    FROM staging_photometry_sources 
    WHERE validation_status = 'valid'
    ON CONFLICT (session_id, pixel_x, pixel_y) DO NOTHING;
    
    GET DIAGNOSTICS sources_count = ROW_COUNT;
    
    -- Mark processed records
    UPDATE staging_photometry_sources 
    SET validation_status = 'processed' 
    WHERE validation_status = 'valid';
    
    RETURN QUERY SELECT sources_count, gaia_count, simbad_count, error_messages;
END;
$$ LANGUAGE plpgsql;

-- Function for database maintenance (VACUUM/ANALYZE)
CREATE OR REPLACE FUNCTION maintain_database()
RETURNS TEXT AS $$
DECLARE
    maintenance_log TEXT := '';
BEGIN
    -- VACUUM ANALYZE main tables
    VACUUM ANALYZE observatory_sessions;
    maintenance_log := maintenance_log || 'VACUUM ANALYZE observatory_sessions; ';
    
    VACUUM ANALYZE photometry_sources;
    maintenance_log := maintenance_log || 'VACUUM ANALYZE photometry_sources; ';
    
    VACUUM ANALYZE gaia_matches;
    maintenance_log := maintenance_log || 'VACUUM ANALYZE gaia_matches; ';
    
    VACUUM ANALYZE simbad_matches;
    maintenance_log := maintenance_log || 'VACUUM ANALYZE simbad_matches; ';
    
    -- VACUUM staging tables (more aggressive for frequent writes)
    VACUUM FULL staging_photometry_sources;
    VACUUM FULL staging_gaia_matches;
    VACUUM FULL staging_simbad_matches;
    maintenance_log := maintenance_log || 'VACUUM FULL staging tables; ';
    
    -- Update statistics
    ANALYZE observatory_sessions;
    ANALYZE photometry_sources;
    ANALYZE gaia_matches;
    ANALYZE simbad_matches;
    maintenance_log := maintenance_log || 'ANALYZE completed; ';
    
    RETURN 'Database maintenance completed: ' || maintenance_log;
END;
$$ LANGUAGE plpgsql;

-- Function to clean old staging data
CREATE OR REPLACE FUNCTION cleanup_staging_data(days_old INTEGER DEFAULT 7)
RETURNS TEXT AS $$
DECLARE
    deleted_count INTEGER;
    cleanup_log TEXT := '';
BEGIN
    -- Clean processed staging records older than specified days
    DELETE FROM staging_photometry_sources 
    WHERE validation_status = 'processed' 
      AND staging_timestamp < NOW() - INTERVAL '1 day' * days_old;
    GET DIAGNOSTICS deleted_count = ROW_COUNT;
    cleanup_log := cleanup_log || 'Deleted ' || deleted_count || ' old photometry staging records; ';
    
    DELETE FROM staging_gaia_matches 
    WHERE validation_status = 'processed' 
      AND staging_timestamp < NOW() - INTERVAL '1 day' * days_old;
    GET DIAGNOSTICS deleted_count = ROW_COUNT;
    cleanup_log := cleanup_log || 'Deleted ' || deleted_count || ' old Gaia staging records; ';
    
    DELETE FROM staging_simbad_matches 
    WHERE validation_status = 'processed' 
      AND staging_timestamp < NOW() - INTERVAL '1 day' * days_old;
    GET DIAGNOSTICS deleted_count = ROW_COUNT;
    cleanup_log := cleanup_log || 'Deleted ' || deleted_count || ' old SIMBAD staging records; ';
    
    RETURN 'Staging cleanup completed: ' || cleanup_log;
END;
$$ LANGUAGE plpgsql;

-- Create a maintenance schedule view
CREATE OR REPLACE VIEW maintenance_status AS
SELECT 
    'observatory_sessions' as table_name,
    pg_size_pretty(pg_total_relation_size('observatory_sessions')) as size,
    (SELECT COUNT(*) FROM observatory_sessions) as row_count,
    (SELECT MAX(processing_timestamp) FROM observatory_sessions) as last_activity
UNION ALL
SELECT 
    'photometry_sources',
    pg_size_pretty(pg_total_relation_size('photometry_sources')),
    (SELECT COUNT(*) FROM photometry_sources),
    (SELECT MAX(measurement_timestamp) FROM photometry_sources)
UNION ALL
SELECT 
    'staging_photometry_sources',
    pg_size_pretty(pg_total_relation_size('staging_photometry_sources')),
    (SELECT COUNT(*) FROM staging_photometry_sources),
    (SELECT MAX(staging_timestamp) FROM staging_photometry_sources);

-- Example Q3C queries for reference:
-- 
-- -- Find sources within 1 degree of coordinates
-- SELECT * FROM complete_source_catalog 
-- WHERE q3c_radial_query(ra, dec, 180.0, 0.0, 1.0);
--
-- -- Find nearest neighbors within 5 arcsec
-- SELECT * FROM complete_source_catalog
-- WHERE q3c_radial_query(ra, dec, 180.0, 0.0, 5.0/3600.0)
-- ORDER BY q3c_dist(ra, dec, 180.0, 0.0);
--
-- -- Cross-match two catalogs
-- SELECT * FROM photometry_sources p1, photometry_sources p2  
-- WHERE p1.source_id != p2.source_id 
-- AND q3c_join(p1.ra, p1.dec, p2.ra, p2.dec, 1.0/3600.0);

COMMIT;