"""
Photometry module for stellar analysis and aperture photometry.
"""

import numpy as np
import pandas as pd
from astropy.io import fits
import astropy.units as u
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std
import matplotlib.pyplot as plt
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus
from photutils.background import Background2D, MedianBackground
from astropy.stats import sigma_clipped_stats
from tkinter import filedialog, messagebox
import logging
import threading
import time
import os
from datetime import datetime
import requests
import random
import json
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from astroquery.utils.tap.core import TapPlus
from astropy.table import Table
import tempfile
import io
from scipy.spatial import cKDTree
import re
import random
from io import BytesIO
import gzip
import subprocess
import uuid
from database import ObservatoryDatabase, DatabaseConfig

logger = logging.getLogger(__name__)


class WSLSettingsManager:
    """Manages WSL astrometry settings and persistence."""
    
    def __init__(self, config_file_path=None):
        """Initialize settings manager with config file path."""
        self.config_file = config_file_path or os.path.join(os.path.dirname(__file__), 'wsl_astrometry_config.json')
        self.default_settings = {
            'wsl_workdir': '$HOME/astrometry-work',
            'wsl_config_file': '/home/burhan/astrometry-4211-4216.cfg',
            'index_file_paths': [],
            'scale_low': 4.7,
            'scale_high': 5.3,
            'scale_units': 'degwidth',
            'downsample': 4,
            'objs': 120,
            'depth': '1-120',
            'uniformize': 10,
            'nsigma': 12,
            'code_tolerance': 0.04,
            'odds_to_solve': '1e6',
            'cpulimit': 300
        }
        
    def load_settings(self):
        """Load settings from config file, return defaults if file doesn't exist."""
        try:
            if os.path.exists(self.config_file):
                with open(self.config_file, 'r') as f:
                    settings = json.load(f)
                    # Merge with defaults to ensure all keys exist
                    merged_settings = self.default_settings.copy()
                    merged_settings.update(settings)
                    return merged_settings
        except Exception as e:
            logger.warning(f"Failed to load WSL settings from {self.config_file}: {e}")
        
        return self.default_settings.copy()
    
    def save_settings(self, settings):
        """Save settings to config file."""
        try:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.config_file), exist_ok=True)
            
            with open(self.config_file, 'w') as f:
                json.dump(settings, f, indent=2)
            logger.info(f"WSL settings saved to {self.config_file}")
            return True
        except Exception as e:
            logger.error(f"Failed to save WSL settings to {self.config_file}: {e}")
            return False
    
    def expand_wsl_path(self, wsl_path):
        """Expand WSL path variables like $HOME and ~/."""
        if not wsl_path:
            return None
            
        try:
            # Use WSL to expand the path
            cmd = ["wsl", "-d", "Ubuntu", "bash", "-c", f"echo '{wsl_path}'"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                expanded_path = result.stdout.strip()
                return expanded_path
            else:
                logger.warning(f"Failed to expand WSL path {wsl_path}: {result.stderr}")
                return wsl_path
        except Exception as e:
            logger.warning(f"Error expanding WSL path {wsl_path}: {e}")
            return wsl_path


class Photometry:
    """
    Class for performing stellar photometry analysis on astronomical images.
    
    Uses fallback strategy for astrometry: tries local solve-field first, then Nova API.
    This provides the best of both worlds - speed of local solving with reliability 
    of cloud-based solving as backup.
    
    Example:
        # Uses fallback strategy (local first, then Nova)
        photometry = Photometry()
        
        # Force only local solve-field (no fallback)
        photometry = Photometry(use_local_astrometry=True)
        photometry.set_astrometry_mode(use_local=True)  # disable fallback
    """
    
    def __init__(self, gui_reference=None, enable_database=True, use_local_astrometry=False):
        """
        Initialize the Photometry class.
        
        Args:
            gui_reference: Reference to the GUI instance for progress updates
            enable_database: Whether to enable database storage
            use_local_astrometry: If False (default), uses fallback strategy (local first, then Nova).
                                 If True, uses only local solve-field (no fallback to Nova).
        """
        self.gui_ref = gui_reference
        self.astrometry_api_key = "oecnptibpffoyzrl"
        self.astrometry_base_url = "https://nova.astrometry.net/api/"
        
        # Local astrometry configuration
        self.use_local_astrometry = use_local_astrometry
        
        # Initialize settings manager and load configuration
        self.settings_manager = WSLSettingsManager()
        self.wsl_settings = self.settings_manager.load_settings()
        
        # WSL paths will be determined dynamically based on actual user and settings
        self.wsl_user = None
        self.wsl_home = None
        self.wsl_astrometry_path = None
        self.wsl_config_path = None
        self.wsl_solve_field_path = None  # Store absolute path to solve-field
        
        # Job tracking for debugging and recovery
        self.last_subid = None
        self.last_jobid = None
        self.last_user_image_id = None
        
        # Store successful API responses for potential reuse
        self.api_data_cache = {}
        
        # HTTP session for authentication and connection reuse
        self.session = self._setup_robust_session()
        self.current_session_key = None
        
        # Initialize WSL environment if using local astrometry
        if self.use_local_astrometry:
            if not self._setup_wsl_environment():
                logger.error("Failed to initialize WSL environment - disabling local astrometry")
                self.use_local_astrometry = False
        
        # Database integration
        self.enable_database = enable_database
        self.database = None
        if enable_database:
            try:
                self.database = ObservatoryDatabase()
                if not self.database.initialize_connection_pool():
                    logger.warning("Database connection failed - continuing without database storage")
                    self.database = None
                    self.enable_database = False
            except Exception as e:
                logger.warning(f"Database initialization failed: {e} - continuing without database storage")
                self.database = None
                self.enable_database = False
        
        # SIMBAD TAP service initialization
        try:
            self.simbad_tap = TapPlus(url="https://simbad.cds.unistra.fr/simbad/sim-tap")
            self.simbad_connection_failures = 0  # Track connection failures
            logger.info("SIMBAD TAP service initialized")
        except Exception as e:
            logger.warning(f"SIMBAD TAP initialization failed: {e}")
            self.simbad_tap = None
            self.simbad_connection_failures = 0
    
    def _setup_robust_session(self):
        """
        Setup a simple HTTP session for astrometry.net API.
        
        Returns:
            requests.Session: Configured session object
        """
        session = requests.Session()
        
        # Keep it simple - just basic session configuration
        session.headers.update({
            'User-Agent': 'Mozilla/5.0 (compatible; DEHO Observatory)',
        })
        
        logger.debug("HTTP session configured for astrometry.net API")
        return session
    
    def _preflight_check_wsl(self):
        """
        Perform preflight checks for WSL environment before attempting to solve.
        
        Returns:
            bool: True if all checks pass, False otherwise
        """
        try:
            logger.info("Performing WSL preflight checks...")
            
            # Check 1: WSL availability
            wsl_check = subprocess.run(
                ["wsl", "--status"], 
                capture_output=True, text=True, timeout=10
            )
            if wsl_check.returncode != 0:
                logger.error("WSL is not available or not running")
                return False
            logger.info("✓ WSL is available")
            
            # Check 2: Ubuntu distro accessibility
            ubuntu_check = subprocess.run(
                ["wsl", "-d", "Ubuntu", "echo", "test"], 
                capture_output=True, text=True, timeout=10
            )
            if ubuntu_check.returncode != 0:
                logger.error("Ubuntu distro is not accessible")
                return False
            logger.info("✓ Ubuntu distro is accessible")
            
            # Check 3: Work directory setup and writability
            if not self.wsl_astrometry_path:
                # Resolve $HOME and set default if missing
                home_result = subprocess.run(
                    ["wsl", "-d", "Ubuntu", "echo", "$HOME"], 
                    capture_output=True, text=True, timeout=10
                )
                if home_result.returncode == 0:
                    home_dir = home_result.stdout.strip()
                    self.wsl_astrometry_path = f"{home_dir}/astrometry-work"
                else:
                    logger.error("Cannot resolve WSL $HOME directory")
                    return False
            
            # Create directory if it doesn't exist
            mkdir_result = subprocess.run(
                ["wsl", "-d", "Ubuntu", "mkdir", "-p", self.wsl_astrometry_path],
                capture_output=True, text=True, timeout=10
            )
            if mkdir_result.returncode != 0:
                logger.error(f"WSL work directory not writable: {self.wsl_astrometry_path} - mkdir failed")
                return False
            
            # Get current user for ownership
            user_result = subprocess.run(
                ["wsl", "-d", "Ubuntu", "whoami"],
                capture_output=True, text=True, timeout=10
            )
            if user_result.returncode == 0:
                current_user = user_result.stdout.strip()
                # Set ownership
                chown_result = subprocess.run(
                    ["wsl", "-d", "Ubuntu", "chown", f"{current_user}:{current_user}", self.wsl_astrometry_path],
                    capture_output=True, text=True, timeout=10
                )
                if chown_result.returncode != 0:
                    logger.warning(f"Failed to set ownership on {self.wsl_astrometry_path}")
            
            # Set permissions
            chmod_result = subprocess.run(
                ["wsl", "-d", "Ubuntu", "chmod", "775", self.wsl_astrometry_path],
                capture_output=True, text=True, timeout=10
            )
            if chmod_result.returncode != 0:
                logger.warning(f"Failed to set permissions on {self.wsl_astrometry_path}")
            
            # Test writability with temp file
            test_file = f"{self.wsl_astrometry_path}/.test_{uuid.uuid4().hex[:8]}"
            touch_result = subprocess.run(
                ["wsl", "-d", "Ubuntu", "touch", test_file],
                capture_output=True, text=True, timeout=10
            )
            if touch_result.returncode != 0:
                logger.error(f"WSL work directory not writable: {self.wsl_astrometry_path}")
                return False
            
            # Clean up test file
            subprocess.run(
                ["wsl", "-d", "Ubuntu", "rm", "-f", test_file],
                capture_output=True, text=True, timeout=10
            )
            
            logger.info(f"✓ Work directory is writable: {self.wsl_astrometry_path}")
            
            # Check 4: solve-field availability
            solve_check = subprocess.run(
                ["wsl", "-d", "Ubuntu", "which", "solve-field"],
                capture_output=True, text=True, timeout=10
            )
            if solve_check.returncode != 0:
                logger.error("solve-field command not found in Ubuntu")
                return False
            self.wsl_solve_field_path = solve_check.stdout.strip()
            logger.info(f"✓ solve-field found at: {self.wsl_solve_field_path}")
            
            # Check 5: Config file readability
            if not self.wsl_config_path:
                logger.warning("No config path set, using default")
                return True
                
            config_check = subprocess.run(
                ["wsl", "-d", "Ubuntu", "test", "-r", self.wsl_config_path],
                capture_output=True, text=True, timeout=10
            )
            if config_check.returncode != 0:
                logger.warning(f"Config file not readable: {self.wsl_config_path} - will use default")
            else:
                logger.info(f"✓ Config file readable: {self.wsl_config_path}")
            
            logger.info("All WSL preflight checks passed!")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("WSL preflight check timeout")
            return False
        except Exception as e:
            logger.error(f"WSL preflight check failed: {e}")
            return False
    
    def _setup_wsl_environment(self):
        """
        Setup WSL environment by detecting the actual user and creating work directories.
        """
        try:
            logger.info("Setting up WSL environment for astrometry solving...")
            
            # Set Ubuntu as default WSL distro
            default_result = subprocess.run(
                ["wsl", "-s", "Ubuntu"], 
                capture_output=True, text=True, timeout=10
            )
            if default_result.returncode == 0:
                logger.info("✓ Ubuntu set as default WSL distro")
            else:
                logger.warning(f"Could not set Ubuntu as default: {default_result.stderr}")
            
            # Detect WSL user
            whoami_result = subprocess.run(
                ["wsl", "-d", "Ubuntu", "whoami"], 
                capture_output=True, text=True, timeout=10
            )
            if whoami_result.returncode == 0:
                self.wsl_user = whoami_result.stdout.strip()
                logger.info(f"Detected WSL user: {self.wsl_user}")
            else:
                logger.error("Failed to detect WSL user")
                return False
            
            # Get HOME directory
            home_result = subprocess.run(
                ["wsl", "-d", "Ubuntu", "echo", "$HOME"], 
                capture_output=True, text=True, timeout=10
            )
            if home_result.returncode == 0:
                self.wsl_home = home_result.stdout.strip()
                logger.info(f"Detected WSL home: {self.wsl_home}")
            else:
                logger.error("Failed to detect WSL home directory")
                return False
            
            # Setup work directory from settings or default
            workdir_setting = self.wsl_settings.get('wsl_workdir', '$HOME/astrometry-work')
            
            # If workdir is not set or empty, use default
            if not workdir_setting:
                workdir_setting = '$HOME/astrometry-work'
            
            # Expand the workdir path (handles $HOME and ~/...)
            self.wsl_astrometry_path = self.settings_manager.expand_wsl_path(workdir_setting)
            if not self.wsl_astrometry_path:
                # Fallback if expansion failed
                self.wsl_astrometry_path = f"{self.wsl_home}/astrometry-work"
            
            logger.info(f"Using WSL work directory: {self.wsl_astrometry_path}")
            
            # Look for config in settings or common locations
            config_setting = self.wsl_settings.get('wsl_config_file')
            config_paths = []
            
            if config_setting:
                expanded_config = self.settings_manager.expand_wsl_path(config_setting)
                if expanded_config:
                    config_paths.append(expanded_config)
            
            # Add fallback paths
            config_paths.extend([
                f"{self.wsl_home}/astrometry-4211-4216.cfg",
                "/etc/astrometry.cfg"
            ])
            
            self.wsl_config_path = "/etc/astrometry.cfg"  # fallback
            for config_path in config_paths:
                check_config = subprocess.run(
                    ["wsl", "-d", "Ubuntu", "test", "-f", config_path],
                    capture_output=True, timeout=5
                )
                if check_config.returncode == 0:
                    self.wsl_config_path = config_path
                    logger.info(f"Found astrometry config: {config_path}")
                    break
            
            # Create work directory with proper permissions
            mkdir_result = subprocess.run(
                ["wsl", "-d", "Ubuntu", "mkdir", "-p", self.wsl_astrometry_path],
                capture_output=True, text=True, timeout=10
            )
            if mkdir_result.returncode != 0:
                logger.error(f"Failed to create work directory: {mkdir_result.stderr}")
                return False
            
            # Set proper permissions (775) and ownership
            chmod_result = subprocess.run(
                ["wsl", "-d", "Ubuntu", "chmod", "775", self.wsl_astrometry_path],
                capture_output=True, text=True, timeout=10
            )
            if chmod_result.returncode != 0:
                logger.warning(f"Failed to set directory permissions: {chmod_result.stderr}")
            
            # Set ownership to current user
            chown_result = subprocess.run(
                ["wsl", "-d", "Ubuntu", "chown", f"{self.wsl_user}:{self.wsl_user}", self.wsl_astrometry_path],
                capture_output=True, text=True, timeout=10
            )
            if chown_result.returncode != 0:
                logger.warning(f"Failed to set directory ownership: {chown_result.stderr}")
            
            logger.info(f"✓ Created WSL work directory: {self.wsl_astrometry_path}")
            
            # Detect solve-field path during setup
            try:
                solve_check = subprocess.run(
                    ["wsl", "-d", "Ubuntu", "which", "solve-field"],
                    capture_output=True, text=True, timeout=10
                )
                if solve_check.returncode == 0:
                    self.wsl_solve_field_path = solve_check.stdout.strip()
                    logger.info(f"✓ solve-field detected at: {self.wsl_solve_field_path}")
                else:
                    logger.warning("solve-field not found in PATH during setup")
                    self.wsl_solve_field_path = "/usr/bin/solve-field"  # fallback
            except Exception as e:
                logger.warning(f"Could not detect solve-field path: {e}")
                self.wsl_solve_field_path = "/usr/bin/solve-field"  # fallback
            
            # Update and persist settings with resolved paths
            self.wsl_settings['wsl_workdir'] = self.wsl_astrometry_path
            self.wsl_settings['wsl_config_file'] = self.wsl_config_path
            self.settings_manager.save_settings(self.wsl_settings)
            
            return True
                
        except Exception as e:
            logger.error(f"Error setting up WSL environment: {str(e)}")
            return False
    
    def _prepare_fits_for_astrometry(self, filepath):
        """
        Prepare FITS file for astrometry.net by creating an optimized version if needed.
        
        Args:
            filepath (str): Original FITS file path
            
        Returns:
            str: Path to astrometry-ready FITS file (may be same as input or new temporary file)
        """
        try:
            # First, validate if current file is already good
            if self._validate_fits_for_astrometry(filepath):
                logger.debug("FITS file is already suitable for astrometry.net")
                return filepath
            
            # If validation failed, try to create a better version using our conversion module
            logger.info("Creating astrometry.net optimized FITS file...")
            
            try:
                from conversion import Conversion
                converter = Conversion()
                
                # Create temporary optimized FITS file
                import tempfile
                temp_dir = tempfile.gettempdir()
                base_name = os.path.splitext(os.path.basename(filepath))[0]
                optimized_path = os.path.join(temp_dir, f"{base_name}_astrometry_opt.fits")
                
                # Use the astrometry-compatible conversion method
                success = converter.create_astrometry_compatible_fits(filepath, optimized_path)
                
                if success and os.path.exists(optimized_path):
                    logger.info(f"Created optimized FITS file: {optimized_path}")
                    return optimized_path
                else:
                    logger.warning("Failed to create optimized FITS file, using original")
                    return filepath
                    
            except ImportError:
                logger.warning("Conversion module not available, using original FITS file")
                return filepath
            except Exception as e:
                logger.warning(f"Error creating optimized FITS file: {e}, using original")
                return filepath
                
        except Exception as e:
            logger.error(f"Error in FITS preparation: {e}")
            return filepath
    
    def photometry(self, filepath, fwhm, threshold):
        """
        Perform aperture photometry on a FITS image.
        
        Args:
            filepath (str): Path to the FITS file
            fwhm (float): Full Width at Half Maximum for star detection
            threshold (float): Detection threshold multiplier
            
        Returns:
            tuple: (sources, photometry_table) or (None, None) if failed
        """
        try:
            logger.info(f"Starting photometry analysis on {filepath}")
            
            # Load FITS file
            with fits.open(filepath, ignore_missing_end=True) as hdulist:
                hdu = hdulist[0]
                image = hdu.data.astype(float)
            
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar2'):
                self.gui_ref.progressbar2['value'] = 25
                self.gui_ref.progressbar2.update()
            
            # Prepare image for analysis
            image -= np.median(image)
            bkg_sigma = mad_std(image)
            
            # Star detection using DAOStarFinder
            daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * bkg_sigma)
            sources = daofind(image)
            
            if sources is None or len(sources) == 0:
                messagebox.showwarning("No Stars Found", "No stellar sources detected in the image.")
                return None, None
            
            # Format source table
            for col in sources.colnames:
                sources[col].info.format = '%.8g'
            
            logger.info(f"Found {len(sources)} stellar sources")
            
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar2'):
                self.gui_ref.progressbar2['value'] = 50
                self.gui_ref.progressbar2.update()
            
            # Perform aperture photometry with background estimation
            positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
            apertures = CircularAperture(positions, r=17.)
            annulus_apertures = CircularAnnulus(positions, r_in=20., r_out=30.)
            
            # Calculate background in annulus
            annulus_masks = annulus_apertures.to_mask(method='center')
            bkg_median = []
            for mask in annulus_masks:
                annulus_data = mask.multiply(image)
                annulus_data_1d = annulus_data[mask.data > 0]
                _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
                bkg_median.append(median_sigclip)
            bkg_median = np.array(bkg_median)
            
            # Perform photometry
            phot_table = aperture_photometry(image, apertures)
            
            # Calculate background flux and subtract it
            aperture_area = apertures.area
            bkg_flux = bkg_median * aperture_area
            phot_table['bkg_flux'] = bkg_flux
            phot_table['flux'] = phot_table['aperture_sum'] - bkg_flux
            
            # Calculate flux errors (simple Poisson + background)
            # Assuming gain = 1 for now
            gain = 1.0
            effective_gain = gain
            flux_err = np.sqrt(phot_table['flux'] / effective_gain + 
                             aperture_area * bkg_sigma**2 / effective_gain)
            phot_table['flux_err'] = flux_err
            
            # Format photometry table
            for col in phot_table.colnames:
                phot_table[col].info.format = '%.8g'
            
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar2'):
                self.gui_ref.progressbar2['value'] = 75
                self.gui_ref.progressbar2.update()
            
            # Prepare FITS file for astrometry.net (may create optimized version)
            astrometry_filepath = self._prepare_fits_for_astrometry(filepath)
            
            # Perform astrometric calibration and catalog cross-matching
            wcs_solution, astrometry_data = self._solve_astrometry(astrometry_filepath, image)
            catalog_data = None
            
            # Normalize calibration data for CSV saving
            normalized_astrometry_data = self._normalize_calibration_data(astrometry_data) if astrometry_data else None
            
            if wcs_solution is not None:
                logger.info("Astrometric calibration successful - will convert pixel to sky coordinates")
                # Convert pixel coordinates to sky coordinates
                sky_coords = self._pixel_to_sky(sources, wcs_solution)
                # Cross-match with catalogs
                catalog_data = self._cross_match_catalogs(sky_coords)
                has_wcs = True
            else:
                logger.warning("No WCS solution available - proceeding without pixel-to-sky conversion")
                logger.info("Using calibration data for field metadata only")
                sky_coords = None
                catalog_data = None
                has_wcs = False
                
            # Update normalized astrometry data with WCS status
            if normalized_astrometry_data is None:
                normalized_astrometry_data = {}
            normalized_astrometry_data['has_wcs'] = has_wcs
                
            # Start CSV saving in a separate thread
            save_thread = threading.Thread(
                target=self._save_results_threaded, 
                args=(sources, phot_table, filepath, fwhm, threshold, wcs_solution, normalized_astrometry_data, sky_coords, catalog_data), 
                daemon=True
            )
            save_thread.start()
            
            # Start plot rendering in a separate thread
            plot_thread = threading.Thread(
                target=self._display_results_threaded, 
                args=(image, apertures), 
                daemon=True
            )
            plot_thread.start()
            
            # Update progress and cleanup UI
            if self.gui_ref:
                if hasattr(self.gui_ref, 'progressbar2'):
                    self.gui_ref.progressbar2['value'] = 100
                    self.gui_ref.progressbar2.update()
                    
                # Close photometry progress window
                if hasattr(self.gui_ref, 'phot_barwindow') and self.gui_ref.phot_barwindow:
                    self.gui_ref.phot_barwindow.destroy()
                    self.gui_ref.phot_barwindow = None
                    
                # Close photometry parameter window
                if hasattr(self.gui_ref, 'photometry_window') and self.gui_ref.photometry_window:
                    self.gui_ref.photometry_window.destroy()
                    self.gui_ref.photometry_window = None
                    
                if hasattr(self.gui_ref, 'statusbar'):
                    self.gui_ref.statusbar['text'] = "Ready"
                    self.gui_ref.statusbar.update()
            
            logger.info("Photometry analysis completed successfully")
            return sources, phot_table
            
        except Exception as e:
            logger.error(f"Error in photometry analysis: {str(e)}")
            messagebox.showerror('Photometry Error', f'Failed to perform photometry: {e}')
            return None, None
    
    def _save_results_threaded(self, sources, phot_table, filepath, fwhm, threshold, wcs_solution, normalized_astrometry_data, sky_coords, catalog_data):
        """
        Save photometry results to CSV files in a separate thread.
        
        Args:
            sources: Detected sources table
            phot_table: Photometry results table
            filepath: Original FITS file path
            fwhm: FWHM parameter used
            threshold: Threshold parameter used
            wcs_solution: WCS solution from astrometry.net
            normalized_astrometry_data: Normalized astrometric calibration data
            sky_coords: Sky coordinates of sources
            catalog_data: Cross-matched catalog data
        """
        try:
            # Use invoke_later to handle GUI operations from thread
            if self.gui_ref and hasattr(self.gui_ref, 'root'):
                self.gui_ref.root.after(0, lambda: self._save_results(sources, phot_table, filepath, fwhm, threshold, wcs_solution, normalized_astrometry_data, sky_coords, catalog_data))
        except Exception as e:
            logger.error(f"Error in threaded save: {str(e)}")
    
    def _save_results(self, sources, phot_table, filepath, fwhm, threshold, wcs_solution, normalized_astrometry_data, sky_coords, catalog_data):
        """
        Save merged photometry results to a single CSV file with metadata.
        
        Args:
            sources: Detected sources table
            phot_table: Photometry results table
            filepath: Original FITS file path
            fwhm: FWHM parameter used
            threshold: Threshold parameter used
            wcs_solution: WCS solution from astrometry.net
            normalized_astrometry_data: Normalized astrometric calibration data
            sky_coords: Sky coordinates of sources
            catalog_data: Cross-matched catalog data
        """
        try:
            # Get save location for merged CSV
            filedest = filedialog.asksaveasfilename(
                title="Save photometry results", 
                defaultextension=".csv",
                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
            )
            
            if not filedest:
                return
            
            # Create merged dataframe
            merged_data = []
            
            for i in range(len(sources)):
                row = {
                    # Metadata
                    'filename': os.path.basename(filepath),
                    'fwhm': fwhm,
                    'threshold': threshold,
                    
                    # Astrometric calibration metadata (normalized)
                    'center_ra': normalized_astrometry_data.get('center_ra', np.nan) if normalized_astrometry_data else np.nan,
                    'center_dec': normalized_astrometry_data.get('center_dec', np.nan) if normalized_astrometry_data else np.nan,
                    'field_width': normalized_astrometry_data.get('field_width_arcmin', np.nan) if normalized_astrometry_data else np.nan,
                    'field_height': normalized_astrometry_data.get('field_height_arcmin', np.nan) if normalized_astrometry_data else np.nan,
                    'pixel_scale': normalized_astrometry_data.get('pixscale', np.nan) if normalized_astrometry_data else np.nan,
                    'orientation': normalized_astrometry_data.get('orientation', np.nan) if normalized_astrometry_data else np.nan,
                    'parity': normalized_astrometry_data.get('parity', np.nan) if normalized_astrometry_data else np.nan,
                    
                    # Source detection data
                    'id': sources['id'][i],
                    'xcentroid': sources['xcentroid'][i],
                    'ycentroid': sources['ycentroid'][i],
                    'sharpness': sources['sharpness'][i],
                    'roundness1': sources['roundness1'][i],
                    'roundness2': sources['roundness2'][i],
                    'npix': sources['npix'][i],
                    'sky': sources['sky'][i],
                    'peak': sources['peak'][i],
                    'flux_detected': sources['flux'][i],
                    'mag': sources['mag'][i],
                    
                    # Sky coordinates (from WCS if available)
                    'ra': sky_coords[i].ra.degree if sky_coords else np.nan,
                    'dec': sky_coords[i].dec.degree if sky_coords else np.nan,
                    'ra_hms': sky_coords[i].ra.to_string(unit=u.hour, sep=':', precision=2) if sky_coords else '',
                    'dec_dms': sky_coords[i].dec.to_string(sep=':', precision=1) if sky_coords else '',
                    
                    # WCS availability status
                    'has_wcs': normalized_astrometry_data.get('has_wcs', False) if normalized_astrometry_data else False,
                    'has_astrometry': normalized_astrometry_data is not None,
                    
                    # Photometry data
                    'aperture_sum': phot_table['aperture_sum'][i],
                    'bkg_flux': phot_table['bkg_flux'][i],
                    'flux': phot_table['flux'][i],
                    'flux_err': phot_table['flux_err'][i],
                    
                    # Derived quantities
                    'snr': phot_table['flux'][i] / phot_table['flux_err'][i] if phot_table['flux_err'][i] > 0 else 0,
                    'mag_calibrated': -2.5 * np.log10(phot_table['flux'][i]) if phot_table['flux'][i] > 0 else np.nan,
                }
                
                # Add catalog cross-match data if available
                if catalog_data and i < len(catalog_data['gaia']):
                    gaia_match = catalog_data['gaia'][i]
                    simbad_match = catalog_data['simbad'][i]
                    
                    row.update({
                        # Gaia DR3 data
                        'gaia_source_id': gaia_match.get('source_id', ''),
                        'gaia_ra': gaia_match.get('ra', np.nan),
                        'gaia_dec': gaia_match.get('dec', np.nan),
                        'phot_g_mean_mag': gaia_match.get('phot_g_mean_mag', np.nan),
                        'phot_bp_mean_mag': gaia_match.get('phot_bp_mean_mag', np.nan),
                        'phot_rp_mean_mag': gaia_match.get('phot_rp_mean_mag', np.nan),
                        'bp_rp': gaia_match.get('bp_rp', np.nan),
                        'parallax': gaia_match.get('parallax', np.nan),
                        'pmra': gaia_match.get('pmra', np.nan),
                        'pmdec': gaia_match.get('pmdec', np.nan),
                        
                        # SIMBAD data (using X-Match field names)
                        'simbad_main_id': simbad_match.get('simbad_main_id', ''),
                        'otype': simbad_match.get('otype', ''),
                        'sp_type': simbad_match.get('sp_type', ''),
                        'rv_value': simbad_match.get('rvz_radvel', np.nan),
                        'distance_result': simbad_match.get('distance_result', np.nan),
                    })
                else:
                    # Add empty catalog columns if no match
                    row.update({
                        'gaia_source_id': '', 'gaia_ra': np.nan, 'gaia_dec': np.nan,
                        'phot_g_mean_mag': np.nan, 'phot_bp_mean_mag': np.nan, 'phot_rp_mean_mag': np.nan,
                        'bp_rp': np.nan, 'parallax': np.nan, 'pmra': np.nan, 'pmdec': np.nan,
                        'simbad_main_id': '', 'otype': '', 'sp_type': '', 'rv_value': np.nan, 'distance_result': np.nan
                    })
                
                merged_data.append(row)
            
            # Convert to DataFrame and save CSV
            df = pd.DataFrame(merged_data)
            df.to_csv(filedest, index=False, float_format='%.8g')
            
            logger.info(f"Merged photometry results saved to {filedest}")
            
            # Store in database if enabled
            if self.enable_database and self.database:
                try:
                    self._store_results_in_database(merged_data, filepath, fwhm, threshold, normalized_astrometry_data, catalog_data)
                except Exception as db_error:
                    logger.error(f"Database storage failed: {db_error}")
                    # Continue with CSV-only save
            
            messagebox.showinfo("Save Complete", f"Photometry results saved successfully!\\n{os.path.basename(filedest)}")
            
        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")
            messagebox.showerror('Save Error', f'Failed to save results: {e}')
    
    def _store_results_in_database(self, merged_data, filepath, fwhm, threshold, normalized_astrometry_data, catalog_data):
        """
        Store photometry results in PostgreSQL database with proper unit conversions.
        
        Args:
            merged_data: List of dictionaries containing all photometry data
            filepath: Original FITS file path
            fwhm: FWHM parameter used
            threshold: Threshold parameter used
            normalized_astrometry_data: Normalized astrometric calibration data
            catalog_data: Cross-matched catalog data
        """
        try:
            logger.info("Starting database storage of photometry results")
            
            # Extract FITS header information for session metadata
            with fits.open(filepath, ignore_missing_end=True) as hdulist:
                header = hdulist[0].header
            
            # Create observatory session record
            session_data = {
                'session_name': f"Photometry_{os.path.splitext(os.path.basename(filepath))[0]}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
                'fits_file_path': filepath,
                'image_center_ra': normalized_astrometry_data.get('center_ra') if normalized_astrometry_data else None,
                'image_center_dec': normalized_astrometry_data.get('center_dec') if normalized_astrometry_data else None,
                'field_of_view_deg': normalized_astrometry_data.get('field_width_arcmin', 0) / 60.0 if normalized_astrometry_data and normalized_astrometry_data.get('field_width_arcmin') else None,
                'pixel_scale_arcsec': normalized_astrometry_data.get('pixscale') if normalized_astrometry_data else None,
                'observation_date': self._parse_fits_date(header.get('DATE-OBS')) if header.get('DATE-OBS') else None,
                'camera_info': header.get('CAMERA', header.get('INSTRUME', 'Unknown')),
                'exposure_time': header.get('EXPTIME'),
                'iso_value': header.get('ISO'),
                'aperture_fnum': header.get('APERTURE'),
                'focal_length_mm': header.get('FOC_LEN'),
                'fwhm_pixels': fwhm,
                'detection_threshold': threshold,
                'total_sources_detected': len(merged_data),
                'astrometry_solved': normalized_astrometry_data.get('has_astrometry', False) if normalized_astrometry_data else False,
                'wcs_available': normalized_astrometry_data.get('has_wcs', False) if normalized_astrometry_data else False,
                'notes': f"Processed with DEHO v2.0. Astrometry: {'Yes' if normalized_astrometry_data else 'No'}"
            }
            
            session_id = self.database.create_session(session_data)
            if not session_id:
                logger.error("Failed to create observatory session record")
                return
            
            logger.info(f"Created database session {session_id}")
            
            # Prepare photometry sources data
            sources_data = []
            gaia_matches_data = []
            simbad_matches_data = []
            
            for i, row in enumerate(merged_data):
                # Photometry source data with proper unit conversions
                source_data = {
                    'pixel_x': row.get('xcentroid'),
                    'pixel_y': row.get('ycentroid'),
                    'ra': row.get('ra') if not np.isnan(row.get('ra', np.nan)) else None,  # degrees
                    'dec': row.get('dec') if not np.isnan(row.get('dec', np.nan)) else None,  # degrees
                    'flux': row.get('flux'),  # counts
                    'flux_err': row.get('flux_err'),  # counts
                    'background_flux': row.get('bkg_flux'),  # counts
                    'aperture_radius': 17.0,  # pixels (fixed aperture size)
                    'fwhm': fwhm,  # pixels
                    'ellipticity': None,  # Not calculated in current implementation
                    'theta': None,  # Not calculated in current implementation
                    'peak_counts': row.get('peak'),
                    'signal_to_noise': row.get('snr'),
                    'saturated': False,  # Would need saturation detection
                    'edge_source': False,  # Would need edge detection
                }
                sources_data.append(source_data)
                
                # Prepare Gaia match data if available
                if row.get('gaia_source_id'):
                    gaia_data = {
                        'source_id': None,  # Will be set after source insertion
                        'gaia_source_id': int(row.get('gaia_source_id', 0)) if row.get('gaia_source_id') else None,
                        'ra': row.get('gaia_ra'),  # degrees
                        'dec': row.get('gaia_dec'),  # degrees  
                        'parallax': row.get('parallax'),  # milliarcseconds (Gaia native unit)
                        'parallax_error': None,  # Not available in current data
                        'pmra': row.get('pmra'),  # mas/year
                        'pmra_error': None,  # Not available
                        'pmdec': row.get('pmdec'),  # mas/year
                        'pmdec_error': None,  # Not available
                        'phot_g_mean_mag': row.get('phot_g_mean_mag'),  # mag
                        'phot_g_mean_flux': None,  # Not available
                        'phot_g_mean_flux_error': None,  # Not available
                        'phot_bp_mean_mag': row.get('phot_bp_mean_mag'),  # mag
                        'phot_bp_mean_flux': None,  # Not available
                        'phot_bp_mean_flux_error': None,  # Not available
                        'phot_rp_mean_mag': row.get('phot_rp_mean_mag'),  # mag
                        'phot_rp_mean_flux': None,  # Not available
                        'phot_rp_mean_flux_error': None,  # Not available
                        'bp_rp': row.get('bp_rp'),  # mag
                        'phot_g_n_obs': None,  # Not available
                        'phot_bp_n_obs': None,  # Not available
                        'phot_rp_n_obs': None,  # Not available
                        'astrometric_excess_noise': None,  # Not available
                        'ruwe': None,  # Not available
                        'match_distance_arcsec': None,  # Would need cross-match distance
                    }
                    gaia_matches_data.append((i, gaia_data))
                
                # Prepare SIMBAD match data if available  
                if row.get('simbad_main_id'):
                    simbad_data = {
                        'source_id': None,  # Will be set after source insertion
                        'main_id': row.get('simbad_main_id'),
                        'ra': row.get('ra') if not np.isnan(row.get('ra', np.nan)) else None,  # degrees
                        'dec': row.get('dec') if not np.isnan(row.get('dec', np.nan)) else None,  # degrees
                        'object_type': row.get('otype'),
                        'spectral_type': row.get('sp_type'),
                        'distance_pc': row.get('distance_result'),  # Convert from various units to pc if needed
                        'radial_velocity_km_s': row.get('rv_value'),  # km/s
                        'match_distance_arcsec': None,  # Would need cross-match distance
                    }
                    simbad_matches_data.append((i, simbad_data))
            
            # Insert photometry sources
            source_ids = self.database.insert_photometry_sources(session_id, sources_data)
            logger.info(f"Inserted {len(source_ids)} photometry sources")
            
            # Insert Gaia matches with correct source IDs
            if gaia_matches_data:
                gaia_final_data = []
                for source_index, gaia_data in gaia_matches_data:
                    if source_index < len(source_ids):
                        gaia_data['source_id'] = source_ids[source_index]
                        gaia_final_data.append(gaia_data)
                
                gaia_count = self.database.insert_gaia_matches(gaia_final_data)
                logger.info(f"Inserted {gaia_count} Gaia matches")
            
            # Insert SIMBAD matches with correct source IDs
            if simbad_matches_data:
                simbad_final_data = []
                for source_index, simbad_data in simbad_matches_data:
                    if source_index < len(source_ids):
                        simbad_data['source_id'] = source_ids[source_index]
                        simbad_final_data.append(simbad_data)
                
                simbad_count = self.database.insert_simbad_matches(simbad_final_data)
                logger.info(f"Inserted {simbad_count} SIMBAD matches")
            
            logger.info(f"Successfully stored complete photometry session {session_id} in database")
            
        except Exception as e:
            logger.error(f"Error storing results in database: {str(e)}")
            raise
    
    def _parse_fits_date(self, date_string):
        """
        Parse FITS DATE-OBS string to datetime object.
        
        Args:
            date_string: FITS date string (various formats)
            
        Returns:
            datetime: Parsed datetime or None if parsing fails
        """
        if not date_string:
            return None
            
        try:
            # Handle various FITS date formats
            if 'T' in date_string:
                # ISO format: 2023-01-15T12:34:56.789
                return datetime.fromisoformat(date_string.replace('Z', '+00:00'))
            else:
                # Date only: 2023-01-15
                return datetime.strptime(date_string, '%Y-%m-%d')
        except Exception as e:
            logger.debug(f"Could not parse FITS date '{date_string}': {e}")
            return None
    
    def _display_results_threaded(self, image, apertures):
        """
        Display results in a separate thread to avoid blocking GUI.
        
        Args:
            image: The astronomical image array
            apertures: CircularAperture objects for detected sources
        """
        try:
            # Add a small delay to ensure CSV saving dialog appears first
            time.sleep(0.5)
            self._display_results(image, apertures)
        except Exception as e:
            logger.error(f"Error in threaded display: {str(e)}")
    
    def _display_results(self, image, apertures):
        """
        Display the image with apertures overlaid.
        
        Args:
            image: The astronomical image array
            apertures: CircularAperture objects for detected sources
        """
        try:
            plt.figure(figsize=(12, 10))
            plt.imshow(image, cmap='gray_r', origin='lower')
            apertures.plot(color='red', lw=1.5, alpha=0.7)
            plt.title('Aperture Photometry Results')
            plt.xlabel('X (pixels)')
            plt.ylabel('Y (pixels)')
            plt.colorbar(label='Counts')
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            logger.error(f"Error displaying results: {str(e)}")
            messagebox.showerror('Display Error', f'Failed to display results: {e}')
    
    def _solve_astrometry(self, filepath, image):
        """
        Solve astrometry with fallback strategy: try local solve-field first, then Nova API.
        
        Args:
            filepath: Path to FITS file
            image: Image data array
            
        Returns:
            tuple: (WCS solution, astrometry_data dict) or (None, None)
        """
        try:
            logger.info("="*60)
            
            if self.use_local_astrometry:
                # Local-only mode (no fallback)
                logger.info("STARTING ASTROMETRIC CALIBRATION WITH LOCAL SOLVE-FIELD ONLY")
                logger.info("="*60)
                logger.info(f"Input file: {filepath}")
                logger.info(f"File size: {os.path.getsize(filepath)} bytes")
                
                # Enforce preflight check before local solve
                if not self._preflight_check_wsl():
                    logger.error("WSL preflight checks failed - local-only mode cannot proceed")
                    return None, None
                
                return self._solve_local_astrometry(filepath, image)
            else:
                # Fallback strategy: try local first, then Nova
                logger.info("STARTING ASTROMETRIC CALIBRATION WITH FALLBACK STRATEGY")
                logger.info("="*60)
                logger.info(f"Input file: {filepath}")
                logger.info(f"File size: {os.path.getsize(filepath)} bytes")
                
                # First attempt: Local solve-field
                logger.info("ATTEMPT 1: Trying local solve-field...")
                
                # Enforce preflight check before local solve
                if not self._preflight_check_wsl():
                    logger.warning("WSL preflight checks failed - skipping to Nova fallback")
                    wcs, astrometry_data = None, None
                else:
                    wcs, astrometry_data = self._solve_local_astrometry(filepath, image)
                
                if wcs is not None and astrometry_data is not None:
                    logger.info("SUCCESS: Local solve-field completed successfully")
                    return wcs, astrometry_data
                
                # Second attempt: Nova API fallback
                logger.info("ATTEMPT 2: Local solve-field failed, falling back to Nova API...")
                return self._solve_nova_astrometry(filepath, image)
            
        except Exception as e:
            logger.error(f"Error in astrometry solving: {str(e)}")
            return None, None
    
    def _solve_nova_astrometry(self, filepath, image):
        """
        Solve astrometry using Nova.astrometry.net API.
        
        Args:
            filepath: Path to FITS file
            image: Image data array
            
        Returns:
            tuple: (WCS solution, astrometry_data dict) or (None, None)
        """
        try:
            logger.info("STARTING ASTROMETRIC CALIBRATION WITH ASTROMETRY.NET API")
            
            # Update progress
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar2'):
                self.gui_ref.progressbar2['value'] = 80
                self.gui_ref.progressbar2.update()
            
            # Step 1: Login to astrometry.net
            logger.info("STEP 1: Logging into Astrometry.net...")
            session_key = self._astrometry_login()
            if not session_key:
                logger.error("STEP 1 FAILED: Could not login to Astrometry.net")
                return None, None
            logger.info(f"STEP 1 SUCCESS: Session key obtained")
            
            # Step 2: Submit image for solving
            logger.info("STEP 2: Submitting image for plate solving...")
            submission_id = self._submit_image_for_solving(filepath, session_key)
            if not submission_id:
                logger.error("STEP 2 FAILED: Could not submit image to Astrometry.net")
                return None, None
            logger.info(f"STEP 2 SUCCESS: Submission ID: {submission_id}")
            
            # Step 3: Wait for solution
            logger.info("STEP 3: Waiting for astrometric solution...")
            job_id = self._wait_for_solution(submission_id, session_key)
            if not job_id:
                logger.error("STEP 3 FAILED: Astrometry solving failed or timed out")
                return None, None
            logger.info(f"STEP 3 SUCCESS: Job completed, Job ID: {job_id}")
            
            # Step 4: Get calibration data and WCS
            logger.info("STEP 4: Retrieving calibration data and WCS solution...")
            
            # Add cool-down period for WCS file generation
            cooldown_time = random.uniform(2.0, 5.0)
            logger.info(f"Waiting {cooldown_time:.1f}s for WCS file generation...")
            time.sleep(cooldown_time)
            
            calibration_data = self._get_calibration_data(job_id)
            
            # Early caching - store data even if WCS retrieval fails
            if calibration_data:
                self._cache_early_data(job_id, calibration_data)
            
            wcs_solution = self._get_wcs_solution(job_id)
            
            # Continue even if WCS fails - we can still use calibration data
            if wcs_solution is None:
                logger.warning("STEP 4 PARTIAL: Could not retrieve WCS solution, but continuing with calibration data")
                if not calibration_data:
                    logger.error("STEP 4 FAILED: No WCS solution and no calibration data available")
                    return None, None
                else:
                    logger.info("STEP 4 SUCCESS: Calibration data available for metadata and analysis")
            else:
                logger.info("STEP 4 SUCCESS: Both WCS solution and calibration data retrieved")
            logger.info("="*60)
            logger.info("ASTROMETRIC CALIBRATION COMPLETED SUCCESSFULLY!")
            logger.info("="*60)
            
            # Cache successful API data for potential future use
            if job_id not in self.api_data_cache:
                self.api_data_cache[job_id] = {
                    'calibration_data': calibration_data,
                    'job_id': job_id,
                    'submission_id': self.last_subid,
                    'user_image_id': self.last_user_image_id,
                    'timestamp': time.time(),
                    'wcs_solution': None
                }
            
            # Update with WCS solution if successful
            self.api_data_cache[job_id]['wcs_solution'] = wcs_solution
            self.api_data_cache[job_id]['has_wcs'] = wcs_solution is not None
            
            return wcs_solution, calibration_data
            
        except Exception as e:
            logger.error(f"Error in Nova astrometric calibration: {str(e)}")
            return None, None
    
    def _astrometry_login(self):
        """Login to Astrometry.net API."""
        try:
            login_url = f"{self.astrometry_base_url}login"
            login_data = {'request-json': json.dumps({"apikey": self.astrometry_api_key})}
            
            logger.debug(f"Attempting login to: {login_url}")
            logger.debug(f"Login data: {login_data}")
            
            response = self.session.post(login_url, data=login_data, timeout=60)
            logger.debug(f"Login response status code: {response.status_code}")
            logger.debug(f"Login response headers: {dict(response.headers)}")
            logger.debug(f"Login response text: {response.text}")
            
            response.raise_for_status()
            
            result = response.json()
            logger.debug(f"Login result: {result}")
            
            if result.get('status') == 'success':
                session_key = result.get('session')
                self.current_session_key = session_key
                # Store session in cookies for automatic authentication
                self.session.cookies.set('astrometry_session', session_key)
                logger.info(f"Astrometry.net login successful, session: {session_key}")
                return session_key
            else:
                logger.error(f"Astrometry.net login failed: {result}")
                return None
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Network error during Astrometry.net login: {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Error logging into Astrometry.net: {str(e)}")
            return None
    
    def parse_submission_body(self, text):
        """
        Parse submission response body with improved validation and fallback strategies.
        
        Args:
            text (str): Response body text
            
        Returns:
            dict|None: Parsed data with jobs, user_images, etc. or None if unparseable
        """
        if not text or not text.strip():
            logger.debug("Empty or whitespace-only response text")
            return None
        
        text = text.strip()
        
        # Check if this looks like JSON before attempting to parse
        is_likely_json = (text.startswith('{') and text.endswith('}')) or (text.startswith('[') and text.endswith(']'))
        
        # Strategy 1: Try standard JSON parsing if it looks like JSON
        if is_likely_json:
            try:
                data = json.loads(text)
                logger.debug("Successfully parsed response as valid JSON")
                
                # Validate that we have the expected structure
                if isinstance(data, dict):
                    has_valid_structure = any(key in data for key in ['jobs', 'user_images', 'status', 'subid'])
                    if has_valid_structure:
                        logger.debug(f"JSON contains expected keys: {list(data.keys())}")
                        return data
                    else:
                        logger.warning(f"JSON parsed but missing expected keys. Available keys: {list(data.keys())}")
                
                return data
                
            except (json.JSONDecodeError, ValueError) as e:
                logger.debug(f"JSON parsing failed despite JSON-like structure: {str(e)}")
        else:
            logger.debug("Response text does not appear to be JSON format")
        
        # Strategy 2: Regex extraction for key fields (for malformed JSON or plain text)
        try:
            result = {}
            
            # Extract jobs array: "jobs": [123, 456] or "jobs":[123]
            jobs_match = re.search(r'"jobs"\s*:\s*\[\s*([0-9,\s]+)\s*\]', text)
            if jobs_match:
                job_ids = [int(x.strip()) for x in jobs_match.group(1).split(',') if x.strip().isdigit()]
                result['jobs'] = job_ids
                logger.info(f"Extracted jobs from text using regex: {job_ids}")
            
            # Extract user_images array: "user_images": [789]
            user_images_match = re.search(r'"user_images"\s*:\s*\[\s*([0-9,\s]+)\s*\]', text)
            if user_images_match:
                user_image_ids = [int(x.strip()) for x in user_images_match.group(1).split(',') if x.strip().isdigit()]
                result['user_images'] = user_image_ids
                logger.info(f"Extracted user_images from text using regex: {user_image_ids}")
            
            # Extract status: "status": "success" 
            status_match = re.search(r'"status"\s*:\s*"([^"]+)"', text)
            if status_match:
                result['status'] = status_match.group(1)
                logger.debug(f"Extracted status: {result['status']}")
            
            # Extract submission ID: "subid": 12345
            subid_match = re.search(r'"subid"\s*:\s*(\d+)', text)
            if subid_match:
                result['subid'] = int(subid_match.group(1))
                logger.debug(f"Extracted subid: {result['subid']}")
            
            # Extract processing timestamps
            processing_started_match = re.search(r'"processing_started"\s*:\s*"([^"]+)"', text)
            if processing_started_match:
                result['processing_started'] = processing_started_match.group(1)
            
            processing_finished_match = re.search(r'"processing_finished"\s*:\s*"([^"]+)"', text)
            if processing_finished_match:
                result['processing_finished'] = processing_finished_match.group(1)
            
            if result:
                logger.info(f"Successfully extracted data via regex fallback: {result}")
                return result
            else:
                logger.warning("Regex extraction found no recognizable patterns")
                logger.debug(f"Failed text content (first 300 chars): {text[:300]}")
                return None
                
        except Exception as e:
            logger.error(f"Error in regex extraction: {str(e)}")
            logger.debug(f"Failed text content (first 300 chars): {text[:300]}")
            return None
    
    def _submit_image_for_solving(self, filepath, session_key):
        """Submit FITS image to Astrometry.net for solving with improved reliability."""
        try:
            submit_url = f"{self.astrometry_base_url}upload"
            
            # Check file existence and size
            if not os.path.exists(filepath):
                logger.error(f"FITS file does not exist: {filepath}")
                return None
                
            file_size = os.path.getsize(filepath)
            logger.info(f"Submitting FITS file: {filepath} (size: {file_size} bytes)")
            
            # Validate file size (astrometry.net has limits)
            max_size = 100 * 1024 * 1024  # 100MB limit
            if file_size > max_size:
                logger.error(f"File too large: {file_size} bytes (max: {max_size} bytes)")
                return None
                
            # Validate FITS file before upload (warnings only)
            try:
                is_valid = self._validate_fits_for_astrometry(filepath)
                if not is_valid:
                    logger.warning("FITS file validation failed, but proceeding anyway")
                else:
                    logger.debug("FITS file validation passed")
            except Exception as e:
                logger.warning(f"FITS validation error (proceeding anyway): {e}")
                
            logger.debug(f"Submit URL: {submit_url}")
            
            # Prepare the submission data - keep it simple and minimal like web interface
            submission_data = {
                "session": session_key,
                "allow_commercial_use": "d",
                "allow_modifications": "d", 
                "publicly_visible": "y",
                "scale_units": "degwidth",
                "scale_type": "ul",
                "scale_lower": 0.1,
                "scale_upper": 10
            }
            
            logger.debug(f"Submission data: {submission_data}")
            
            # Retry logic with exponential backoff and rate limiting
            max_retries = 5  # Increased retries
            base_delay = 3   # Longer initial delay
            
            for attempt in range(max_retries):
                try:
                    # Add rate limiting delay
                    if attempt > 0:
                        delay = base_delay * (2 ** (attempt - 1)) + random.uniform(0, 2)
                        logger.info(f"Rate limiting delay: {delay:.1f}s before attempt {attempt + 1}")
                        time.sleep(delay)
                        
                    logger.info(f"Upload attempt {attempt + 1}/{max_retries}")
                    
                    # Use the proven working format from debug testing
                    with open(filepath, 'rb') as f:
                        files = {
                            'file': (os.path.basename(filepath), f, 'application/fits'),
                            'request-json': (None, json.dumps(submission_data), 'application/json')
                        }
                        
                        logger.info("Starting file upload to Astrometry.net...")
                        response = self.session.post(
                            submit_url, 
                            files=files, 
                            timeout=600  # 10 minutes timeout
                        )
                    
                    logger.debug(f"Upload response status code: {response.status_code}")
                    logger.debug(f"Upload response headers: {dict(response.headers)}")
                    logger.debug(f"Upload response text: {response.text}")
                    
                    response.raise_for_status()
                    
                    result = response.json()
                    logger.debug(f"Upload result: {result}")
                    
                    if result.get('status') == 'success':
                        subid = result.get('subid')
                        # Store submission ID immediately
                        self.last_subid = subid
                        
                        # Also store any jobs or user_images if already available
                        if 'jobs' in result and result['jobs']:
                            self.last_jobid = result['jobs'][0]
                            logger.info(f"Job ID available immediately: {self.last_jobid}")
                        
                        if 'user_images' in result and result['user_images']:
                            self.last_user_image_id = result['user_images'][0]
                            logger.info(f"User image ID available immediately: {self.last_user_image_id}")
                        
                        logger.info(f"Image submitted successfully, submission ID: {subid}")
                        return subid
                    elif result.get('status') == 'error':
                        error_msg = result.get('errormessage', 'Unknown error')
                        logger.error(f"API error: {error_msg}")
                        
                        # Don't retry for certain errors
                        if 'no json' in error_msg.lower() or 'invalid' in error_msg.lower():
                            logger.error("Non-retryable error, aborting")
                            return None
                        
                        # Retry for other errors
                        if attempt < max_retries - 1:
                            delay = base_delay * (2 ** attempt)  # Exponential backoff
                            logger.warning(f"Retrying in {delay} seconds...")
                            time.sleep(delay)
                            continue
                        else:
                            logger.error(f"Unexpected response: {result}")
                            if attempt < max_retries - 1:
                                delay = base_delay * (2 ** attempt)
                                logger.warning(f"Retrying in {delay} seconds...")
                                time.sleep(delay)
                                continue
                        
                        return None
                        
                except requests.exceptions.RequestException as e:
                    logger.error(f"Network error during upload attempt {attempt + 1}: {str(e)}")
                    if attempt < max_retries - 1:
                        delay = base_delay * (2 ** attempt)
                        logger.warning(f"Retrying in {delay} seconds...")
                        time.sleep(delay)
                    continue
                except Exception as e:
                    logger.error(f"Unexpected error during upload: {str(e)}")
                    if attempt < max_retries - 1:
                        delay = base_delay * (2 ** attempt)  
                        logger.warning(f"Retrying in {delay} seconds...")
                        time.sleep(delay)
                    continue
            
            logger.error(f"All {max_retries} upload attempts failed")
            return None
                    
        except Exception as e:
            logger.error(f"Fatal error during image submission: {str(e)}")
            return None
    
    def _validate_fits_for_astrometry(self, filepath):
        """
        Validate FITS file for astrometry.net compatibility.
        
        Args:
            filepath (str): Path to FITS file
            
        Returns:
            bool: True if file is valid for astrometry.net
        """
        try:
            logger.debug(f"Validating FITS file for astrometry.net: {filepath}")
            
            with fits.open(filepath) as hdul:
                if len(hdul) == 0:
                    logger.error("FITS file has no HDUs")
                    return False
                    
                header = hdul[0].header
                data = hdul[0].data
                
                # Check essential headers
                required_headers = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2']
                for req_header in required_headers:
                    if req_header not in header:
                        logger.error(f"Missing required header: {req_header}")
                        return False
                
                # Check image dimensions
                if header['NAXIS'] != 2:
                    logger.error(f"Image must be 2D, found NAXIS={header['NAXIS']}")
                    return False
                    
                width = header['NAXIS1']
                height = header['NAXIS2']
                
                # Check minimum size (astrometry.net needs sufficient resolution)
                if width < 100 or height < 100:
                    logger.error(f"Image too small: {width}x{height} (minimum 100x100)")
                    return False
                    
                # Check maximum size (performance considerations)
                max_pixels = 8000 * 8000
                if width * height > max_pixels:
                    logger.warning(f"Large image: {width}x{height} may take longer to solve")
                    
                # Check data type
                if data is None:
                    logger.error("FITS file contains no image data")
                    return False
                    
                # Check for reasonable data range
                if data.size == 0:
                    logger.error("FITS file contains empty image data")
                    return False
                    
                data_min, data_max = np.nanmin(data), np.nanmax(data)
                if np.isnan(data_min) or np.isnan(data_max):
                    logger.error("FITS file contains all NaN values")
                    return False
                    
                if data_min == data_max:
                    logger.error("FITS file contains constant values (no variation)")
                    return False
                    
                # Check for conflicting WCS headers (can confuse astrometry.net)
                problematic_wcs = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 
                                 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
                existing_wcs = [h for h in problematic_wcs if h in header]
                
                if existing_wcs:
                    logger.warning(f"FITS file contains WCS headers that may interfere: {existing_wcs}")
                    logger.warning("Consider using conversion.py to create a clean FITS file")
                
                logger.debug(f"FITS validation passed: {width}x{height}, range: {data_min:.1f} to {data_max:.1f}")
                return True
                
        except Exception as e:
            logger.error(f"FITS validation error: {str(e)}")
            return False
    
    def _wait_for_solution(self, submission_id, session_key, max_wait=1200):  # Increased to 20 minutes
        """Wait for astrometry solution to complete with intelligent polling and robust parsing."""
        try:
            start_time = time.time()
            self.last_subid = submission_id
            
            logger.info(f"Waiting for solution, submission ID: {submission_id}")
            logger.info(f"Maximum wait time: {max_wait} seconds ({max_wait/60:.1f} minutes)")
            
            # Initial delay - wait for job to be queued
            initial_delay = random.randint(25, 35)  # 25-35s with jitter
            logger.info(f"Waiting {initial_delay} seconds for job to be queued...")
            time.sleep(initial_delay)
            
            check_count = 0
            last_status = None
            consecutive_failures = 0
            max_consecutive_failures = 8  # More lenient - only for true dead ends
            discovered_job_id = None
            
            while time.time() - start_time < max_wait:
                check_count += 1
                elapsed_time = time.time() - start_time
                
                logger.info(f"Status check #{check_count}, elapsed: {elapsed_time:.1f}s")
                
                # If we have a job_id, switch to job polling exclusively
                if discovered_job_id:
                    return self._poll_job_until_complete(discovered_job_id, start_time, max_wait, last_status)
                
                # Poll submission status to discover job_id
                submission_data = self._get_submission_status_robust(submission_id)
                
                if submission_data is None:
                    consecutive_failures += 1
                    if consecutive_failures >= max_consecutive_failures:
                        logger.error(f"Too many consecutive failures ({consecutive_failures}), aborting")
                        return None
                    
                    # Exponential backoff for failures
                    backoff_delay = min(30, 5 * (2 ** min(consecutive_failures, 4)))
                    logger.warning(f"Failed to get submission status, waiting {backoff_delay}s...")
                    time.sleep(backoff_delay)
                    continue
                
                # Reset failure counter on successful parse (even from text/plain)
                consecutive_failures = 0
                
                jobs = submission_data.get('jobs', [])
                user_images = submission_data.get('user_images', [])
                
                # Store discovered IDs for debugging
                if jobs:
                    discovered_job_id = jobs[0]
                    self.last_jobid = discovered_job_id
                    logger.info(f"DISCOVERED job_id: {discovered_job_id} - switching to job polling")
                    continue  # Next iteration will poll the job directly
                
                if user_images:
                    self.last_user_image_id = user_images[0]
                    logger.info(f"DISCOVERED user_image_id: {user_images[0]}")
                
                # Log submission processing status
                processing_started = submission_data.get('processing_started')
                processing_finished = submission_data.get('processing_finished')
                
                if processing_started and not processing_finished:
                    logger.info("SUBMISSION: Processing started, waiting for job creation...")
                elif processing_finished:
                    logger.info("SUBMISSION: Processing finished")
                else:
                    logger.info("SUBMISSION: Still being processed...")
                
                # Intelligent sleep interval
                sleep_time = self._calculate_submission_poll_interval(elapsed_time)
                logger.debug(f"Next submission check in {sleep_time}s...")
                time.sleep(sleep_time)
            
            logger.error(f"TIMEOUT: Astrometry solving timed out after {max_wait}s ({max_wait/60:.1f} minutes)")
            return None
            
        except Exception as e:
            logger.error(f"Fatal error while waiting for solution: {str(e)}")
            return None
    
    def _poll_job_until_complete(self, job_id, start_time, max_wait, last_status):
        """Poll job status until completion."""
        logger.info(f"JOB POLLING: Now monitoring job {job_id} directly")
        
        while time.time() - start_time < max_wait:
            elapsed_time = time.time() - start_time
            
            # Get job status with robust parsing
            job_data = self._get_job_status_robust(job_id)
            
            if job_data is None:
                logger.warning("Failed to get job status, retrying...")
                time.sleep(10)
                continue
            
            status = job_data.get('status', 'unknown')
            
            # Log status changes
            if status != last_status:
                logger.info(f"JOB STATUS CHANGE: {last_status} -> {status} (elapsed: {elapsed_time:.1f}s)")
                last_status = status
                
                # Update GUI if available
                if self.gui_ref:
                    self._update_gui_status(status, elapsed_time)
            
            if status == 'success':
                logger.info(f"JOB SUCCESS: Job {job_id} completed successfully!")
                return job_id
            elif status == 'failure':
                logger.error(f"JOB FAILED: Job {job_id} failed. Details: {job_data}")
                return None
            
            # Status logging
            logger.info(f"JOB STATUS: {status}")
            
            # Adaptive polling based on status
            sleep_time = self._calculate_job_poll_interval(elapsed_time, status)
            time.sleep(sleep_time)
        
        logger.error(f"JOB TIMEOUT: Job {job_id} timed out after {max_wait}s")
        return None
    
    def _update_gui_status(self, status, elapsed_time):
        """Update GUI with current processing status."""
        if hasattr(self.gui_ref, 'statusbar'):
            status_msg = f"Astrometry: {status} ({elapsed_time:.0f}s)"
            self.gui_ref.statusbar.config(text=status_msg)
            self.gui_ref.statusbar.update()
    
    def _calculate_submission_poll_interval(self, elapsed_time):
        """Calculate polling interval for submission status."""
        # More frequent initially, then back off
        if elapsed_time < 60:    return 8   # First minute: 8s
        elif elapsed_time < 180: return 15  # Next 2 minutes: 15s
        else:                    return 25  # After 3 minutes: 25s
    
    def _calculate_job_poll_interval(self, elapsed_time, status):
        """Calculate polling interval for job status."""
        if status == 'solving':           return 12  # Active solving: 12s
        elif status in ['queued']:        return 20  # Queued: 20s  
        elif status == 'processing':      return 15  # Processing: 15s
        elif elapsed_time < 300:          return 10  # First 5min: 10s
        else:                             return 18  # After 5min: 18s
    
    def _get_submission_status_robust(self, submission_id, max_retries=3):
        """Get submission status with robust parsing and retry logic."""
        status_url = f"{self.astrometry_base_url}submissions/{submission_id}"
        
        for attempt in range(max_retries):
            try:
                # Add explicit Accept header and session authentication
                headers = {"Accept": "application/json"}
                params = {}
                if self.current_session_key:
                    params['session'] = self.current_session_key
                else:
                    logger.warning(f"No session key available for submission status check (attempt {attempt + 1})")
                    
                response = self.session.get(status_url, headers=headers, params=params, timeout=45)
                
                # Handle 403 Forbidden (authentication issues)
                if response.status_code == 403:
                    logger.warning(f"HTTP 403 Forbidden for submission {submission_id} - attempting re-authentication")
                    if attempt < max_retries - 1 and self._reauthenticate():
                        logger.info("Re-authentication successful, retrying submission status check")
                        time.sleep(1)  # Brief delay before retry
                        continue
                    else:
                        logger.error("Re-authentication failed or no retries left")
                        response.raise_for_status()
                
                # Handle other HTTP errors with exponential backoff + jitter
                if response.status_code in [429, 500, 502, 503, 504]:
                    if attempt < max_retries - 1:
                        delay = (2 ** attempt) + random.uniform(0.1, 0.5)  # Add jitter
                        logger.warning(f"HTTP {response.status_code}, retrying in {delay:.1f}s (attempt {attempt + 1})")
                        time.sleep(delay)
                        continue
                
                response.raise_for_status()
                
                # Try response.json() first, fallback to robust parsing
                try:
                    return response.json()
                except (json.JSONDecodeError, ValueError):
                    logger.debug("Standard response.json() failed, trying robust parsing")
                
                # Use our robust parser
                parsed_data = self.parse_submission_body(response.text)
                if parsed_data:
                    if 'jobs' in parsed_data or 'user_images' in parsed_data:
                        logger.info(f"RECOVERY: Extracted data from non-JSON response: {parsed_data}")
                    return parsed_data
                    
                # If parsing failed, log and retry
                logger.warning(f"Failed to parse response on attempt {attempt + 1}")
                logger.debug(f"Response text (first 300 chars): {response.text[:300]}")
                
                if attempt < max_retries - 1:
                    delay = (2 ** attempt) + random.uniform(0.1, 0.3)
                    time.sleep(delay)
                    continue
                    
            except requests.exceptions.RequestException as e:
                logger.warning(f"Network error getting submission status (attempt {attempt + 1}): {str(e)}")
                if attempt < max_retries - 1:
                    delay = (2 ** attempt) + random.uniform(0.2, 0.8)
                    time.sleep(delay)
                    continue
        
        logger.error("Failed to get submission status after all retries")
        return None
    
    def _get_job_status_robust(self, job_id, max_retries=3):
        """Get job status with robust parsing and retry logic."""
        job_url = f"{self.astrometry_base_url}jobs/{job_id}"
        
        for attempt in range(max_retries):
            try:
                # Add explicit Accept header and session authentication
                headers = {"Accept": "application/json"}
                params = {}
                if self.current_session_key:
                    params['session'] = self.current_session_key
                else:
                    logger.warning(f"No session key available for job status check (attempt {attempt + 1})")
                    
                response = self.session.get(job_url, headers=headers, params=params, timeout=40)
                
                # Handle 403 Forbidden (authentication issues)
                if response.status_code == 403:
                    logger.warning(f"HTTP 403 Forbidden for job {job_id} - attempting re-authentication")
                    if attempt < max_retries - 1 and self._reauthenticate():
                        logger.info("Re-authentication successful, retrying job status check")
                        time.sleep(1)  # Brief delay before retry
                        continue
                    else:
                        logger.error("Re-authentication failed or no retries left")
                        response.raise_for_status()
                
                # Handle other HTTP errors with exponential backoff + jitter
                if response.status_code in [429, 500, 502, 503, 504]:
                    if attempt < max_retries - 1:
                        delay = (2 ** attempt) + random.uniform(0.1, 0.4)
                        logger.warning(f"HTTP {response.status_code} for job {job_id}, retrying in {delay:.1f}s")
                        time.sleep(delay)
                        continue
                
                response.raise_for_status()
                
                # Try response.json() first
                try:
                    return response.json()
                except (json.JSONDecodeError, ValueError):
                    logger.debug(f"Standard job response.json() failed for {job_id}, trying text parsing")
                
                # Fallback: try to extract status from text response
                text = response.text
                if text:
                    # Look for status in various formats
                    status_match = re.search(r'"status"\s*:\s*"([^"]+)"', text)
                    if status_match:
                        status = status_match.group(1)
                        logger.info(f"RECOVERY: Extracted job status from text: {status}")
                        return {"status": status}
                    
                    # Try other common patterns
                    if "success" in text.lower():
                        return {"status": "success"}
                    elif "failure" in text.lower() or "failed" in text.lower():
                        return {"status": "failure"}
                    elif "solving" in text.lower():
                        return {"status": "solving"}
                    elif "queued" in text.lower():
                        return {"status": "queued"}
                
                logger.warning(f"Could not parse job status from response on attempt {attempt + 1}")
                logger.debug(f"Job response text (first 300 chars): {text[:300]}")
                
                if attempt < max_retries - 1:
                    delay = (2 ** attempt) + random.uniform(0.1, 0.3)
                    time.sleep(delay)
                    continue
                    
            except requests.exceptions.RequestException as e:
                logger.warning(f"Network error getting job status (attempt {attempt + 1}): {str(e)}")
                if attempt < max_retries - 1:
                    delay = (2 ** attempt) + random.uniform(0.2, 0.6)
                    time.sleep(delay)
                    continue
        
        logger.error(f"Failed to get job status for {job_id} after all retries")
        return None
    
    def _calculate_poll_interval(self, elapsed_time, current_status):
        """Calculate intelligent polling interval based on elapsed time and status."""
        # Faster polling during critical phases, slower during long waits
        if current_status == 'solving':
            return 15  # Check every 15s while actively solving
        elif current_status in ['queued', 'processing']:
            return 10 if elapsed_time < 300 else 20  # 10s first 5min, then 20s
        else:
            # Progressive backoff: 10s → 20s → 30s
            if elapsed_time < 180:    # First 3 minutes
                return 10
            elif elapsed_time < 600:  # Next 7 minutes  
                return 20
            else:                     # After 10 minutes
                return 30
    
    def _get_calibration_data(self, job_id, max_retries=3):
        """Get calibration data from solved job with retry logic."""
        
        for attempt in range(max_retries):
            try:
                calib_url = f"{self.astrometry_base_url}jobs/{job_id}/calibration/"
                logger.debug(f"Getting calibration data from: {calib_url} (attempt {attempt + 1}/{max_retries})")
                
                response = self.session.get(calib_url, timeout=60)
                logger.debug(f"Calibration response: {response.status_code}")
                
                # Handle HTTP errors with retry
                if response.status_code in [429, 500, 502, 503, 504]:
                    if attempt < max_retries - 1:
                        delay = 2 ** attempt
                        logger.warning(f"HTTP {response.status_code} for calibration data, retrying in {delay}s")
                        time.sleep(delay)
                        continue
                
                response.raise_for_status()
                
                # Validate response content
                if not response.text.strip():
                    logger.warning("Empty calibration response")
                    if attempt < max_retries - 1:
                        time.sleep(2 ** attempt)
                        continue
                    else:
                        logger.error("Empty calibration response after all retries")
                        return {}
                
                try:
                    result = response.json()
                    logger.info(f"Retrieved calibration data successfully")
                    logger.debug(f"Calibration data: {result}")
                    return result
                except json.JSONDecodeError as e:
                    logger.warning(f"Failed to parse calibration JSON (attempt {attempt + 1}): {str(e)}")
                    if attempt < max_retries - 1:
                        time.sleep(2 ** attempt)
                        continue
                    else:
                        logger.error("Failed to parse calibration JSON after all retries")
                        return {}
                
            except requests.exceptions.RequestException as e:
                logger.warning(f"Network error getting calibration data (attempt {attempt + 1}): {str(e)}")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                else:
                    logger.error(f"Network error getting calibration data after {max_retries} attempts: {str(e)}")
                    return {}
                    
            except Exception as e:
                logger.warning(f"Error getting calibration data (attempt {attempt + 1}): {str(e)}")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                else:
                    logger.error(f"Error getting calibration data after {max_retries} attempts: {str(e)}")
                    return {}
        
        logger.error(f"Failed to get calibration data after {max_retries} attempts")
        return {}
    
    def get_last_astrometry_data(self):
        """
        Get the most recent successful astrometry data for external access.
        
        Returns:
            dict: Dictionary containing job_id, user_image_id, submission_id, etc.
                 Returns None if no successful data is available.
        """
        if not self.api_data_cache:
            logger.warning("No cached astrometry data available")
            return None
        
        # Get the most recent entry
        latest_job_id = max(self.api_data_cache.keys(), key=lambda k: self.api_data_cache[k]['timestamp'])
        latest_data = self.api_data_cache[latest_job_id].copy()
        
        logger.info(f"Retrieved cached astrometry data for job {latest_job_id}")
        return latest_data
    
    def get_astrometry_urls(self, job_id=None):
        """
        Generate URLs for downloading various astrometry products.
        
        Args:
            job_id: Specific job ID to use. If None, uses the most recent successful job.
            
        Returns:
            dict: Dictionary of URLs for different astrometry products
        """
        if job_id is None:
            if self.last_jobid:
                job_id = self.last_jobid
            else:
                logger.error("No job ID available for generating URLs")
                return {}
        
        base_url = f"{self.astrometry_base_url}jobs/{job_id}"
        urls = {
            'calibration': f"{base_url}/calibration/",
            'wcs_file': f"https://nova.astrometry.net/wcs_file/{job_id}",  # Correct WCS FITS endpoint
            'wcs_header': f"{base_url}/wcsheader/",  # Header fallback endpoint
            'annotated_image': f"{base_url}/annotated_display/",
            'new_fits': f"{base_url}/new_fits_file/",
            'rdls': f"{base_url}/rdls/",
            'corr': f"{base_url}/corr/",
            # Legacy API endpoint (deprecated)
            'wcs_api': f"{base_url}/wcs/"
        }
        
        logger.info(f"Generated astrometry URLs for job {job_id}")
        return urls
    
    def _get_wcs_solution(self, job_id, max_retries=5):
        """Get WCS solution using the correct FITS endpoint with comprehensive handling."""
        
        for attempt in range(max_retries):
            try:
                # Try primary WCS FITS file endpoint (correct Nova endpoint)
                wcs_result = self._fetch_wcs_fits_file(job_id, attempt, max_retries)
                if wcs_result is not None:
                    return wcs_result
                
                # If primary endpoint fails, try API header endpoint as fallback
                logger.info(f"WCS FITS file failed, trying header fallback (attempt {attempt + 1})")
                wcs_result = self._try_wcs_header_endpoint(job_id, attempt, max_retries)
                if wcs_result is not None:
                    return wcs_result
                
                # Exponential backoff with jitter for WCS file generation timing
                if attempt < max_retries - 1:
                    # Exponential backoff: 2s, 4s, 8s, 16s, 32s
                    delay_base = 2 ** (attempt + 1)
                    delay = delay_base + random.uniform(0.5, 2.0)
                    logger.warning(f"All WCS endpoints failed, retrying in {delay:.1f}s (attempt {attempt + 1}/{max_retries})")
                    time.sleep(delay)
                
            except Exception as e:
                logger.error(f"Unexpected error in WCS retrieval (attempt {attempt + 1}): {str(e)}")
                if attempt < max_retries - 1:
                    delay = 2 ** (attempt + 1)
                    time.sleep(delay)
                    continue
        
        logger.error(f"Failed to get WCS solution after {max_retries} attempts from all endpoints")
        return None
    
    def _fetch_wcs_fits_file(self, job_id, attempt, max_retries):
        """Fetch WCS from the correct Nova FITS endpoint with binary handling and gzip support."""
        try:
            # Use the correct Nova endpoint for WCS FITS files
            wcs_fits_url = f"https://nova.astrometry.net/wcs_file/{job_id}"
            logger.debug(f"Fetching WCS FITS from: {wcs_fits_url} (attempt {attempt + 1}/{max_retries})")
            
            # Prepare headers for binary FITS data
            headers = {
                'Accept': 'application/fits, application/octet-stream, */*',
                'Accept-Encoding': 'gzip, deflate',
                'User-Agent': 'AstrometryPythonClient/1.0'
            }
            
            # Add session authentication
            params = {}
            if self.current_session_key:
                params['session'] = self.current_session_key
            
            # Make request with redirect handling
            response = self.session.get(
                wcs_fits_url,
                headers=headers,
                params=params,
                timeout=120,
                allow_redirects=False  # Handle redirects manually
            )
            
            logger.debug(f"WCS FITS response: status={response.status_code}, url={response.url}")
            logger.debug(f"Response headers: {dict(response.headers)}")
            
            # Handle redirects
            if response.status_code in [302, 303, 307, 308]:
                redirect_url = response.headers.get('Location', '')
                logger.info(f"WCS FITS request redirected to: {redirect_url}")
                
                # Check for login redirects
                if 'login' in redirect_url.lower() or 'auth' in redirect_url.lower():
                    logger.warning("WCS FITS redirected to login page - attempting re-auth")
                    if attempt < max_retries - 1 and self._reauthenticate():
                        params['session'] = self.current_session_key
                
                # Follow redirect
                response = self.session.get(
                    redirect_url,
                    headers=headers,
                    params=params,
                    timeout=120
                )
            
            # Check for HTTP errors
            if response.status_code == 404:
                logger.warning(f"WCS FITS file not found for job {job_id} (may not be ready yet)")
                return None
            elif response.status_code == 429:
                logger.warning("Rate limited by Nova.astrometry.net")
                return None
            
            response.raise_for_status()
            
            # Process the binary response
            return self._process_fits_binary_response(response, job_id)
            
        except requests.exceptions.RequestException as e:
            logger.warning(f"Network error fetching WCS FITS file: {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Error fetching WCS FITS file: {str(e)}")
            return None
    
    def _process_fits_binary_response(self, response, job_id):
        """Process binary FITS response with gzip decompression and validation."""
        try:
            content = response.content
            content_length = len(content)
            content_type = response.headers.get('content-type', '').lower()
            content_encoding = response.headers.get('content-encoding', '').lower()
            
            logger.debug(f"FITS content: length={content_length}, type={content_type}, encoding={content_encoding}")
            
            # Check for minimum file size
            if content_length < 100:
                logger.warning(f"WCS FITS content too small ({content_length} bytes)")
                return None
            
            # Handle gzip compression
            if content_encoding == 'gzip' or self._is_gzipped(content):
                logger.debug("Decompressing gzipped WCS FITS content")
                try:
                    content = gzip.decompress(content)
                    logger.debug(f"Decompressed size: {len(content)} bytes")
                except Exception as e:
                    logger.warning(f"Failed to decompress gzipped content: {str(e)}")
                    return None
            
            # Check for HTML error pages
            if content_type and 'text/html' in content_type:
                return self._handle_fits_html_response(content, job_id)
            
            # Parse the FITS file
            return self._parse_fits_content(content)
            
        except Exception as e:
            logger.error(f"Error processing FITS binary response: {str(e)}")
            return None
    
    def _is_gzipped(self, content):
        """Check if content is gzip compressed using magic bytes."""
        return len(content) >= 2 and content[:2] == b'\x1f\x8b'
    
    def _handle_fits_html_response(self, content, job_id):
        """Handle HTML responses from FITS endpoint."""
        try:
            html_text = content.decode('utf-8', errors='ignore')[:200].lower()
            logger.warning(f"WCS FITS endpoint returned HTML content")
            logger.debug(f"HTML snippet: {html_text[:150]}")
            
            if 'login' in html_text or 'sign in' in html_text:
                logger.error("WCS FITS endpoint requires authentication")
            elif 'not found' in html_text or '404' in html_text:
                logger.warning(f"WCS FITS file not found for job {job_id}")
            elif 'rate limit' in html_text:
                logger.error("Rate limited by Nova.astrometry.net")
            else:
                logger.warning("Unknown HTML response from FITS endpoint")
                
            return None
            
        except Exception as e:
            logger.debug(f"Error analyzing HTML response: {str(e)}")
            return None
    
    def _parse_fits_content(self, content):
        """Safely parse FITS content with proper validation."""
        try:
            # Use BytesIO for in-memory FITS parsing
            with BytesIO(content) as fits_buffer:
                with fits.open(fits_buffer, ignore_missing_end=True, ignore_missing_simple=True) as hdulist:
                    if len(hdulist) == 0:
                        logger.warning("FITS file contains no HDUs")
                        return None
                    
                    header = hdulist[0].header
                    logger.debug(f"FITS header contains {len(header)} keys")
                    
                    # Validate essential WCS keys
                    required_keys = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2']
                    missing_keys = [key for key in required_keys if key not in header]
                    
                    if missing_keys:
                        logger.warning(f"FITS missing required WCS keys: {missing_keys}")
                        return None
                    
                    # Check for CD matrix or CDELT scaling
                    has_cd_matrix = any(key in header for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'])
                    has_cdelt = any(key in header for key in ['CDELT1', 'CDELT2'])
                    
                    if not (has_cd_matrix or has_cdelt):
                        logger.warning("FITS missing CD matrix or CDELT scaling information")
                        return None
                    
                    # Log key WCS parameters for debugging
                    wcs_keys = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1', 'CDELT2', 'CTYPE1', 'CTYPE2']
                    logger.debug("WCS parameters:")
                    for key in wcs_keys:
                        if key in header:
                            logger.debug(f"  {key}: {header[key]}")
                    
                    # Create WCS object
                    wcs = WCS(header)
                    
                    # Validate WCS object
                    if not wcs.has_celestial:
                        logger.warning("WCS object does not have celestial coordinates")
                        return None
                    
                    logger.info("Successfully parsed WCS from FITS file")
                    return wcs
                    
        except Exception as e:
            logger.warning(f"Failed to parse FITS content: {str(e)}")
            return None
    
    
    def _try_wcs_header_endpoint(self, job_id, attempt, max_retries):
        """Try fallback WCS header endpoint that returns plain text."""
        try:
            wcs_header_url = f"{self.astrometry_base_url}jobs/{job_id}/wcsheader/"
            logger.debug(f"Trying WCS header fallback: {wcs_header_url}")
            
            headers = {'Accept': 'text/plain, application/json;q=0.8'}
            params = {}
            if self.current_session_key:
                params['session'] = self.current_session_key
            else:
                logger.warning("No session key available for WCS header request")
            
            response = self.session.get(
                wcs_header_url,
                headers=headers,
                params=params,
                timeout=90,
                allow_redirects=False
            )
            
            # Handle redirects
            if response.status_code in [302, 303]:
                redirect_url = response.headers.get('Location', '')
                logger.info(f"WCS header request redirected, following: {redirect_url}")
                response = self.session.get(redirect_url, headers=headers, params=params, timeout=90)
            
            response.raise_for_status()
            
            # Parse as text header
            if response.text.strip():
                wcs = self._parse_text_header(response.text)
                if wcs:
                    logger.info("Successfully obtained WCS from header fallback endpoint")
                    return wcs
            
            return None
            
        except Exception as e:
            logger.warning(f"Error accessing WCS header fallback endpoint: {str(e)}")
            return None
    
    def _handle_html_response(self, response, job_id, attempt=0, max_retries=3):
        """Analyze HTML response for authentication or rate limiting issues."""
        html_content = response.text[:500].lower()  # First 500 chars for analysis
        
        logger.warning(f"WCS endpoint returned HTML content")
        logger.debug(f"Response status: {response.status_code}")
        logger.debug(f"Response URL: {response.url}")
        logger.debug(f"HTML snippet (first 150 chars): {html_content[:150]}")
        
        # Detect common issues and determine action
        if 'login' in html_content or 'sign in' in html_content or 'auth' in html_content:
            logger.warning("WCS request appears to need authentication - will attempt re-auth")
            if attempt < max_retries - 1:  # Only try re-auth if we have retries left
                return 'retry_with_new_auth'
            else:
                logger.error("Authentication failed after retries")
        elif 'rate limit' in html_content or 'too many requests' in html_content:
            logger.error("Rate limited by astrometry.net server")
        elif 'not found' in html_content or '404' in html_content:
            logger.error(f"WCS data not found for job {job_id}")
        elif 'error' in html_content or 'exception' in html_content:
            logger.warning("Server error in HTML response")
        else:
            logger.warning("Unknown HTML response - possibly transient server issue")
        
        return None
    
    def _parse_binary_fits(self, content):
        """Legacy method - delegates to the new comprehensive FITS parser."""
        return self._parse_fits_content(content)
    
    def _parse_text_header(self, text_content):
        """Try to parse WCS from text header format with forgiving END handling."""
        try:
            if not text_content or not text_content.strip():
                return None
            
            text_content = text_content.strip()
            
            # Check if this looks like FITS header text
            if 'CRVAL1' not in text_content and 'RA' not in text_content:
                logger.debug("Text content does not appear to contain WCS information")
                return None
            
            # First attempt: try parsing as-is
            try:
                header = fits.Header.fromstring(text_content, sep='\n')
            except Exception as e:
                # If END is missing, append it and retry
                if 'END' not in text_content or not text_content.endswith('END'):
                    logger.debug("Adding END to FITS header text")
                    text_content_with_end = text_content + '\nEND\n'
                    try:
                        header = fits.Header.fromstring(text_content_with_end, sep='\n')
                    except Exception as e2:
                        logger.debug(f"Failed to parse text header even with END: {str(e2)}")
                        return None
                else:
                    logger.debug(f"Failed to parse text header: {str(e)}")
                    return None
            
            # Verify WCS keys
            required_keys = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2']
            missing_keys = [key for key in required_keys if key not in header]
            
            if missing_keys:
                logger.debug(f"Text header missing WCS keys: {missing_keys}")
                return None
            
            wcs = WCS(header)
            logger.info("Successfully parsed WCS from text header")
            return wcs
            
        except Exception as e:
            logger.debug(f"Failed to parse text header: {str(e)}")
            return None
    
    def _cache_early_data(self, job_id, calibration_data):
        """Cache API data early, even if WCS retrieval fails."""
        try:
            self.api_data_cache[job_id] = {
                'calibration_data': calibration_data,
                'job_id': job_id,
                'submission_id': self.last_subid,
                'user_image_id': self.last_user_image_id,
                'timestamp': time.time(),
                'wcs_solution': None  # Will be updated if WCS is successful
            }
            logger.debug(f"Cached early astrometry data for job {job_id}")
        except Exception as e:
            logger.warning(f"Failed to cache early data: {str(e)}")
    
    def _reauthenticate(self):
        """Re-authenticate with astrometry.net API."""
        try:
            logger.info("Attempting to re-authenticate with astrometry.net")
            old_session_key = self.current_session_key
            self.current_session_key = None  # Clear old session
            
            new_session_key = self._astrometry_login()
            if new_session_key:
                logger.info("Re-authentication successful")
                return True
            else:
                logger.error("Re-authentication failed")
                self.current_session_key = old_session_key  # Restore old session as fallback
                return False
        except Exception as e:
            logger.error(f"Error during re-authentication: {str(e)}")
            return False
    
    def _solve_local_astrometry(self, filepath, image):
        """
        Solve astrometry using local solve-field command.
        
        Args:
            filepath: Path to FITS file
            image: Image data array
            
        Returns:
            tuple: (WCS solution, astrometry_data dict) or (None, None)
        """
        try:
            import subprocess
            import shutil
            from astropy.wcs import WCS
            from astropy.io import fits
            
            logger.info("="*60)
            logger.info("STARTING LOCAL ASTROMETRY WITH SOLVE-FIELD")
            logger.info("="*60)
            logger.info(f"Input file: {filepath}")
            
            # Update progress
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar2'):
                self.gui_ref.progressbar2['value'] = 80
                self.gui_ref.progressbar2.update()
            
            # Check WSL environment is setup and perform preflight checks
            if not self.wsl_astrometry_path:
                logger.error("WSL environment not properly initialized")
                return None, None
            
            # Perform comprehensive preflight checks
            if not self._preflight_check_wsl():
                logger.error("WSL preflight checks failed - skipping to Nova API")
                return None, None
            
            filename = os.path.basename(filepath)
            wsl_fits_path = f"{self.wsl_astrometry_path}/{filename}"
            
            # Convert Windows path to WSL format for copying
            win_path_for_wsl = filepath.replace("\\", "/").replace("C:", "/mnt/c")
            
            # Copy file to WSL work directory (should be writable by current user)
            copy_cmd = ["wsl", "-d", "Ubuntu", "cp", win_path_for_wsl, wsl_fits_path]
            logger.info(f"Copying FITS file inside WSL: {win_path_for_wsl} -> {wsl_fits_path}")
            copy_result = subprocess.run(copy_cmd, capture_output=True, text=True)
            
            if copy_result.returncode != 0:
                logger.error(f"Failed to copy file to WSL: {copy_result.stderr}")
                return None, None
            
            # Set proper permissions on the copied file
            chmod_cmd = ["wsl", "-d", "Ubuntu", "chmod", "664", wsl_fits_path]
            chmod_result = subprocess.run(chmod_cmd, capture_output=True, text=True)
            if chmod_result.returncode != 0:
                logger.warning(f"Failed to set file permissions: {chmod_result.stderr}")
            
            # Validate required fields before constructing command
            distro = "Ubuntu"
            solve_binary = self.wsl_solve_field_path or "/usr/bin/solve-field"
            workdir = self.wsl_astrometry_path
            input_path = wsl_fits_path
            config_path = self.wsl_config_path or "/home/burhan/astrometry-4211-4216.cfg"
            
            # Using proven working solve-field parameters (hardcoded for reliability)
            
            # Log the WSL config path for debugging
            logger.info("WSL config path (linux): %r", self.wsl_config_path)
            
            # Validate all required fields
            validation_errors = []
            if not distro:
                validation_errors.append("WSL distro not specified")
            if not solve_binary:
                validation_errors.append("solve-field binary path not found")
            if not workdir:
                validation_errors.append("WSL work directory not set")
            if not input_path:
                validation_errors.append("Input FITS path not set")
            if not config_path:
                validation_errors.append("Astrometry config path not set")
            
            if validation_errors:
                error_msg = "solve-field validation failed: " + ", ".join(validation_errors)
                logger.error(error_msg)
                return None, None
            
            # Log validated configuration
            logger.info("solve-field validation passed:")
            logger.info("  Distro: %s", distro)
            logger.info("  Binary: %s", solve_binary)
            logger.info("  Workdir: %s", workdir)
            logger.info("  Input: %s", input_path)
            logger.info("  Config: %s", config_path)
            
            # Construct solve-field command with validated parameters (using proven working settings)
            solve_cmd = [
                "wsl", "-d", distro,            # e.g., "Ubuntu"
                solve_binary,                   # e.g., "/usr/bin/solve-field"
                input_path,                     # e.g., "/home/burhan/astrometry-work/v2.fits"

                "--config", config_path,        # e.g., "/home/burhan/astrometry-4211-4216.cfg"
                "--dir", workdir,               # e.g., "/home/burhan/astrometry-work"

                "--overwrite", "--no-plots", "-v",

                # Matches the working solve:
                "--scale-units", "degwidth", "--scale-low", "4.7", "--scale-high", "5.3",
                "--downsample", "4",
                "--objs", "120",
                "--depth", "1-120",
                "--uniformize", "10",
                "--resort",
                "--nsigma", "12",
                "--code-tolerance", "0.04",
                "--odds-to-solve", "1e6",
                "--parity", "pos",
                "--crpix-center",
                "--cpulimit", "300",
            ]
            
            logger.info(f"Running solve-field command: {' '.join(solve_cmd)}")
            
            # Update status
            if self.gui_ref and hasattr(self.gui_ref, 'statusbar'):
                self.gui_ref.statusbar.config(text="Local Astrometry: Running solve-field...")
                self.gui_ref.statusbar.update()
            
            # Run solve-field (no Python timeout - solve-field manages its own with --cpulimit)
            result = subprocess.run(
                solve_cmd,
                capture_output=True,
                text=True
            )
            
            if result.returncode == 0:
                logger.info("solve-field completed successfully")
                
                # Look for the .new file (solved FITS with WCS) in WSL
                base_name = os.path.splitext(filename)[0]
                wcs_file_wsl = f"{self.wsl_astrometry_path}/{base_name}.new"
                
                # Check if the .new file exists in WSL
                check_cmd = ["wsl", "-d", "Ubuntu", "test", "-f", wcs_file_wsl]
                check_result = subprocess.run(check_cmd, capture_output=True)
                
                if check_result.returncode == 0:
                    logger.info(f"Found WCS solution file in WSL: {wcs_file_wsl}")
                    
                    # Copy the solved file back to Windows temp directory for processing
                    temp_dir = os.path.dirname(filepath)
                    temp_wcs_file = os.path.join(temp_dir, f"{base_name}.new")
                    
                    # Convert Windows path for WSL
                    win_temp_path = temp_wcs_file.replace("\\", "/").replace("C:", "/mnt/c")
                    copy_back_cmd = ["wsl", "-d", "Ubuntu", "cp", wcs_file_wsl, win_temp_path]
                    subprocess.run(copy_back_cmd, check=True)
                    
                    # Read WCS from solved file
                    with fits.open(temp_wcs_file) as hdul:
                        header = hdul[0].header
                        wcs = WCS(header)
                        
                        # Extract astrometry data from header
                        astrometry_data = self._extract_local_astrometry_data(header)
                        
                        logger.info("Local astrometry solution successful!")
                        logger.info(f"Center RA: {astrometry_data.get('center_ra', 'N/A')}")
                        logger.info(f"Center Dec: {astrometry_data.get('center_dec', 'N/A')}")
                        logger.info(f"Pixel scale: {astrometry_data.get('pixscale', 'N/A')} arcsec/pixel")
                        
                        # Clean up temp file
                        try:
                            os.remove(temp_wcs_file)
                        except:
                            pass
                        
                        # Update progress
                        if self.gui_ref and hasattr(self.gui_ref, 'progressbar2'):
                            self.gui_ref.progressbar2['value'] = 95
                            self.gui_ref.progressbar2.update()
                        
                        # Update status
                        if self.gui_ref and hasattr(self.gui_ref, 'statusbar'):
                            self.gui_ref.statusbar.config(text="Local Astrometry: Complete")
                            self.gui_ref.statusbar.update()
                        
                        return wcs, astrometry_data
                else:
                    logger.error("solve-field completed but no .new file found")
            else:
                logger.error(f"solve-field failed with return code: {result.returncode}")
                logger.error(f"stdout: {result.stdout}")
                logger.error(f"stderr: {result.stderr}")
            
            return None, None
            
        except Exception as e:
            logger.error(f"Error in local astrometry: {str(e)}")
            return None, None
    
    def _extract_local_astrometry_data(self, header):
        """
        Extract astrometry data from a FITS header solved by solve-field.
        
        Args:
            header: FITS header from solved file
            
        Returns:
            dict: Normalized astrometry data
        """
        try:
            astrometry_data = {}
            
            # Center coordinates
            if 'CRVAL1' in header and 'CRVAL2' in header:
                astrometry_data['center_ra'] = float(header['CRVAL1'])
                astrometry_data['center_dec'] = float(header['CRVAL2'])
            
            # Pixel scale (convert from degrees to arcseconds)
            if 'CD1_1' in header and 'CD2_2' in header:
                cd1_1 = abs(float(header['CD1_1']))
                cd2_2 = abs(float(header['CD2_2']))
                pixscale = (cd1_1 + cd2_2) / 2.0 * 3600.0  # degrees to arcsec
                astrometry_data['pixscale'] = pixscale
            
            # Field dimensions
            if 'NAXIS1' in header and 'NAXIS2' in header and 'pixscale' in astrometry_data:
                width_pixels = int(header['NAXIS1'])
                height_pixels = int(header['NAXIS2'])
                pixscale = astrometry_data['pixscale']
                
                field_width_arcmin = (width_pixels * pixscale) / 60.0
                field_height_arcmin = (height_pixels * pixscale) / 60.0
                
                astrometry_data['field_width_arcmin'] = field_width_arcmin
                astrometry_data['field_height_arcmin'] = field_height_arcmin
            
            # Orientation and parity (if available)
            if 'CD1_2' in header and 'CD2_1' in header and 'CD1_1' in header and 'CD2_2' in header:
                cd1_1 = abs(float(header['CD1_1']))
                cd1_2 = float(header['CD1_2'])
                cd2_1 = float(header['CD2_1'])
                cd2_2 = abs(float(header['CD2_2']))
                orientation = np.degrees(np.arctan2(cd1_2, cd2_1))
                astrometry_data['orientation'] = orientation
                
                # Parity (sign of determinant)
                det = cd1_1 * cd2_2 - cd1_2 * cd2_1
                astrometry_data['parity'] = 1 if det > 0 else -1
            
            # Add metadata
            astrometry_data['has_wcs'] = True
            astrometry_data['has_astrometry'] = True
            astrometry_data['solver'] = 'local_solve_field'
            
            logger.info("Successfully extracted local astrometry data")
            return astrometry_data
            
        except Exception as e:
            logger.error(f"Error extracting local astrometry data: {str(e)}")
            return {'has_wcs': False, 'has_astrometry': False, 'solver': 'local_solve_field'}
    
    def set_astrometry_mode(self, use_local=True):
        """
        Set the astrometry solving mode.
        
        Args:
            use_local (bool): If True, use ONLY local solve-field (no fallback).
                            If False, use fallback strategy (local first, then Nova API).
        """
        self.use_local_astrometry = use_local
        mode = "local solve-field only (no fallback)" if use_local else "fallback strategy (local first, then Nova API)"
        logger.info(f"Astrometry mode set to: {mode}")
    
    def get_astrometry_mode(self):
        """
        Get the current astrometry solving mode.
        
        Returns:
            str: Current astrometry mode description
        """
        return "local solve-field only (no fallback)" if self.use_local_astrometry else "fallback strategy (local first, then Nova API)"
    
    def update_wsl_settings(self, new_settings):
        """
        Update WSL settings and save to config file.
        
        Args:
            new_settings (dict): Dictionary containing new settings to update
            
        Returns:
            bool: True if settings were successfully updated and saved
        """
        try:
            # Update settings
            self.wsl_settings.update(new_settings)
            
            # Expand paths if needed
            if 'wsl_workdir' in new_settings:
                expanded_workdir = self.settings_manager.expand_wsl_path(new_settings['wsl_workdir'])
                if expanded_workdir:
                    self.wsl_astrometry_path = expanded_workdir
                    self.wsl_settings['wsl_workdir'] = expanded_workdir
            
            if 'wsl_config_file' in new_settings:
                expanded_config = self.settings_manager.expand_wsl_path(new_settings['wsl_config_file'])
                if expanded_config:
                    self.wsl_config_path = expanded_config
                    self.wsl_settings['wsl_config_file'] = expanded_config
            
            # Save settings
            success = self.settings_manager.save_settings(self.wsl_settings)
            if success:
                logger.info("WSL settings updated successfully")
            
            return success
            
        except Exception as e:
            logger.error(f"Failed to update WSL settings: {e}")
            return False
    
    def get_wsl_settings(self):
        """
        Get current WSL settings.
        
        Returns:
            dict: Current WSL settings
        """
        return self.wsl_settings.copy()
    
    
    def _normalize_calibration_data(self, calibration_data):
        """Normalize calibration data keys for CSV saving."""
        if not calibration_data:
            return {}
        
        try:
            normalized = {}
            
            # Map standard calibration keys
            if 'ra' in calibration_data:
                normalized['center_ra'] = float(calibration_data['ra'])
            if 'dec' in calibration_data:
                normalized['center_dec'] = float(calibration_data['dec'])
            
            # Convert arcsec to arcmin for field dimensions
            if 'width_arcsec' in calibration_data:
                normalized['field_width_arcmin'] = float(calibration_data['width_arcsec']) / 60.0
            if 'height_arcsec' in calibration_data:
                normalized['field_height_arcmin'] = float(calibration_data['height_arcsec']) / 60.0
            
            # Direct mappings
            if 'pixscale' in calibration_data:
                normalized['pixscale'] = float(calibration_data['pixscale'])
            if 'orientation' in calibration_data:
                normalized['orientation'] = float(calibration_data['orientation'])
            if 'parity' in calibration_data:
                normalized['parity'] = int(calibration_data['parity'])
            
            logger.debug(f"Normalized calibration data: {normalized}")
            return normalized
            
        except Exception as e:
            logger.warning(f"Error normalizing calibration data: {str(e)}")
            return {}

    def _pixel_to_sky(self, sources, wcs):
        """Convert pixel coordinates to sky coordinates."""
        try:
            if wcs is None:
                return None
                
            x_pix = sources['xcentroid']
            y_pix = sources['ycentroid']
            
            # Convert to sky coordinates
            sky_coords = wcs.pixel_to_world(x_pix, y_pix)
            
            logger.info(f"Converted {len(sources)} sources to sky coordinates")
            return sky_coords
            
        except Exception as e:
            logger.error(f"Error converting to sky coordinates: {str(e)}")
            return None
    
    def _cross_match_catalogs(self, sky_coords):
        """Cross-match sources with Gaia and SIMBAD catalogs using direct TAP queries."""
        try:
            if sky_coords is None:
                return None
                
            logger.info(f"Starting catalog cross-matching for {len(sky_coords)} sources")
            
            # Calculate field boundaries with some padding
            ras = [coord.ra.degree for coord in sky_coords]
            decs = [coord.dec.degree for coord in sky_coords]
            
            ra_min, ra_max = min(ras) - 0.01, max(ras) + 0.01
            dec_min, dec_max = min(decs) - 0.01, max(decs) + 0.01
            
            # Calculate field center and radius for circular queries
            center_ra = (ra_min + ra_max) / 2
            center_dec = (dec_min + dec_max) / 2
            
            # Calculate radius to encompass the field (with some padding)
            max_separation = max(
                np.sqrt((ra - center_ra)**2 + (dec - center_dec)**2) 
                for ra, dec in zip(ras, decs)
            )
            search_radius = max_separation + 0.02  # Add padding
            
            logger.info(f"Field center: RA={center_ra:.4f}, DEC={center_dec:.4f}")
            logger.info(f"Search radius: {search_radius:.4f} degrees ({search_radius*3600:.1f} arcsec)")
            
            # Update progress
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar2'):
                self.gui_ref.progressbar2['value'] = 85
                self.gui_ref.progressbar2.update()
            
            # STEP 1: Optimized Gaia query (single async CIRCLE with RUWE/magnitude filtering)
            gaia_matches_dict = self._query_gaia_optimized(sky_coords, field_radius_deg=search_radius)
            logger.info(f"Optimized Gaia matching: {len(gaia_matches_dict)}/{len(sky_coords)} sources matched")
            
            # Convert optimized matches to the expected format
            gaia_matches = []
            for i in range(len(sky_coords)):
                if i in gaia_matches_dict:
                    gaia_matches.append(gaia_matches_dict[i])
                else:
                    # Empty match for unmatched sources
                    gaia_matches.append({
                        'gaia_source_id': '',
                        'gaia_ra': np.nan,
                        'gaia_dec': np.nan,
                        'phot_g_mean_mag': np.nan,
                        'phot_bp_mean_mag': np.nan,
                        'phot_rp_mean_mag': np.nan,
                        'bp_rp': np.nan,
                        'parallax': np.nan,
                        'pmra': np.nan,
                        'pmdec': np.nan,
                        'ruwe': np.nan,
                        'match_distance_arcsec': np.nan
                    })
            
            # Update progress
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar2'):
                self.gui_ref.progressbar2['value'] = 90
                self.gui_ref.progressbar2.update()
            
            # STEP 2: Query SIMBAD using efficient X-Match approach
            logger.info(f"Performing SIMBAD X-Match for {len(sky_coords)} sources...")
            simbad_matches = self._simbad_xmatch(sky_coords)
            
            # Update progress
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar2'):
                self.gui_ref.progressbar2['value'] = 95
                self.gui_ref.progressbar2.update()
            
            logger.info("Catalog cross-matching completed successfully")
            return {
                'gaia': gaia_matches,
                'simbad': simbad_matches  # Direct TAP queries for each source
            }
            
        except Exception as e:
            logger.error(f"Error in catalog cross-matching: {str(e)}")
            return None
    
    def _reinitialize_simbad_tap(self):
        """Reinitialize SIMBAD TAP connection after connection failures."""
        try:
            logger.info("Reinitializing SIMBAD TAP connection...")
            self.simbad_tap = TapPlus(url="https://simbad.cds.unistra.fr/simbad/sim-tap")
            self.simbad_connection_failures = 0
            logger.info("SIMBAD TAP connection reinitialized successfully")
            return True
        except Exception as e:
            logger.warning(f"SIMBAD TAP reinitialization failed: {e}")
            self.simbad_tap = None
            return False
    
    def _query_simbad_tap(self, ra, dec, radius_arcsec=5.0, max_retries=4):
        """
        Query SIMBAD using TAP service for a single coordinate with enhanced connection error prevention.
        
        Args:
            ra: Right ascension in degrees
            dec: Declination in degrees
            radius_arcsec: Search radius in arcseconds
            max_retries: Maximum number of retry attempts (increased default)
            
        Returns:
            Dictionary with SIMBAD data or empty values if no match
        """
        # Default empty result
        empty_result = {
            'main_id': '',
            'otype': '',
            'sp_type': '',
            'rv_value': np.nan,
            'distance_result': np.nan,
            'match_distance_arcsec': np.nan,
            'ra': np.nan,
            'dec': np.nan
        }
        
        if self.simbad_tap is None:
            logger.warning("SIMBAD TAP service not available")
            return empty_result
        
        # Convert radius to degrees
        radius_deg = radius_arcsec / 3600.0
        
        # ADQL query for SIMBAD using confirmed available columns
        adql_query = f"""
        SELECT TOP 1 
            main_id, 
            otype,
            sp_type,
            ra,
            dec,
            DISTANCE(
                POINT('ICRS', {ra}, {dec}),
                POINT('ICRS', ra, dec)
            ) AS match_distance_deg
        FROM basic 
        WHERE 1=CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra}, {dec}, {radius_deg})
        )
        ORDER BY match_distance_deg ASC
        """
        
        # Retry logic with exponential backoff and preventive measures
        for attempt in range(max_retries + 1):
            try:
                # Preventive measure: Add conservative delay before retry attempts
                if attempt > 0:
                    base_delay = 3.0  # Start with longer base delay for retries
                    exponential_factor = 2 ** (attempt - 1)
                    jitter = random.uniform(-0.5, 0.5)  # ±0.5s jitter
                    delay = (base_delay * exponential_factor) + jitter  # 3s, 6s, 12s, 24s (with jitter)
                    logger.info(f"Pre-query delay with jitter: {delay:.1f}s (attempt {attempt + 1})")
                    time.sleep(delay)
                
                # Additional preventive measure: Check connection health before important queries
                if attempt > 0:
                    try:
                        # Test with a minimal query to check connection health
                        test_query = "SELECT TOP 1 main_id FROM basic WHERE main_id LIKE '* alf%'"
                        test_job = self.simbad_tap.launch_job(test_query)
                        test_results = test_job.get_results()
                        logger.debug("Connection health check passed")
                    except Exception as health_e:
                        logger.warning(f"Connection health check failed: {health_e}")
                        # Try to reinitialize connection if health check fails
                        if not self._reinitialize_simbad_tap():
                            raise Exception("Connection health check failed and reinitialization unsuccessful")
                
                # Execute the actual query
                logger.debug(f"Executing SIMBAD query for RA={ra:.6f}, Dec={dec:.6f} (attempt {attempt + 1})")
                job = self.simbad_tap.launch_job(adql_query)
                results = job.get_results()
                
                if results is not None and len(results) > 0:
                    row = results[0]
                    
                    result = {
                        'main_id': str(row['main_id']) if row['main_id'] is not np.ma.masked else '',
                        'otype': str(row['otype']) if 'otype' in row.colnames and row['otype'] is not np.ma.masked else '',
                        'sp_type': str(row['sp_type']) if 'sp_type' in row.colnames and row['sp_type'] is not np.ma.masked else '',
                        'rv_value': np.nan,  # Not available in basic table
                        'distance_result': np.nan,  # Not available in basic table
                        'match_distance_arcsec': float(row['match_distance_deg']) if 'match_distance_deg' in row.colnames and row['match_distance_deg'] is not np.ma.masked else np.nan,
                        'ra': float(row['ra']) if 'ra' in row.colnames and row['ra'] is not np.ma.masked else np.nan,
                        'dec': float(row['dec']) if 'dec' in row.colnames and row['dec'] is not np.ma.masked else np.nan
                    }
                    
                    # Convert match distance from degrees to arcseconds
                    if not np.isnan(result['match_distance_arcsec']):
                        result['match_distance_arcsec'] *= 3600.0
                    
                    return result
                else:
                    return empty_result
                    
            except Exception as e:
                # Check for specific connection errors that should trigger retry
                error_msg = str(e).lower()
                is_connection_error = any(keyword in error_msg for keyword in [
                    'connection', 'timeout', 'forcibly closed', 'network', 
                    'temporarily unavailable', 'service unavailable', 'http 503',
                    'http 502', 'http 500'
                ])
                
                if attempt < max_retries and is_connection_error:
                    self.simbad_connection_failures += 1
                    logger.warning(f"SIMBAD TAP query failed (attempt {attempt + 1}/{max_retries + 1}) for RA={ra:.6f}, Dec={dec:.6f}: {e}")
                    
                    # Progressive connection recovery strategy
                    if self.simbad_connection_failures >= 3:
                        logger.info("Multiple connection failures detected, attempting connection recovery...")
                        
                        # Try to reinitialize the connection
                        if self._reinitialize_simbad_tap():
                            logger.info("SIMBAD connection reinitialized successfully")
                        else:
                            logger.warning("SIMBAD connection reinitialization failed")
                            
                        # Add extra delay after reinitialization
                        recovery_delay = 3.0
                        logger.info(f"Connection recovery delay: {recovery_delay}s")
                        time.sleep(recovery_delay)
                    
                    # The delay is now handled at the start of the next iteration
                    continue
                else:
                    # Final failure or non-retryable error
                    if is_connection_error:
                        self.simbad_connection_failures += 1
                    logger.warning(f"SIMBAD TAP query failed for RA={ra}, Dec={dec}: {e}")
                    return empty_result
        
        # Should not reach here, but just in case
        return empty_result
    
    def _query_simbad_batch(self, sky_coords, radius_arcsec=5.0):
        """
        Query SIMBAD for all sources using a more accurate batch approach that mimics individual cone searches.
        
        Args:
            sky_coords: List of SkyCoord objects for all sources
            radius_arcsec: Search radius in arcseconds
            
        Returns:
            List of dictionaries with SIMBAD data for each source (in same order)
        """
        # Create empty results array matching input size
        empty_result = {
            'main_id': '',
            'otype': '',
            'sp_type': '',
            'rv_value': np.nan,
            'distance_result': np.nan,
            'match_distance_arcsec': np.nan,
            'ra': np.nan,
            'dec': np.nan
        }
        
        simbad_matches = [empty_result.copy() for _ in sky_coords]
        
        if self.simbad_tap is None:
            logger.warning("SIMBAD TAP service not available")
            return simbad_matches
        
        if len(sky_coords) == 0:
            return simbad_matches
        
        try:
            logger.info(f"Executing improved batch SIMBAD query for {len(sky_coords)} sources...")
            
            # Convert radius to degrees
            radius_deg = radius_arcsec / 3600.0
            
            # Extract coordinates
            ras = [coord.ra.degree for coord in sky_coords]
            decs = [coord.dec.degree for coord in sky_coords]
            
            # Since UNION is not supported by SIMBAD ADQL, use field-based approach for all batch sizes
            # but optimize the field size calculation for better accuracy
            
            logger.info("Using optimized field approach for batch query")
            
            # Calculate field boundaries more carefully
            if len(sky_coords) <= 20:  # For small batches, use tight field with small padding
                padding_factor = 1.5
            elif len(sky_coords) <= 50:  # Medium batches get moderate padding
                padding_factor = 2.0
            else:  # Large batches need generous padding
                padding_factor = 3.0
            
            # Calculate bounding box with optimized padding
            ra_min, ra_max = min(ras) - radius_deg*padding_factor, max(ras) + radius_deg*padding_factor
            dec_min, dec_max = min(decs) - radius_deg*padding_factor, max(decs) + radius_deg*padding_factor
            
            # Handle RA wrap-around at 0/360 degrees
            if ra_max - ra_min > 180:  # Likely crossing 0/360 boundary
                logger.warning("Field crosses RA 0/360 boundary - using special handling")
                if ra_min < 180:
                    ra_condition = f"(ra >= {ra_min} OR ra <= {ra_max - 360})"
                else:
                    ra_condition = f"(ra >= {ra_min - 360} OR ra <= {ra_max})"
            else:
                ra_condition = f"ra BETWEEN {ra_min} AND {ra_max}"
            
            adql_query = f"""
            SELECT 
                main_id, 
                otype,
                sp_type,
                ra,
                dec
            FROM basic 
            WHERE {ra_condition}
            AND dec BETWEEN {dec_min} AND {dec_max}
            """
            
            logger.info("Executing SIMBAD batch query...")
            
            # Execute the batch query
            job = self.simbad_tap.launch_job(adql_query)
            results = job.get_results()
            
            if results is None or len(results) == 0:
                logger.info("No SIMBAD objects found in batch query")
                return simbad_matches
            
            logger.info(f"Retrieved {len(results)} SIMBAD objects from batch query")
            
            # Process field-based results - need coordinate matching
            # Convert results to objects list
            simbad_objects = []
            for row in results:
                obj = {
                    'main_id': str(row['main_id']) if row['main_id'] is not np.ma.masked else '',
                    'otype': str(row['otype']) if 'otype' in row.colnames and row['otype'] is not np.ma.masked else '',
                    'sp_type': str(row['sp_type']) if 'sp_type' in row.colnames and row['sp_type'] is not np.ma.masked else '',
                    'ra': float(row['ra']) if 'ra' in row.colnames and row['ra'] is not np.ma.masked else np.nan,
                    'dec': float(row['dec']) if 'dec' in row.colnames and row['dec'] is not np.ma.masked else np.nan
                }
                if not (np.isnan(obj['ra']) or np.isnan(obj['dec'])):
                    simbad_objects.append(obj)
            
            logger.info(f"Matching {len(simbad_objects)} SIMBAD objects to {len(sky_coords)} sources")
            
            # Match each source to nearest SIMBAD object
            for i, coord in enumerate(sky_coords):
                source_ra = coord.ra.degree
                source_dec = coord.dec.degree
                
                best_match = None
                best_distance = float('inf')
                
                # Find nearest SIMBAD object within radius
                for simbad_obj in simbad_objects:
                    # Calculate proper angular separation using haversine formula
                    ra1 = np.radians(source_ra)
                    dec1 = np.radians(source_dec)
                    ra2 = np.radians(simbad_obj['ra'])
                    dec2 = np.radians(simbad_obj['dec'])
                    
                    # Haversine formula for great circle distance
                    delta_ra = ra2 - ra1
                    delta_dec = dec2 - dec1
                    
                    a = np.sin(delta_dec/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(delta_ra/2)**2
                    c = 2 * np.arcsin(np.sqrt(a))
                    distance_arcsec = np.degrees(c) * 3600.0
                    
                    if distance_arcsec <= radius_arcsec and distance_arcsec < best_distance:
                        best_match = simbad_obj.copy()
                        best_distance = distance_arcsec
                
                # Store the best match for this source
                if best_match:
                    simbad_matches[i] = {
                        'main_id': best_match['main_id'],
                        'otype': best_match['otype'],
                        'sp_type': best_match['sp_type'],
                        'rv_value': np.nan,
                        'distance_result': np.nan,
                        'match_distance_arcsec': best_distance,
                        'ra': best_match['ra'],
                        'dec': best_match['dec']
                    }
            
            successful_matches = len([m for m in simbad_matches if m['main_id']])
            logger.info(f"Improved batch SIMBAD matching completed: {successful_matches}/{len(sky_coords)} matches found")
            
            return simbad_matches
            
        except Exception as e:
            logger.error(f"Batch SIMBAD query failed: {e}")
            
            # Fallback to individual queries for smaller batches
            if len(sky_coords) <= 10:
                logger.info("Falling back to individual SIMBAD queries for small batch...")
                return self._query_simbad_fallback(sky_coords, radius_arcsec)
            else:
                logger.warning("Batch SIMBAD query failed and batch too large for individual fallback")
                return simbad_matches
    
    def _simbad_xmatch(self, sky_coords, radius_arcsec=5.0):
        """
        Perform efficient SIMBAD X-Match using single bulk TAP upload.
        Much faster than individual queries - replaces the slow individual SIMBAD loop.
        
        Args:
            sky_coords: List of SkyCoord objects
            radius_arcsec: Search radius in arcseconds
            
        Returns:
            List of SIMBAD match dictionaries with X-Match field names, one per input coordinate
        """
        try:
            if not sky_coords or len(sky_coords) == 0:
                return []
                
            logger.info(f"Starting SIMBAD bulk X-Match for {len(sky_coords)} sources (radius={radius_arcsec} arcsec)")
            
            # Use only the TAP-upload bulk X-Match approach
            matches = self._simbad_xmatch_bulk(sky_coords, radius_arcsec)
            
            successful_matches = len([m for m in matches if m.get('simbad_main_id')])
            logger.info(f"SIMBAD X-Match completed: {successful_matches}/{len(sky_coords)} sources matched")
            return matches
            
        except Exception as e:
            logger.error(f"SIMBAD X-Match failed: {e}")
            logger.info("Falling back to empty results")
            return self._create_empty_simbad_matches(len(sky_coords))
    
    def _simbad_xmatch_bulk(self, sky_coords, radius_arcsec=5.0):
        """
        Perform true SIMBAD X-Match using single bulk TAP upload.
        
        Args:
            sky_coords: List of SkyCoord objects
            radius_arcsec: Search radius in arcseconds
            
        Returns:
            List of SIMBAD match dictionaries with X-Match field names
        """
        try:
            from astropy.table import Table
            import numpy as np
            
            # Create astropy table directly with proper format (no CSV file needed)
            coords_table = Table()
            coords_table['source_id'] = np.arange(len(sky_coords), dtype=int)
            # Round coordinates to avoid decimal parsing issues in SIMBAD TAP
            coords_table['ra'] = np.array([round(coord.ra.degree, 6) for coord in sky_coords], dtype=float)  # Ensure degrees (0-360)
            coords_table['dec'] = np.array([round(coord.dec.degree, 6) for coord in sky_coords], dtype=float)  # Ensure degrees (-90 to +90)
            
            # Set column metadata to help TAP understand the format
            coords_table['ra'].unit = 'deg'
            coords_table['dec'].unit = 'deg'
            
            logger.info(f"Created upload table with {len(coords_table)} coordinates in degrees")
            
            # Build lean ADQL X-Match query with correct field names and aliases
            # Use radius in arcseconds and convert with division in ADQL to avoid decimal parsing issues
            radius_arcsec_int = int(radius_arcsec)
            radius_deg = float(radius_arcsec) / 3600.0  # e.g., 5" -> 0.0013888889

            xmatch_query = (
                "SELECT "
                "  u.source_id AS source_id, "
                "  u.ra  AS input_ra, "
                "  u.dec AS input_dec, "
                "  b.main_id AS simbad_main_id, "
                "  b.ra, b.dec, b.otype, b.sp_type, b.rvz_radvel, "
                "  b.pmra, b.pmdec, b.plx_value AS parallax, "
                "  DISTANCE(POINT('ICRS', u.ra, u.dec), POINT('ICRS', b.ra, b.dec)) AS distance_result "
                "FROM TAP_UPLOAD.\"upload_table\" AS u "
                "JOIN basic AS b "
                "  ON 1 = CONTAINS(POINT('ICRS', b.ra, b.dec), "
                f"                  CIRCLE('ICRS', u.ra, u.dec, {radius_deg:.8f})) "
                "ORDER BY source_id ASC, distance_result ASC"
            )
            
            logger.debug(f"ADQL query: {xmatch_query.strip()}")
            logger.debug(f"Upload table columns: {coords_table.colnames}")
            
            # Use astroquery SIMBAD TAP upload API - pass table as keyword argument
            from astroquery.simbad import Simbad
            results = Simbad.query_tap(xmatch_query, upload_table=coords_table)
            
            if results is None or len(results) == 0:
                logger.info("No X-Match results returned from TAP upload")
                return self._create_empty_simbad_matches(len(sky_coords))
            
            logger.info(f"TAP X-Match returned {len(results)} matches")
            
            # Process X-Match results and merge with input coordinates
            return self._process_tap_xmatch_results(results, len(sky_coords))
            
        except Exception as e:
            logger.error(f"SIMBAD TAP X-Match failed: {e}")
            logger.info("Falling back to single cone search + local matching")
            return self._simbad_cone_search_fallback(sky_coords, radius_arcsec)
    
    def _simbad_cone_search_fallback(self, sky_coords, radius_arcsec=5.0):
        """
        Fallback: Single SIMBAD cone search for whole field + local nearest-neighbor matching.
        No uploads, no MIME, no batching, no ADQL - just one query + astropy matching.
        
        Args:
            sky_coords: List of SkyCoord objects
            radius_arcsec: Search radius in arcseconds
            
        Returns:
            List of SIMBAD match dictionaries
        """
        try:
            import numpy as np
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            
            logger.info(f"SIMBAD cone search fallback for {len(sky_coords)} sources")
            
            # Smart partitioning: if sources are too spread out, query each individually
            # This is better than a single huge cone that returns irrelevant objects
            ras = [coord.ra.degree for coord in sky_coords]
            decs = [coord.dec.degree for coord in sky_coords]
            
            ra_span = max(ras) - min(ras)
            dec_span = max(decs) - min(decs)
            field_span = max(ra_span, dec_span)
            
            logger.debug(f"Field span: RA={ra_span:.3f}°, Dec={dec_span:.3f}°, max={field_span:.3f}°")
            
            # If field is too large (>5°), fall back to individual queries
            if field_span > 5.0:
                logger.info(f"Field span {field_span:.1f}° too large for cone search, using individual queries")
                return self._create_empty_simbad_matches(len(sky_coords))
            
            # For reasonably sized fields, do individual small cone searches
            all_simbad_objects = []
            for i, coord in enumerate(sky_coords):
                search_radius_deg = (radius_arcsec + 60) / 3600.0  # Add 60" padding for nearby objects
                
                individual_query = f"""
                SELECT main_id, ra, dec, otype, sp_type, rvz_radvel, pmra, pmdec, plx_value
                FROM basic
                WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {coord.ra.degree}, {coord.dec.degree}, {search_radius_deg}))
                """
                
                try:
                    job = self.simbad_tap.launch_job(individual_query)
                    objects = job.get_results()
                    if objects is not None and len(objects) > 0:
                        all_simbad_objects.extend(objects)
                        logger.debug(f"Source {i+1}: found {len(objects)} SIMBAD objects nearby")
                except Exception as e:
                    logger.debug(f"Source {i+1}: individual cone query failed: {e}")
            
            # Remove duplicates by main_id
            seen_ids = set()
            unique_objects = []
            for obj in all_simbad_objects:
                if obj['main_id'] not in seen_ids:
                    seen_ids.add(obj['main_id'])
                    unique_objects.append(obj)
            
            simbad_objects = unique_objects
            logger.info(f"Retrieved {len(simbad_objects)} unique SIMBAD objects from {len(sky_coords)} individual cone searches")
            
            if simbad_objects is None or len(simbad_objects) == 0:
                logger.info("No SIMBAD objects found in field")
                return self._create_empty_simbad_matches(len(sky_coords))
            
            logger.info(f"Retrieved {len(simbad_objects)} SIMBAD objects from field")
            
            # Debug: Show some sample coordinates 
            logger.debug(f"Sample input coords: {[(c.ra.degree, c.dec.degree) for c in sky_coords[:3]]}")
            logger.debug(f"Sample SIMBAD coords: {[(obj['ra'], obj['dec']) for obj in simbad_objects[:3]]}")
            
            # Local nearest-neighbor matching using astropy
            matches = self._local_nearest_neighbor_match(sky_coords, simbad_objects, radius_arcsec)
            
            successful_matches = len([m for m in matches if m['simbad_main_id']])
            logger.info(f"Local matching: {successful_matches}/{len(sky_coords)} sources matched")
            
            return matches
            
        except Exception as e:
            logger.error(f"SIMBAD cone search fallback failed: {e}")
            return self._create_empty_simbad_matches(len(sky_coords))
    
    def _local_nearest_neighbor_match(self, sky_coords, simbad_objects, radius_arcsec):
        """
        Perform local nearest-neighbor matching between input coordinates and SIMBAD objects.
        Pure astropy coordinate matching - no server-side processing.
        """
        try:
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            import numpy as np
            
            # Convert SIMBAD objects to SkyCoord for efficient matching
            simbad_coords = SkyCoord(
                ra=[float(obj['ra']) for obj in simbad_objects] * u.deg,
                dec=[float(obj['dec']) for obj in simbad_objects] * u.deg,
                frame='icrs'
            )
            
            matches = []
            radius = radius_arcsec * u.arcsec
            
            # For each input coordinate, find nearest SIMBAD object within radius
            for i, input_coord in enumerate(sky_coords):
                # Calculate separations to all SIMBAD objects
                separations = input_coord.separation(simbad_coords)
                
                # Debug: Show minimum separation for first few sources
                if i < 3:
                    min_sep = np.min(separations)
                    logger.debug(f"Source {i+1} at ({input_coord.ra.degree:.3f}, {input_coord.dec.degree:.3f}): min_sep={min_sep.to(u.arcsec):.1f} vs radius={radius:.1f}")
                
                # Find objects within search radius
                within_radius = separations <= radius
                
                if np.any(within_radius):
                    # Find the closest match
                    closest_idx = np.argmin(separations)
                    closest_separation = separations[closest_idx]
                    
                    if closest_separation <= radius:
                        # Format the match
                        simbad_obj = simbad_objects[closest_idx]
                        match_dict = {
                            'simbad_main_id': str(simbad_obj['main_id']) if simbad_obj['main_id'] is not None else '',
                            'ra': float(simbad_obj['ra']) if simbad_obj['ra'] is not None else np.nan,
                            'dec': float(simbad_obj['dec']) if simbad_obj['dec'] is not None else np.nan,
                            'otype': str(simbad_obj['otype']) if simbad_obj['otype'] is not None else '',
                            'sp_type': str(simbad_obj['sp_type']) if simbad_obj['sp_type'] is not None else '',
                            'rvz_radvel': float(simbad_obj['rvz_radvel']) if simbad_obj['rvz_radvel'] is not None else np.nan,
                            'pmra': float(simbad_obj['pmra']) if simbad_obj['pmra'] is not None else np.nan,
                            'pmdec': float(simbad_obj['pmdec']) if simbad_obj['pmdec'] is not None else np.nan,
                            'parallax': float(simbad_obj['plx_value']) if simbad_obj['plx_value'] is not None else np.nan,
                            'separation_arcsec': closest_separation.arcsec,
                            'distance_pc': (1000.0 / float(simbad_obj['plx_value'])) if (simbad_obj['plx_value'] is not None and float(simbad_obj['plx_value']) > 0) else np.nan
                        }
                        matches.append(match_dict)
                    else:
                        matches.append(self._create_empty_simbad_match())
                else:
                    matches.append(self._create_empty_simbad_match())
            
            return matches
            
        except Exception as e:
            logger.error(f"Local nearest-neighbor matching failed: {e}")
            return self._create_empty_simbad_matches(len(sky_coords))
    
    # Removed old per-object processing functions - using only TAP bulk upload
    
    # Removed old astroquery formatting function - using only TAP bulk upload
    
    # Removed legacy CSV uploader - using direct astropy Table upload
    
    def _process_tap_xmatch_results(self, results, num_sources):
        """Process TAP X-Match results and create match list with one entry per source."""
        # Initialize empty matches for all sources
        matches = self._create_empty_simbad_matches(num_sources)
        
        # Group results by source_id and take the closest match (results are already sorted by distance_result)
        for row in results:
            source_id = int(row['source_id'])
            
            if source_id < len(matches):
                # Convert TAP X-Match fields with standardized units
                distance_result = float(row['distance_result']) if row['distance_result'] is not None else np.nan
                separation_arcsec = distance_result * 3600.0 if not np.isnan(distance_result) else np.nan
                
                match_dict = {
                    'simbad_main_id': str(row['simbad_main_id']) if row['simbad_main_id'] is not None else '',
                    'ra': float(row['ra']) if row['ra'] is not None else np.nan,
                    'dec': float(row['dec']) if row['dec'] is not None else np.nan,
                    'otype': str(row['otype']) if row['otype'] is not None else '',
                    'sp_type': str(row['sp_type']) if row['sp_type'] is not None else '',
                    'rvz_radvel': float(row['rvz_radvel']) if row['rvz_radvel'] is not None else np.nan,
                    'pmra': float(row['pmra']) if row['pmra'] is not None else np.nan,
                    'pmdec': float(row['pmdec']) if row['pmdec'] is not None else np.nan,
                    'parallax': float(row['parallax']) if row['parallax'] is not None else np.nan,
                    'separation_arcsec': separation_arcsec,
                    'distance_pc': (1000.0 / float(row['parallax'])) if (row['parallax'] is not None and float(row['parallax']) > 0) else np.nan
                }
                
                # Only update if this is the first match or closer than existing match
                if not matches[source_id]['simbad_main_id'] or separation_arcsec < matches[source_id].get('separation_arcsec', float('inf')):
                    matches[source_id] = match_dict
        
        return matches
    
    # Removed fallback function - using only TAP bulk upload
    
    def _create_empty_simbad_match(self):
        """Create an empty SIMBAD match dictionary with standardized field names."""
        return {
            'simbad_main_id': '',
            'ra': np.nan,
            'dec': np.nan,
            'otype': '',
            'sp_type': '',
            'rvz_radvel': np.nan,
            'pmra': np.nan,
            'pmdec': np.nan,
            'parallax': np.nan,
            'separation_arcsec': np.nan,
            'distance_pc': np.nan
        }
    
    def _simbad_xmatch_batch(self, batch_coords, radius_arcsec, start_index):
        """
        Process a single batch of coordinates using efficient area query + filtering.
        """
        if not batch_coords:
            return []
        
        # Calculate area boundaries for the batch
        ras = [coord.ra.degree for coord in batch_coords]
        decs = [coord.dec.degree for coord in batch_coords]
        
        ra_min, ra_max = min(ras) - 0.02, max(ras) + 0.02  # Add padding
        dec_min, dec_max = min(decs) - 0.02, max(decs) + 0.02
        
        
        # Query SIMBAD for all objects in the region (simplified schema)
        region_query = f"""
        SELECT b.main_id, b.ra, b.dec, b.pmra, b.pmdec, b.plx_value as parallax,
               otypes.otype as object_type
        FROM basic b
        LEFT JOIN otypes ON b.oid = otypes.oidref
        WHERE b.ra BETWEEN {ra_min} AND {ra_max}
          AND b.dec BETWEEN {dec_min} AND {dec_max}
        """
        
        job = self.simbad_tap.launch_job(region_query)
        results = job.get_results()
        
        if results is None or len(results) == 0:
            logger.debug(f"No SIMBAD objects found in region RA: {ra_min:.3f}-{ra_max:.3f}, Dec: {dec_min:.3f}-{dec_max:.3f}")
            return self._create_empty_simbad_matches(len(batch_coords))
        
        
        # Match each coordinate with the nearest SIMBAD object within radius
        batch_matches = []
        radius_deg = radius_arcsec / 3600.0
        
        for i, coord in enumerate(batch_coords):
            best_match = None
            best_separation = float('inf')
            matches_within_radius = 0
            
            for row in results:
                simbad_ra = float(row['ra'])
                simbad_dec = float(row['dec'])
                
                # Calculate proper spherical separation accounting for coordinate system
                # Use astropy SkyCoord for accurate separation calculation
                simbad_coord = SkyCoord(ra=simbad_ra*u.deg, dec=simbad_dec*u.deg, frame='icrs')
                separation_angle = coord.separation(simbad_coord)
                separation = separation_angle.degree
                
                if separation <= radius_deg:
                    matches_within_radius += 1
                    if separation < best_separation:
                        best_match = row
                        best_separation = separation
            
            if best_match is not None:
                match_dict = self._format_simbad_match_from_row(best_match, best_separation * 3600)
                batch_matches.append(match_dict)
            else:
                batch_matches.append(self._create_empty_simbad_match())
        
        return batch_matches
    
    def _format_simbad_match_from_row(self, row, separation_arcsec):
        """Format a SIMBAD result row into match dictionary format."""
        try:
            return {
                'main_id': str(row['main_id']) if row['main_id'] is not None else '',
                'ra': float(row['ra']) if row['ra'] is not None else np.nan,
                'dec': float(row['dec']) if row['dec'] is not None else np.nan,
                'object_type': str(row['object_type']) if row['object_type'] is not None else '',
                'v_mag': np.nan,  # Magnitude data not included in simplified query
                'b_mag': np.nan,  # Magnitude data not included in simplified query
                'pmra': float(row['pmra']) if row['pmra'] is not None else np.nan,
                'pmdec': float(row['pmdec']) if row['pmdec'] is not None else np.nan,
                'parallax': float(row['parallax']) if row['parallax'] is not None else np.nan,
                'separation_arcsec': float(separation_arcsec),
                'distance_pc': (1000.0 / float(row['parallax'])) if (row['parallax'] is not None and float(row['parallax']) > 0) else np.nan
            }
        except Exception as e:
            logger.warning(f"Error formatting SIMBAD match: {e}")
            return self._create_empty_simbad_match()
    
    def _process_xmatch_results(self, results, num_sources):
        """
        Process X-Match results and create a match list with one entry per source.
        Takes the best match (closest) for each source.
        """
        # Initialize empty matches for all sources
        simbad_matches = self._create_empty_simbad_matches(num_sources)
        
        # Group results by source_id and take the best match for each
        current_source_id = None
        best_match = None
        
        for row in results:
            source_id = int(row['source_id'])
            separation = float(row['separation_arcsec'])
            
            # If this is a new source or a better match for current source
            if source_id != current_source_id:
                # Save previous best match if we had one
                if current_source_id is not None and best_match is not None:
                    simbad_matches[current_source_id] = self._format_simbad_match(best_match)
                
                # Start new source
                current_source_id = source_id
                best_match = row
            else:
                # Same source - check if this is a better match
                if separation < float(best_match['separation_arcsec']):
                    best_match = row
        
        # Don't forget the last match
        if current_source_id is not None and best_match is not None:
            simbad_matches[current_source_id] = self._format_simbad_match(best_match)
        
        return simbad_matches
    
    def _format_simbad_match(self, row):
        """Format a SIMBAD X-Match result row into the expected match dictionary format."""
        try:
            return {
                'main_id': str(row['main_id']) if row['main_id'] is not None else '',
                'ra': float(row['simbad_ra']) if row['simbad_ra'] is not None else np.nan,
                'dec': float(row['simbad_dec']) if row['simbad_dec'] is not None else np.nan,
                'object_type': str(row['object_type']) if row['object_type'] is not None else '',
                'v_mag': float(row['v_mag']) if row['v_mag'] is not None else np.nan,
                'b_mag': float(row['b_mag']) if row['b_mag'] is not None else np.nan,
                'pmra': float(row['pmra']) if row['pmra'] is not None else np.nan,
                'pmdec': float(row['pmdec']) if row['pmdec'] is not None else np.nan,
                'parallax': float(row['parallax']) if row['parallax'] is not None else np.nan,
                'separation_arcsec': float(row['separation_arcsec']),
                'distance_pc': (1000.0 / float(row['parallax'])) if (row['parallax'] is not None and float(row['parallax']) > 0) else np.nan
            }
        except Exception as e:
            logger.warning(f"Error formatting SIMBAD match: {e}")
            return self._create_empty_simbad_match()
    
    def _create_empty_simbad_matches(self, count):
        """Create a list of empty SIMBAD match dictionaries."""
        return [self._create_empty_simbad_match() for _ in range(count)]
    
    def _query_simbad_fallback(self, sky_coords, radius_arcsec=5.0):
        """
        Fallback to individual SIMBAD queries when batch query fails.
        Only used for small batches to avoid connection overload.
        """
        simbad_matches = []
        
        for i, coord in enumerate(sky_coords):
            ra = coord.ra.degree  
            dec = coord.dec.degree
            
            # Use the individual query method with retry logic
            simbad_match = self._query_simbad_tap(ra, dec, radius_arcsec, max_retries=2)
            simbad_matches.append(simbad_match)
            
            # Longer delay between individual queries
            if i < len(sky_coords) - 1:
                time.sleep(1.0)  # 1 second delay for fallback mode
        
        logger.info(f"Fallback SIMBAD queries completed for {len(sky_coords)} sources")
        return simbad_matches
    
    def _query_gaia_optimized(self, detected_sources_coords, field_radius_deg=0.25):
        """
        Optimized Gaia query: Single async CIRCLE query with RUWE/magnitude filtering.
        Replaces 189 tiny cone searches with one field query + local matching.
        
        Args:
            detected_sources_coords: List of SkyCoord objects for detected sources
            field_radius_deg: Field radius in degrees for single query
            
        Returns:
            dict: Mapping from source index to Gaia data, or None if query fails
        """
        try:
            if not detected_sources_coords:
                return {}
            
            logger.info(f"Starting optimized Gaia query for {len(detected_sources_coords)} sources")
            
            # Calculate field center from detected sources
            ras = [coord.ra.degree for coord in detected_sources_coords]
            decs = [coord.dec.degree for coord in detected_sources_coords]
            
            center_ra = np.mean(ras)
            center_dec = np.mean(decs)
            
            # Adaptive radius based on source distribution
            ra_span = max(ras) - min(ras)
            dec_span = max(decs) - min(decs)
            auto_radius = max(ra_span, dec_span) / 2.0 + 0.05  # Add 0.05° padding
            query_radius = max(field_radius_deg, auto_radius)
            
            logger.info(f"Field center: ({center_ra:.4f}°, {center_dec:.4f}°), radius: {query_radius:.4f}°")
            
            # Single optimized ADQL query with RUWE and magnitude filtering
            query = f"""
            SELECT 
                source_id, ra, dec, 
                phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
                bp_rp, parallax, pmra, pmdec, ruwe,
                DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {center_ra}, {center_dec})) AS center_distance
            FROM gaiadr3.gaia_source 
            WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {center_ra}, {center_dec}, {query_radius})) = 1
                AND phot_g_mean_mag < 19.0
                AND phot_g_mean_mag IS NOT NULL
                AND (ruwe IS NULL OR ruwe < 1.4)
                AND (phot_bp_mean_mag IS NOT NULL AND phot_rp_mean_mag IS NOT NULL)
                AND parallax_over_error > 3.0
            ORDER BY phot_g_mean_mag ASC
            """
            
            logger.debug(f"Optimized Gaia query: {query}")
            
            # Execute single async query
            start_time = time.time()
            job = Gaia.launch_job_async(query)
            gaia_sources = job.get_results()
            query_time = time.time() - start_time
            
            if gaia_sources is None or len(gaia_sources) == 0:
                logger.warning("No Gaia sources returned from optimized query")
                return {}
            
            logger.info(f"Optimized Gaia query completed in {query_time:.2f}s: {len(gaia_sources)} sources")
            
            # Convert Gaia results to SkyCoord for efficient matching
            gaia_coords = SkyCoord(
                ra=gaia_sources['ra'],
                dec=gaia_sources['dec'],
                unit='deg'
            )
            
            # Perform local nearest-neighbor matching
            matches = {}
            match_radius = 3.0 * u.arcsec  # 3 arcsec matching radius
            
            start_match_time = time.time()
            for i, source_coord in enumerate(detected_sources_coords):
                # Find closest Gaia source
                separations = source_coord.separation(gaia_coords)
                min_sep_idx = separations.argmin()
                min_separation = separations[min_sep_idx]
                
                if min_separation <= match_radius:
                    gaia_row = gaia_sources[min_sep_idx]
                    matches[i] = {
                        'gaia_source_id': str(gaia_row['source_id']),
                        'gaia_ra': float(gaia_row['ra']),
                        'gaia_dec': float(gaia_row['dec']),
                        'phot_g_mean_mag': float(gaia_row['phot_g_mean_mag']) if gaia_row['phot_g_mean_mag'] is not np.ma.masked else np.nan,
                        'phot_bp_mean_mag': float(gaia_row['phot_bp_mean_mag']) if gaia_row['phot_bp_mean_mag'] is not np.ma.masked else np.nan,
                        'phot_rp_mean_mag': float(gaia_row['phot_rp_mean_mag']) if gaia_row['phot_rp_mean_mag'] is not np.ma.masked else np.nan,
                        'bp_rp': float(gaia_row['bp_rp']) if gaia_row['bp_rp'] is not np.ma.masked else np.nan,
                        'parallax': float(gaia_row['parallax']) if gaia_row['parallax'] is not np.ma.masked else np.nan,
                        'pmra': float(gaia_row['pmra']) if gaia_row['pmra'] is not np.ma.masked else np.nan,
                        'pmdec': float(gaia_row['pmdec']) if gaia_row['pmdec'] is not np.ma.masked else np.nan,
                        'ruwe': float(gaia_row['ruwe']) if gaia_row['ruwe'] is not np.ma.masked else np.nan,
                        'match_distance_arcsec': float(min_separation.arcsec)
                    }
            
            match_time = time.time() - start_match_time
            logger.info(f"Local matching completed in {match_time:.2f}s: {len(matches)}/{len(detected_sources_coords)} matched")
            
            return matches
            
        except Exception as e:
            logger.error(f"Optimized Gaia query failed: {e}")
            return {}

    def _query_gaia_field(self, center_ra, center_dec, radius_deg):
        """Query Gaia DR3 for sources in the field with overload protection and tiling."""
        try:
            # Configuration for query limits
            MAX_RADIUS_DEG = 0.6  # Maximum cone search radius
            MAX_ROWS = 150000     # TOP limit for queries
            MAGNITUDE_LIMIT = 18.5  # Brightness limit to reduce results
            TILE_SIZE_DEG = 0.5   # Size of individual tiles for large fields
            
            logger.info(f"Starting Gaia field query: center=({center_ra:.4f}, {center_dec:.4f}), radius={radius_deg:.4f}°")
            
            # Check if field is too large and needs tiling
            if radius_deg <= MAX_RADIUS_DEG:
                # Single query for small fields
                return self._query_gaia_single(center_ra, center_dec, radius_deg, MAX_ROWS, MAGNITUDE_LIMIT)
            else:
                # Tile large fields into smaller queries
                logger.info(f"Field radius {radius_deg:.4f}° exceeds max {MAX_RADIUS_DEG}°, using tiled approach")
                return self._query_gaia_tiled(center_ra, center_dec, radius_deg, TILE_SIZE_DEG, MAX_ROWS, MAGNITUDE_LIMIT)
                
        except Exception as e:
            logger.error(f"Error in Gaia field query: {str(e)}")
            return None
    
    def _query_gaia_single(self, center_ra, center_dec, radius_deg, max_rows, mag_limit, max_retries=3):
        """Single Gaia query with retry logic for server errors."""
        
        for attempt in range(max_retries):
            try:
                # Calculate retry parameters (smaller radius/stricter limits on retry)
                retry_radius = radius_deg * (0.8 ** attempt)  # Reduce radius by 20% each retry
                retry_mag_limit = mag_limit - (0.5 * attempt)  # Stricter magnitude limit
                retry_max_rows = int(max_rows * (0.7 ** attempt))  # Fewer rows
                
                logger.debug(f"Gaia query attempt {attempt + 1}: radius={retry_radius:.4f}°, mag_limit={retry_mag_limit:.1f}, max_rows={retry_max_rows}")
                
                # Optimized ADQL query with RUWE filter, magnitude cuts, and minimal columns
                query = f"""
                SELECT TOP {retry_max_rows}
                    source_id, ra, dec, 
                    phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
                    bp_rp, parallax, pmra, pmdec, ruwe
                FROM gaiadr3.gaia_source 
                WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {center_ra}, {center_dec}, {retry_radius})) = 1
                    AND phot_g_mean_mag < {retry_mag_limit}
                    AND phot_g_mean_mag IS NOT NULL
                    AND (ruwe IS NULL OR ruwe < 1.4)
                    AND (phot_bp_mean_mag IS NOT NULL AND phot_rp_mean_mag IS NOT NULL)
                ORDER BY phot_g_mean_mag ASC
                """
                
                logger.debug(f"Gaia ADQL query: {query}")
                job = Gaia.launch_job_async(query)
                results = job.get_results()
                
                if len(results) == 0:
                    logger.warning(f"No Gaia sources found (attempt {attempt + 1})")
                    if attempt == max_retries - 1:
                        return None
                    continue
                
                logger.info(f"Gaia query successful: {len(results)} sources retrieved")
                return results
                
            except Exception as e:
                error_msg = str(e).lower()
                logger.warning(f"Gaia query attempt {attempt + 1} failed: {str(e)}")
                
                # Check for server overload indicators
                if '500' in error_msg or 'server error' in error_msg or 'timeout' in error_msg or 'overload' in error_msg:
                    if attempt < max_retries - 1:
                        delay = (2 ** attempt) + random.uniform(0.5, 2.0)
                        logger.info(f"Server error detected, retrying with smaller parameters in {delay:.1f}s")
                        time.sleep(delay)
                        continue
                else:
                    # Non-retryable error
                    logger.error(f"Non-retryable Gaia query error: {str(e)}")
                    return None
        
        logger.error("All Gaia query attempts failed")
        return None
    
    def _query_gaia_tiled(self, center_ra, center_dec, radius_deg, tile_size, max_rows_per_tile, mag_limit):
        """Query Gaia using tiling approach for large fields."""
        try:
            # Calculate grid of tiles to cover the circular field
            tiles = self._calculate_field_tiles(center_ra, center_dec, radius_deg, tile_size)
            logger.info(f"Tiling large field into {len(tiles)} tiles")
            
            all_results = []
            successful_tiles = 0
            
            for i, (tile_ra, tile_dec, tile_radius) in enumerate(tiles):
                logger.debug(f"Querying tile {i+1}/{len(tiles)}: center=({tile_ra:.4f}, {tile_dec:.4f}), radius={tile_radius:.4f}°")
                
                # Query this tile
                tile_results = self._query_gaia_single(tile_ra, tile_dec, tile_radius, max_rows_per_tile, mag_limit)
                
                if tile_results is not None and len(tile_results) > 0:
                    all_results.append(tile_results)
                    successful_tiles += 1
                    logger.debug(f"Tile {i+1} returned {len(tile_results)} sources")
                else:
                    logger.debug(f"Tile {i+1} returned no sources")
                
                # Small delay between tile queries to avoid overwhelming server
                if i < len(tiles) - 1:
                    time.sleep(0.5)
            
            if successful_tiles == 0:
                logger.warning("No tiles returned Gaia sources")
                return None
            
            # Combine all results
            if len(all_results) == 1:
                combined_results = all_results[0]
            else:
                from astropy.table import vstack
                combined_results = vstack(all_results)
            
            # Remove duplicates (sources appearing in multiple tiles)
            combined_results = self._remove_duplicate_sources(combined_results, 'source_id')
            
            logger.info(f"Tiled Gaia query complete: {len(combined_results)} unique sources from {successful_tiles}/{len(tiles)} tiles")
            return combined_results
            
        except Exception as e:
            logger.error(f"Error in tiled Gaia query: {str(e)}")
            return None
    
    def _calculate_field_tiles(self, center_ra, center_dec, radius_deg, tile_size):
        """Calculate grid of tiles to cover a circular field."""
        tiles = []
        
        # Simple square tiling that covers the circular area
        # Calculate bounding box
        ra_min = center_ra - radius_deg
        ra_max = center_ra + radius_deg
        dec_min = center_dec - radius_deg
        dec_max = center_dec + radius_deg
        
        # Generate tile centers
        ra_current = ra_min
        while ra_current <= ra_max:
            dec_current = dec_min
            while dec_current <= dec_max:
                # Check if tile center is within the original circular field
                separation = np.sqrt((ra_current - center_ra)**2 + (dec_current - center_dec)**2)
                
                if separation <= radius_deg + (tile_size * 0.5):  # Include tiles that overlap the edge
                    # Use smaller radius for individual tiles, but ensure coverage
                    tile_radius = min(tile_size * 0.6, radius_deg)  # 0.6 factor for overlap
                    tiles.append((ra_current, dec_current, tile_radius))
                
                dec_current += tile_size * 0.8  # 80% step for overlap
            ra_current += tile_size * 0.8  # 80% step for overlap
        
        return tiles
    
    def _remove_duplicate_sources(self, table, id_column):
        """Remove duplicate sources based on ID column."""
        try:
            # Convert to pandas for easy duplicate removal, then back to astropy
            df = table.to_pandas()
            df_unique = df.drop_duplicates(subset=[id_column], keep='first')
            
            # Convert back to astropy table
            unique_table = Table.from_pandas(df_unique)
            
            logger.debug(f"Removed {len(table) - len(unique_table)} duplicate sources")
            return unique_table
            
        except Exception as e:
            logger.warning(f"Error removing duplicates, returning original table: {str(e)}")
            return table
    
    def _query_vizier_catalogs_field(self, center_ra, center_dec, radius_deg):
        """Query VizieR B catalogs for sources in the field (NO SIMBAD)."""
        try:
            # Configure VizieR to query multiple catalogs
            vizier = Vizier(row_limit=1000)
            vizier.columns = ['*']  # Get all available columns
            
            # Convert to SkyCoord for region query
            center_coord = SkyCoord(center_ra, center_dec, unit='deg')
            radius = radius_deg * u.deg
            
            logger.debug(f"VizieR B catalogs query: center={center_coord}, radius={radius}")
            
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
            
            all_results = []
            
            for catalog in catalogs_to_try:
                try:
                    logger.debug(f"Querying VizieR catalog: {catalog}")
                    tables = vizier.query_region(center_coord, radius=radius, catalog=catalog)
                    
                    if tables and len(tables) > 0:
                        # Get first catalog table
                        table = tables[0]
                        if len(table) > 0:
                            # Add catalog source information
                            table['catalog_source'] = [catalog] * len(table)
                            all_results.append(table)
                            logger.debug(f"Found {len(table)} sources in {catalog}")
                
                except Exception as cat_error:
                    logger.debug(f"Catalog {catalog} query failed: {cat_error}")
                    continue
            
            if not all_results:
                logger.warning("No VizieR catalog sources found in field")
                return None
            
            # Combine results from all catalogs
            from astropy.table import vstack
            combined_results = vstack(all_results)
            
            logger.info(f"VizieR B catalogs query returned {len(combined_results)} total sources from {len(all_results)} catalogs")
            
            # Convert coordinates to standard format
            combined_results = self._convert_vizier_coordinates(combined_results)
            
            return combined_results
            
        except Exception as e:
            logger.error(f"Error querying VizieR B catalogs: {str(e)}")
            return None
    
    def _convert_vizier_coordinates(self, vizier_table):
        """Convert VizieR B catalog RA/DEC with robust column detection and format handling."""
        try:
            c = vizier_table.colnames
            logger.debug(f"Available VizieR columns: {c}")
            
            # Determine coordinate columns - VizieR uses different naming conventions
            ra_key = None
            dec_key = None
            
            # VizieR B catalog coordinate column preferences
            for potential_ra in ['_RA', 'RA_ICRS', 'RAJ2000', 'RA_deg', 'RA_d', 'ra', 'RA']:
                if potential_ra in c:
                    ra_key = potential_ra
                    break
            
            for potential_dec in ['_DE', 'DE_ICRS', 'DEJ2000', 'DEC_deg', 'DEC_d', 'dec', 'DEC']:
                if potential_dec in c:
                    dec_key = potential_dec
                    break
            
            if ra_key is None or dec_key is None:
                logger.error(f"Could not find RA/DEC columns in VizieR data. Available: {c}")
                return vizier_table
            
            logger.debug(f"Using VizieR coordinate columns: {ra_key}, {dec_key}")
            
            # VizieR typically returns coordinates in degrees already
            try:
                # Check if data is numeric (VizieR usually provides degrees directly)
                ra_sample = vizier_table[ra_key][0]
                if isinstance(ra_sample, (int, float)) and not isinstance(ra_sample, str):
                    # Already numeric degrees
                    vizier_table['RA_deg'] = vizier_table[ra_key]
                    vizier_table['DEC_deg'] = vizier_table[dec_key]
                    logger.debug(f"Used numeric VizieR coordinates directly (degrees)")
                else:
                    # Parse as coordinate strings 
                    try:
                        # Try degrees first (most common for VizieR)
                        coords = SkyCoord(vizier_table[ra_key], vizier_table[dec_key], unit=(u.deg, u.deg))
                        vizier_table['RA_deg'] = coords.ra.deg
                        vizier_table['DEC_deg'] = coords.dec.deg
                        logger.debug(f"Converted {len(vizier_table)} VizieR coordinates from degree strings")
                    except Exception:
                        # Fallback to hourangle/degree (less common for VizieR)
                        try:
                            coords = SkyCoord(vizier_table[ra_key], vizier_table[dec_key], unit=(u.hourangle, u.deg))
                            vizier_table['RA_deg'] = coords.ra.deg
                            vizier_table['DEC_deg'] = coords.dec.deg
                            logger.debug(f"Converted {len(vizier_table)} VizieR coordinates from hourangle/deg strings")
                        except Exception as coord_error:
                            logger.error(f"Failed to parse VizieR coordinates: {coord_error}")
                            # Use direct assignment as fallback
                            try:
                                vizier_table['RA_deg'] = [float(ra) for ra in vizier_table[ra_key]]
                                vizier_table['DEC_deg'] = [float(dec) for dec in vizier_table[dec_key]]
                                logger.debug("Used direct float conversion for VizieR coordinates")
                            except Exception:
                                logger.error("All coordinate conversion methods failed")
                                return vizier_table
            except Exception as parse_error:
                logger.warning(f"Failed to parse VizieR coordinates with SkyCoord, trying direct conversion: {str(parse_error)}")
                # Final fallback: try direct numeric conversion
                try:
                    vizier_table['RA_deg'] = [float(ra) for ra in vizier_table[ra_key]]
                    vizier_table['DEC_deg'] = [float(dec) for dec in vizier_table[dec_key]]
                    logger.debug("Used direct numeric conversion as fallback for VizieR")
                except Exception as final_error:
                    logger.error(f"All VizieR coordinate parsing methods failed: {str(final_error)}")
                    return vizier_table
            
            # Validate the converted coordinates and log statistics
            if 'RA_deg' in vizier_table.colnames and 'DEC_deg' in vizier_table.colnames:
                valid_coords = sum(1 for ra, dec in zip(vizier_table['RA_deg'], vizier_table['DEC_deg']) 
                                 if not (np.isnan(float(ra)) or np.isnan(float(dec))))
                
                logger.info(f"Successfully converted {valid_coords}/{len(vizier_table)} VizieR B catalog coordinates")
                if valid_coords > 0:
                    logger.debug(f"RA range: {np.nanmin(vizier_table['RA_deg']):.4f} to {np.nanmax(vizier_table['RA_deg']):.4f}")
                    logger.debug(f"DEC range: {np.nanmin(vizier_table['DEC_deg']):.4f} to {np.nanmax(vizier_table['DEC_deg']):.4f}")
                
            return vizier_table
            
        except Exception as e:
            logger.error(f"Error in VizieR B catalog coordinate conversion: {str(e)}")
            logger.debug(f"VizieR table columns: {vizier_table.colnames}")
            logger.debug(f"Table length: {len(vizier_table)}")
            return vizier_table
    
    def _nearest_neighbor_match(self, source_coords, catalog, catalog_type):
        """Perform nearest-neighbor matching between sources and catalog with proper coordinate handling."""
        try:
            if catalog is None or len(catalog) == 0:
                logger.debug(f"No {catalog_type} catalog data available for matching")
                return [{}] * len(source_coords)
                
            # Extract coordinates from catalog with proper handling
            if catalog_type == 'gaia':
                catalog_ras = catalog['ra']
                catalog_decs = catalog['dec']
                logger.debug(f"Using Gaia numeric coordinates: {len(catalog)} sources")
            else:  # SIMBAD
                # Try to find coordinate columns with tolerance to various names
                ra_col = None
                dec_col = None
                
                # Try to find RA column with fallback order
                for ra_candidate in ['RA_deg', 'ra', 'RA', '_RAJ2000', 'RAJ2000']:
                    if ra_candidate in catalog.colnames:
                        ra_col = ra_candidate
                        break
                
                # Try to find DEC column with fallback order  
                for dec_candidate in ['DEC_deg', 'dec', 'DEC', '_DEJ2000', 'DEJ2000']:
                    if dec_candidate in catalog.colnames:
                        dec_col = dec_candidate
                        break
                
                if ra_col is None or dec_col is None:
                    logger.error(f"Could not find SIMBAD coordinate columns. Available: {catalog.colnames}")
                    return [{}] * len(source_coords)
                
                catalog_ras = catalog[ra_col]
                catalog_decs = catalog[dec_col]
                logger.debug(f"Using SIMBAD coordinates from {ra_col}/{dec_col}: {len(catalog)} sources")
            
            # Validate coordinate data
            valid_coords = []
            valid_indices = []
            
            for i, (ra, dec) in enumerate(zip(catalog_ras, catalog_decs)):
                try:
                    # Check if coordinates are valid numbers
                    ra_val = float(ra)
                    dec_val = float(dec)
                    
                    if not (np.isnan(ra_val) or np.isnan(dec_val)):
                        valid_coords.append([ra_val, dec_val])
                        valid_indices.append(i)
                except (ValueError, TypeError):
                    logger.debug(f"Invalid {catalog_type} coordinates at index {i}: RA={ra}, DEC={dec}")
                    continue
            
            if len(valid_coords) == 0:
                logger.warning(f"No valid coordinates found in {catalog_type} catalog")
                return [{}] * len(source_coords)
            
            # Add guardrails for KD-tree building
            try:
                catalog_coords = np.array(valid_coords)
                
                # Assert we have numeric arrays before building KD-tree
                if not np.issubdtype(catalog_coords.dtype, np.number):
                    logger.error(f"Non-numeric coordinates detected for {catalog_type} KD-tree")
                    logger.error(f"Available {catalog_type} columns: {catalog.colnames}")
                    if len(valid_coords) > 0:
                        logger.error(f"First 5 coordinate samples: {valid_coords[:5]}")
                    return [{}] * len(source_coords)
                
                if catalog_coords.shape[1] != 2:
                    logger.error(f"Invalid coordinate array shape for {catalog_type}: {catalog_coords.shape}")
                    return [{}] * len(source_coords)
                
                # Check for any remaining NaN or infinite values
                finite_mask = np.isfinite(catalog_coords).all(axis=1)
                if not finite_mask.all():
                    n_invalid = (~finite_mask).sum()
                    logger.warning(f"Removing {n_invalid} non-finite coordinates from {catalog_type}")
                    catalog_coords = catalog_coords[finite_mask]
                    valid_indices = [vi for vi, mask in zip(valid_indices, finite_mask) if mask]
                    
                    if len(catalog_coords) == 0:
                        logger.error(f"No finite coordinates remaining for {catalog_type} KD-tree")
                        return [{}] * len(source_coords)
                
                # Build KD-Tree for efficient nearest-neighbor search
                tree = cKDTree(catalog_coords)
                
                logger.info(f"Successfully built {catalog_type} KD-tree with {len(catalog_coords)} coordinates")
                
            except Exception as kdtree_error:
                logger.error(f"Failed to build KD-tree for {catalog_type}: {str(kdtree_error)}")
                logger.error(f"Available {catalog_type} columns: {catalog.colnames}")
                if len(valid_coords) > 0:
                    logger.error(f"First 5 coordinate samples: {valid_coords[:5]}")
                return [{}] * len(source_coords)
            
            matches = []
            match_radius_deg = 5.0 / 3600.0  # 5 arcsec in degrees
            
            logger.debug(f"Performing nearest-neighbor matching: {len(source_coords)} sources vs {len(valid_coords)} {catalog_type} entries")
            
            for coord in source_coords:
                source_pos = np.array([coord.ra.degree, coord.dec.degree])
                
                # Find nearest neighbor
                distance, tree_index = tree.query(source_pos)
                
                # Check if match is within tolerance
                if distance <= match_radius_deg:
                    # Get the original catalog index
                    catalog_index = valid_indices[tree_index]
                    
                    if catalog_type == 'gaia':
                        match_data = self._extract_gaia_data(catalog[catalog_index])
                    else:
                        match_data = self._extract_vizier_catalog_data(catalog[catalog_index])
                    
                    # Add match distance for debugging
                    match_data['match_distance_arcsec'] = distance * 3600.0
                    matches.append(match_data)
                else:
                    matches.append({})  # No match found
            
            num_matches = sum(1 for m in matches if m)
            logger.info(f"Found {num_matches}/{len(source_coords)} matches in {catalog_type} (within 5 arcsec)")
            return matches
            
        except Exception as e:
            logger.error(f"Error in nearest-neighbor matching for {catalog_type}: {str(e)}")
            logger.debug(f"Catalog type: {catalog_type}, catalog shape: {len(catalog) if catalog is not None else 'None'}")
            if catalog is not None:
                logger.debug(f"Catalog columns: {catalog.colnames}")
            return [{}] * len(source_coords)
    
    def _extract_gaia_data(self, gaia_row):
        """Extract data from Gaia catalog row."""
        try:
            return {
                'source_id': str(gaia_row['source_id']),
                'ra': float(gaia_row['ra']),
                'dec': float(gaia_row['dec']),
                'phot_g_mean_mag': float(gaia_row['phot_g_mean_mag']) if gaia_row['phot_g_mean_mag'] is not None else np.nan,
                'phot_bp_mean_mag': float(gaia_row['phot_bp_mean_mag']) if gaia_row['phot_bp_mean_mag'] is not None else np.nan,
                'phot_rp_mean_mag': float(gaia_row['phot_rp_mean_mag']) if gaia_row['phot_rp_mean_mag'] is not None else np.nan,
                'bp_rp': float(gaia_row['bp_rp']) if gaia_row['bp_rp'] is not None else np.nan,
                'parallax': float(gaia_row['parallax']) if gaia_row['parallax'] is not None else np.nan,
                'pmra': float(gaia_row['pmra']) if gaia_row['pmra'] is not None else np.nan,
                'pmdec': float(gaia_row['pmdec']) if gaia_row['pmdec'] is not None else np.nan,
            }
        except Exception as e:
            logger.warning(f"Error extracting Gaia data: {str(e)}")
            return {}
    
    def _extract_vizier_catalog_data(self, vizier_row):
        """Extract data from VizieR B catalog row with robust column existence checks."""
        try:
            # Use safer table.colnames access
            try:
                available_cols = vizier_row.table.colnames
            except AttributeError:
                # Fallback for different astropy versions
                available_cols = vizier_row.colnames
            
            extracted_data = {}
            
            # Helper function to safely extract field with fallback
            def safe_extract(field_names, default_value, converter=None):
                """Safely extract a field with multiple possible column names."""
                if isinstance(field_names, str):
                    field_names = [field_names]
                
                for field_name in field_names:
                    if field_name in available_cols:
                        try:
                            value = vizier_row[field_name]
                            if value is not None and value != '':
                                if converter:
                                    return converter(value)
                                return value
                        except (ValueError, TypeError) as e:
                            logger.debug(f"Failed to convert VizieR field {field_name}: {e}")
                            continue
                return default_value
            
            # Extract main identifier (different catalogs use different column names)
            extracted_data['main_id'] = safe_extract([
                'Name', 'MAIN_ID', 'HIP', 'HD', 'HR', 'Source', 'recno', 'source_id'
            ], '', str).strip()
            
            # Extract object type - VizieR catalogs use different naming
            extracted_data['otype'] = safe_extract([
                'VarType', 'Type', 'Class', 'OType', 'OTYPE'
            ], '', str).strip()
            
            # Extract spectral type from various VizieR catalogs
            extracted_data['sp_type'] = safe_extract([
                'SpType', 'SpT', 'MK', 'Spectrum', 'SP_TYPE', 'SPTYPE'
            ], '', str).strip()
            
            # Extract distance - handle parallax conversion if needed
            distance_result = safe_extract(['Dist', 'Distance', 'distance'], np.nan, float)
            if np.isnan(distance_result):
                # Try parallax conversion
                plx = safe_extract(['Plx', 'parallax'], np.nan, float)
                if not np.isnan(plx) and plx > 0:
                    distance_result = 1000.0 / plx  # Convert mas to pc
            extracted_data['distance_result'] = distance_result
            
            # Extract radial velocity from VizieR catalogs
            extracted_data['rv_value'] = safe_extract([
                'RV', 'Vrad', 'RadVel', 'RVZ_RADVEL', 'rv_value'
            ], np.nan, float)
            
            # Get catalog source information
            extracted_data['catalog_source'] = safe_extract(['catalog_source'], '', str)
            
            logger.debug(f"Extracted VizieR data: main_id='{extracted_data.get('main_id', '')}', "
                        f"otype='{extracted_data.get('otype', '')}', "
                        f"catalog='{extracted_data.get('catalog_source', '')}', "
                        f"available_cols={len(available_cols)}")
            return extracted_data
            
        except Exception as e:
            logger.warning(f"Error extracting VizieR catalog data: {str(e)}")
            logger.debug(f"Available VizieR columns: {getattr(vizier_row, 'colnames', 'unknown')}")
            # Return default structure on error
            return {
                'main_id': '',
                'otype': '',
                'sp_type': '',
                'distance_result': np.nan,
                'rv_value': np.nan
            }