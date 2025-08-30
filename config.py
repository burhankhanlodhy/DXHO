"""
Configuration module for the Digital Exoplanet Hunting Observatory.
"""

import logging
import os
from datetime import datetime


def setup_logging():
    """
    Set up logging configuration for the application with UTF-8 encoding support.
    """
    # Create logs directory if it doesn't exist
    log_dir = "logs"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    # Create log filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = os.path.join(log_dir, f"deho_{timestamp}.log")
    
    # Configure logging with UTF-8 encoding support
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename, encoding='utf-8'),
            logging.StreamHandler()  # Console handler (uses system default encoding)
        ]
    )
    
    logger = logging.getLogger(__name__)
    logger.info("Logging initialized with UTF-8 encoding support")
    logger.info(f"Log file: {log_filename}")


# Application configuration
APP_CONFIG = {
    'name': 'Digital Exoplanet Hunting Observatory',
    'version': '2.0',
    'window_size': '720x480',
    'window_position': '+270+150',
    'about_text': {
        'name': 'Digital Exoplanet Hunting Observatory',
        'version': '2.0',
        'developers': [
            'M. Burhan Khan Lodhi',
            'Shomyal Khan', 
            'Mushahid Zafar Jafri'
        ],
        'institution': 'University of Management and Technology, Sialkot'
    }
}

# File type configurations
FILE_TYPES = {
    'images': [
        ("PNG files", "*.png"),
        ("TIFF files", "*.tiff;*.tif"),
        ("JPEG files", "*.jpg;*.jpeg"),
        ("BMP files", "*.bmp"),
        ("GIF files", "*.gif"),
        ("All supported images", "*.png;*.tiff;*.tif;*.jpg;*.jpeg;*.bmp;*.gif"),
        ("All files", "*.*")
    ],
    'fits': [
        ("FITS files", "*.fits"),
        ("All files", "*.*")
    ],
    'raw': [
        ("Canon RAW (CR2/CR3)", "*.cr2;*.cr3"),
        ("Nikon RAW (NEF)", "*.nef"),
        ("Sony RAW (ARW)", "*.arw"),
        ("Adobe DNG", "*.dng"),
        ("Fujifilm RAW (RAF)", "*.raf"),
        ("Olympus RAW (ORF)", "*.orf"),
        ("Panasonic RAW (RW2)", "*.rw2"),
        ("Pentax RAW (PEF)", "*.pef"),
        ("Samsung RAW (SRW)", "*.srw"),
        ("All RAW formats", "*.cr2;*.cr3;*.nef;*.arw;*.dng;*.raf;*.orf;*.rw2;*.pef;*.srw"),
        ("All files", "*.*")
    ],
    'csv': [
        ("CSV files", "*.csv"),
        ("All files", "*.*")
    ]
}

# Default analysis parameters
DEFAULT_PARAMS = {
    'fwhm': 5.0,
    'threshold': 3.0,
    'aperture_radius': 17.0
}