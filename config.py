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
        ("TIFF files", "*.tiff"),
        ("JPEG files", "*.jpg"),
        ("GIF files", "*.gif"),
        ("All files", "*.*")
    ],
    'fits': [
        ("FITS files", "*.fits"),
        ("All files", "*.*")
    ],
    'raw': [
        ("CR2 files", "*.cr2"),
        ("NEF files", "*.nef"),
        ("ARW files", "*.arw"),
        ("DNG files", "*.dng"),
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