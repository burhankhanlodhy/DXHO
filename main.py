#!/usr/bin/env python3
"""
Digital Exoplanet Hunting Observatory
Main application entry point.

This is the main entry point for the astronomical image processing application.
It initializes logging, creates the GUI, and starts the application.

Authors:
    - M. Burhan Khan Lodhi
    - Shomyal Khan
    - Mushahid Zafar Jafri

Institution: University of Management and Technology, Sialkot
"""

import sys
import os
import tkinter as tk
from tkinter import messagebox
import logging

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from config import setup_logging, APP_CONFIG
    from gui import MainGUI
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure all required modules are in the same directory.")
    sys.exit(1)


def check_dependencies():
    """
    Check if all required dependencies are available.
    
    Returns:
        bool: True if all dependencies are available, False otherwise
    """
    required_modules = [
        'numpy', 'pandas', 'astropy', 'photutils', 
        'matplotlib', 'rawpy', 'exifread', 'PIL'
    ]
    
    missing_modules = []
    
    for module in required_modules:
        try:
            __import__(module)
        except ImportError:
            missing_modules.append(module)
    
    if missing_modules:
        error_msg = (
            f"Missing required modules: {', '.join(missing_modules)}\\n\\n"
            "Please install them using:\\n"
            "pip install -r requirements.txt"
        )
        
        # Try to show GUI error if tkinter is available
        try:
            root = tk.Tk()
            root.withdraw()  # Hide main window
            messagebox.showerror("Missing Dependencies", error_msg)
            root.destroy()
        except:
            print(error_msg)
        
        return False
    
    return True


def main():
    """
    Main application entry point.
    """
    try:
        # Check dependencies first
        if not check_dependencies():
            return 1
        
        # Set up logging
        setup_logging()
        logger = logging.getLogger(__name__)
        
        logger.info("="*50)
        logger.info(f"Starting {APP_CONFIG['name']} v{APP_CONFIG['version']}")
        logger.info("="*50)
        
        # Create and run the GUI application
        app = MainGUI()
        app.run()
        
        logger.info("Application terminated normally")
        return 0
        
    except Exception as e:
        # Log the error if logging is set up
        try:
            logger = logging.getLogger(__name__)
            logger.critical(f"Critical error in main: {str(e)}", exc_info=True)
        except:
            print(f"Critical error: {e}")
        
        # Try to show GUI error
        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("Critical Error", f"A critical error occurred:\\n\\n{str(e)}")
            root.destroy()
        except:
            pass
        
        return 1
    
    except KeyboardInterrupt:
        try:
            logger = logging.getLogger(__name__)
            logger.info("Application interrupted by user")
        except:
            print("Application interrupted by user")
        return 0


if __name__ == "__main__":
    sys.exit(main())