#!/usr/bin/env python3
"""
Digital Exoplanet Hunting Observatory v2.0
Improved version with better error handling, logging, and user experience.

Authors:
- M. Burhan Khan Lodhi
- Shomyal Khan
- Mushahid Zafar Jafri

Institution: University of Management and Technology, Sialkot
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std
import matplotlib.pyplot as plt
from photutils.aperture import aperture_photometry, CircularAperture
import rawpy
import exifread
from PIL import Image
import tkinter as tk
from tkinter import filedialog, messagebox
import tkinter.ttk as tkk
import logging
import threading
import os
from datetime import datetime

# Set up logging
def setup_logging():
    log_dir = "logs"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = os.path.join(log_dir, f"deho_{timestamp}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

logger = setup_logging()

# Application configuration
APP_CONFIG = {
    'name': 'Digital Exoplanet Hunting Observatory',
    'version': '2.0',
    'developers': ['M. Burhan Khan Lodhi', 'Shomyal Khan', 'Mushahid Zafar Jafri'],
    'institution': 'University of Management and Technology, Sialkot'
}

FILE_TYPES = {
    'images': [("PNG files", "*.png"), ("TIFF files", "*.tiff"), ("JPEG files", "*.jpg"), ("GIF files", "*.gif"), ("All files", "*.*")],
    'fits': [("FITS files", "*.fits"), ("All files", "*.*")],
    'raw': [("CR2 files", "*.cr2"), ("NEF files", "*.nef"), ("ARW files", "*.arw"), ("DNG files", "*.dng"), ("All files", "*.*")],
    'csv': [("CSV files", "*.csv"), ("All files", "*.*")]
}

class GUI:
    """Improved GUI class with better error handling and user experience."""
    
    def __init__(self):
        """Initialize the GUI with improved error handling."""
        self.root = None
        self.statusbar = None
        self.progressbar = None
        self.progressbar2 = None
        self.progress_window = None
        self.photometry_window = None
        
        # File paths
        self.current_image_path = None
        self.current_fits_path = None
        self.current_raw_path = None
        self.fits_files = []
        
        # Menu references for enabling/disabling
        self.file_menu = None
        self.process_menu = None
        self.analysis_menu = None
        
        # Button references
        self.convert_btn = None
        self.stack_btn = None
        self.photometry_btn = None
        
        logger.info('GUI initialized')

    def open(self):
        """Open an image file with improved error handling."""
        try:
            filepath = filedialog.askopenfilename(
                title="Select Image File",
                filetypes=FILE_TYPES['images']
            )
            if filepath:
                self.current_image_path = filepath
                self._update_status(f"Image loaded: {os.path.basename(filepath)}")
                logger.info(f"Image loaded: {filepath}")
                self._update_ui_state()
                self._add_info_message(f"Image file loaded: {os.path.basename(filepath)}")
        except Exception as e:
            logger.error(f"Error opening image: {str(e)}")
            messagebox.showerror('Error', f'Failed to load image: {e}')

    def open_fits(self):
        """Open multiple FITS files for stacking."""
        try:
            filepaths = filedialog.askopenfilenames(
                title="Select FITS Files",
                filetypes=FILE_TYPES['fits']
            )
            if filepaths:
                self.fits_files = filepaths
                self._update_status(f"Selected {len(filepaths)} FITS files")
                logger.info(f"Selected {len(filepaths)} FITS files for stacking")
                self._add_info_message(f"Selected {len(filepaths)} FITS files for stacking")
        except Exception as e:
            logger.error(f"Error selecting FITS files: {str(e)}")
            messagebox.showerror('Error', f'Failed to select FITS files: {e}')

    def raw(self):
        """Open a RAW file with improved file type support."""
        try:
            filepath = filedialog.askopenfilename(
                title="Select RAW File",
                filetypes=FILE_TYPES['raw']
            )
            if filepath:
                self.current_raw_path = filepath
                self._update_status(f"RAW file ready: {os.path.basename(filepath)}")
                logger.info(f"RAW file selected: {filepath}")
                self._update_ui_state()
                self._add_info_message(f"RAW file loaded: {os.path.basename(filepath)}")
        except Exception as e:
            logger.error(f"Error selecting RAW file: {str(e)}")
            messagebox.showerror('Error', f'Failed to select RAW file: {e}')

    def fit(self):
        """Open a single FITS file for analysis."""
        try:
            filepath = filedialog.askopenfilename(
                title="Select FITS File",
                filetypes=FILE_TYPES['fits']
            )
            if filepath:
                self.current_fits_path = filepath
                self._update_status(f"FITS file ready: {os.path.basename(filepath)}")
                logger.info(f"FITS file selected: {filepath}")
                self._update_ui_state()
                self._add_info_message(f"FITS file loaded: {os.path.basename(filepath)}")
        except Exception as e:
            logger.error(f"Error selecting FITS file: {str(e)}")
            messagebox.showerror('Error', f'Failed to select FITS file: {e}')

    def save(self):
        """Get save location for converted files."""
        try:
            self.savepath = filedialog.asksaveasfilename(
                title="Save FITS File",
                defaultextension=".fits",
                filetypes=FILE_TYPES['fits']
            )
            if self.savepath:
                logger.info(f"Save path selected: {self.savepath}")
        except Exception as e:
            logger.error(f"Error selecting save path: {str(e)}")
            messagebox.showerror('Error', f'Failed to select save location: {e}')

    def fits(self):
        """Select FITS file for photometry analysis."""
        self.fit()

    def fitsopen(self):
        """Load and display FITS file."""
        if not self.current_fits_path:
            messagebox.showwarning('No FITS File', 'Please select a FITS file first.')
            return
        
        try:
            F1.load_fits(self.current_fits_path)
            self._update_status("FITS file displayed")
        except Exception as e:
            logger.error(f"Error loading FITS: {str(e)}")
            messagebox.showerror('Error', f'Failed to load FITS file: {e}')

    def conversion(self):
        """Start RAW to FITS conversion process."""
        if not self.current_raw_path:
            messagebox.showwarning('No RAW File', 'Please select a RAW file first.')
            self.raw()
            return
        
        if not hasattr(self, 'savepath') or not self.savepath:
            self.save()
            if not hasattr(self, 'savepath') or not self.savepath:
                return
        
        self._update_status("Converting RAW to FITS...")
        self.progress()

    def photometry(self):
        """Start photometry analysis with improved dialog."""
        if not self.current_fits_path:
            messagebox.showwarning('No FITS File', 'Please select a FITS file first.')
            self.fit()
            if not self.current_fits_path:
                return
        
        self._show_photometry_dialog()

    def _show_photometry_dialog(self):
        """Show improved photometry parameters dialog."""
        self.photometry_window = tk.Toplevel(self.root)
        self.photometry_window.title("Aperture Photometry Parameters")
        self.photometry_window.geometry("380x250")
        self.photometry_window.resizable(False, False)
        
        # Make window modal
        self.photometry_window.transient(self.root)
        self.photometry_window.grab_set()
        
        # Center the window
        self.photometry_window.update_idletasks()
        x = (self.photometry_window.winfo_screenwidth() // 2) - (self.photometry_window.winfo_width() // 2)
        y = (self.photometry_window.winfo_screenheight() // 2) - (self.photometry_window.winfo_height() // 2)
        self.photometry_window.geometry(f"+{x}+{y}")
        
        # Main frame
        main_frame = tkk.Frame(self.photometry_window, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Instructions
        tkk.Label(main_frame, text="Configure photometry parameters:", font=('Arial', 10, 'bold')).pack(pady=(0, 15))
        
        # FWHM parameter
        fwhm_frame = tkk.Frame(main_frame)
        fwhm_frame.pack(fill=tk.X, pady=5)
        tkk.Label(fwhm_frame, text="FWHM (pixels):").pack(side=tk.LEFT)
        self.fwhm_var = tk.StringVar(value="5.0")
        fwhm_entry = tkk.Entry(fwhm_frame, textvariable=self.fwhm_var, width=15)
        fwhm_entry.pack(side=tk.RIGHT)
        
        # Threshold parameter
        threshold_frame = tkk.Frame(main_frame)
        threshold_frame.pack(fill=tk.X, pady=5)
        tkk.Label(threshold_frame, text="Threshold (sigma):").pack(side=tk.LEFT)
        self.threshold_var = tk.StringVar(value="3.0")
        threshold_entry = tkk.Entry(threshold_frame, textvariable=self.threshold_var, width=15)
        threshold_entry.pack(side=tk.RIGHT)
        
        # Help text
        help_text = "FWHM: Typical values 3-8 pixels\\nThreshold: Typical values 2-5 sigma"
        tkk.Label(main_frame, text=help_text, font=('Arial', 8), foreground='gray').pack(pady=10)
        
        # Buttons
        button_frame = tkk.Frame(main_frame)
        button_frame.pack(pady=15)
        tkk.Button(button_frame, text="Start Analysis", command=self.perform).pack(side=tk.LEFT, padx=5)
        tkk.Button(button_frame, text="Cancel", command=self.photometry_window.destroy).pack(side=tk.LEFT, padx=5)
        
        # Focus on first entry
        fwhm_entry.focus()

    def perform(self):
        """Perform photometry analysis with validation."""
        try:
            fwhm_value = float(self.fwhm_var.get())
            threshold_value = float(self.threshold_var.get())
            
            # Validate parameters
            if fwhm_value <= 0 or threshold_value <= 0:
                messagebox.showerror('Invalid Parameters', 'FWHM and threshold must be positive values.')
                return
            
            if fwhm_value < 1:
                messagebox.showwarning('Low FWHM', 'FWHM less than 1 pixel may not detect stars properly.')
                if not messagebox.askyesno('Continue?', 'Continue with FWHM < 1?'):
                    return
            
            if fwhm_value > 50:
                messagebox.showwarning('High FWHM', 'FWHM greater than 50 pixels is unusually large.')
                if not messagebox.askyesno('Continue?', 'Continue with large FWHM?'):
                    return
            
            if threshold_value > 20:
                messagebox.showwarning('High Threshold', 'Threshold greater than 20 sigma is very high and may miss faint stars.')
                if not messagebox.askyesno('Continue?', 'Continue with high threshold?'):
                    return
            
            self.photometry_window.destroy()
            self._update_status("Performing photometry analysis...")
            self.photometry_progress()
            
            # Run photometry in separate thread to prevent GUI freezing
            def photometry_thread():
                try:
                    result = P1.photometry(self.current_fits_path, fwhm_value, threshold_value)
                    if result:
                        self._update_status("Photometry analysis completed")
                        logger.info("Photometry analysis completed successfully")
                except Exception as e:
                    logger.error(f"Error in photometry thread: {str(e)}")
                    messagebox.showerror('Analysis Error', f'Photometry failed: {e}')
            
            threading.Thread(target=photometry_thread, daemon=True).start()
            
        except ValueError:
            messagebox.showerror('Invalid Input', 'Please enter valid numeric values.')
        except Exception as e:
            logger.error(f"Error in photometry setup: {str(e)}")
            messagebox.showerror('Error', f'Failed to start photometry: {e}')

    def stack(self):
        """Stack FITS images with improved error handling."""
        if not self.fits_files:
            self.open_fits()
            if not self.fits_files:
                messagebox.showinfo('Stack Images', 'Please select multiple FITS files to stack.')
                return
        
        if len(self.fits_files) < 2:
            messagebox.showwarning('Insufficient Files', 'Please select at least 2 FITS files for stacking.')
            self.open_fits()
            if len(self.fits_files) < 2:
                return
        
        try:
            self._update_status(f"Stacking {len(self.fits_files)} images...")
            data = F1.mean_fits(self.fits_files)
            if data is not None:
                # Display the stacked result
                plt.figure(figsize=(12, 10))
                vmin, vmax = np.percentile(data, [1, 99])
                plt.imshow(data.T, cmap='viridis', vmin=vmin, vmax=vmax, origin='lower')
                plt.colorbar(label='Mean Counts')
                plt.title(f'Stacked Image ({len(self.fits_files)} frames)')
                plt.xlabel('X (pixels)')
                plt.ylabel('Y (pixels)')
                plt.tight_layout()
                plt.show()
                
                self._update_status("Image stacking completed")
                logger.info(f"Successfully stacked {len(self.fits_files)} images")
                
                # Ask if user wants to save
                if messagebox.askyesno("Save Stacked Image", "Would you like to save the stacked image?"):
                    save_path = filedialog.asksaveasfilename(
                        title="Save Stacked FITS",
                        defaultextension=".fits",
                        filetypes=FILE_TYPES['fits']
                    )
                    if save_path:
                        hdu = fits.PrimaryHDU(data=data)
                        hdu.writeto(save_path, overwrite=True)
                        messagebox.showinfo("Saved", f"Stacked image saved to {save_path}")
        except Exception as e:
            logger.error(f"Error in image stacking: {str(e)}")
            messagebox.showerror('Stacking Error', f'Failed to stack images: {e}')

    def about(self):
        """Show improved about dialog."""
        about_window = tk.Toplevel(self.root)
        about_window.title("About")
        about_window.geometry("420x320")
        about_window.resizable(False, False)
        
        # Center the window
        about_window.update_idletasks()
        x = (about_window.winfo_screenwidth() // 2) - (about_window.winfo_width() // 2)
        y = (about_window.winfo_screenheight() // 2) - (about_window.winfo_height() // 2)
        about_window.geometry(f"+{x}+{y}")
        
        # Content frame
        content_frame = tkk.Frame(about_window, padding="20")
        content_frame.pack(fill=tk.BOTH, expand=True)
        
        # Application name and version
        tkk.Label(content_frame, text=APP_CONFIG['name'], 
                 font=('Arial', 14, 'bold')).pack(pady=(0, 5))
        tkk.Label(content_frame, text=f"Version {APP_CONFIG['version']}", 
                 font=('Arial', 10)).pack(pady=(0, 20))
        
        # Developers
        tkk.Label(content_frame, text="Developed by:", 
                 font=('Arial', 11, 'bold')).pack(pady=(0, 10))
        for dev in APP_CONFIG['developers']:
            tkk.Label(content_frame, text=f"• {dev}", 
                     font=('Arial', 9)).pack()
        
        # Institution
        tkk.Label(content_frame, text=APP_CONFIG['institution'], 
                 font=('Arial', 9), foreground='gray').pack(pady=20)
        
        # Features
        features_text = ("Features: RAW Processing • Image Stacking • Aperture Photometry\\n"
                        "Supports: CR2, NEF, ARW, DNG, FITS formats")
        tkk.Label(content_frame, text=features_text, 
                 font=('Arial', 8), foreground='navy').pack(pady=10)
        
        # Close button
        tkk.Button(content_frame, text="Close", 
                  command=about_window.destroy).pack(pady=15)

    def help(self):
        """Show comprehensive help dialog."""
        help_window = tk.Toplevel(self.root)
        help_window.title("User Guide")
        help_window.geometry("650x550")
        
        # Create text widget with scrollbar
        text_frame = tkk.Frame(help_window)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        help_text = tk.Text(text_frame, wrap=tk.WORD, font=('Consolas', 9))
        scrollbar = tkk.Scrollbar(text_frame, orient=tk.VERTICAL, command=help_text.yview)
        help_text.configure(yscrollcommand=scrollbar.set)
        
        help_content = f"""
{APP_CONFIG['name']} v{APP_CONFIG['version']} - User Guide

═══════════════════════════════════════════════════════════════

1. OPENING FILES
   • File → Open Image: Load standard image formats (PNG, TIFF, JPEG)
   • File → Open Raw File: Load astronomical RAW images
   • File → Load Fits: Load FITS astronomical images
   
   Supported RAW formats: CR2, NEF, ARW, DNG
   Supported standard formats: PNG, TIFF, JPEG, GIF

2. RAW TO FITS CONVERSION
   • Open a RAW file first
   • Select Process → Convert to FITS
   • Choose save location
   • EXIF metadata is preserved in FITS header
   
   Benefits:
   - Preserves full dynamic range
   - Maintains astronomical metadata
   - Compatible with analysis software

3. IMAGE STACKING
   • Select Process → Stack Images
   • Choose multiple FITS files of the same target
   • Images are mean-combined to reduce noise
   • Optionally save the stacked result
   
   Tips:
   - Use images of the same exposure settings
   - Ensure proper alignment before stacking
   - More frames = better noise reduction

4. APERTURE PHOTOMETRY
   • Load a FITS file for analysis
   • Select Analysis → Aperture Photometry
   • Configure parameters:
     * FWHM: Stellar profile width (typical: 3-8 pixels)
     * Threshold: Detection sensitivity (typical: 2-5 sigma)
   • Results displayed and optionally saved as single CSV file
   
   Output file contains:
   - Source positions (X, Y coordinates)
   - Aperture photometry measurements
   - Background statistics
   - Combined in one convenient CSV file

5. TIPS FOR BEST RESULTS
   • Use dark frames to reduce noise
   • Take multiple exposures for stacking
   • Adjust photometry parameters for your conditions
   • Check log files for processing details
   • Save work in FITS format for analysis

6. TROUBLESHOOTING
   • Check log files in logs/ directory
   • Ensure sufficient disk space
   • Verify file permissions
   • Close other memory-intensive applications

For technical support, contact the development team.
        """
        
        help_text.insert(tk.END, help_content)
        help_text.config(state=tk.DISABLED)
        
        help_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    def show_image(self):
        """Display current image or RAW file with improved viewer."""
        # Check if we have either a regular image or RAW file
        filepath = self.current_image_path or self.current_raw_path
        if not filepath:
            messagebox.showwarning('No Image', 'No image or RAW file is currently loaded.\\n\\nPlease use File > Open Image or File > Open Raw File first.')
            return
        
        try:
            img_window = tk.Toplevel(self.root)
            img_window.title(f"Image Viewer: {os.path.basename(filepath)}")
            
            # Create scrollable canvas
            canvas = tk.Canvas(img_window, bg='black')
            v_scrollbar = tkk.Scrollbar(img_window, orient=tk.VERTICAL, command=canvas.yview)
            h_scrollbar = tkk.Scrollbar(img_window, orient=tk.HORIZONTAL, command=canvas.xview)
            canvas.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
            
            # Load and display image (handle RAW files)
            if filepath.lower().endswith(('.cr2', '.nef', '.arw', '.dng')):
                # Handle RAW files
                with rawpy.imread(filepath) as raw:
                    rgb = raw.postprocess(no_auto_bright=True, use_auto_wb=False, gamma=None)
                pil_image = Image.fromarray(rgb)
            else:
                # Handle regular image files
                pil_image = Image.open(filepath)
            # Resize if too large
            max_size = (1200, 900)
            original_size = pil_image.size
            pil_image.thumbnail(max_size, Image.Resampling.LANCZOS)
            
            # Convert to PhotoImage - handle different formats
            try:
                # For processed images (especially RAW), use PIL conversion
                from PIL import ImageTk
                photo = ImageTk.PhotoImage(pil_image, master=canvas)
            except Exception:
                # Fallback: try direct loading for regular image files
                if not filepath.lower().endswith(('.cr2', '.nef', '.arw', '.dng')):
                    photo = tk.PhotoImage(master=canvas, file=filepath)
                else:
                    raise
            canvas.create_image(0, 0, image=photo, anchor=tk.NW)
            canvas.configure(scrollregion=canvas.bbox("all"))
            
            # Pack elements
            canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
            h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
            
            # Status bar for image info
            info_frame = tkk.Frame(img_window)
            info_frame.pack(side=tk.BOTTOM, fill=tk.X)
            info_text = f"Original: {original_size[0]}x{original_size[1]} | Displayed: {pil_image.size[0]}x{pil_image.size[1]}"
            tkk.Label(info_frame, text=info_text).pack(pady=2)
            
            # Keep reference to prevent garbage collection
            canvas.image = photo
            
            img_window.geometry("900x700")
            
        except Exception as e:
            logger.error(f"Error displaying image: {str(e)}")
            messagebox.showerror('Display Error', f'Failed to display image: {e}')

    def progress(self):
        """Show improved progress window for conversion."""
        self.progress_window = tk.Toplevel(self.root)
        self.progress_window.title("RAW to FITS Conversion")
        self.progress_window.geometry("450x120")
        self.progress_window.resizable(False, False)
        
        # Make window modal
        self.progress_window.transient(self.root)
        self.progress_window.grab_set()
        
        # Center the window
        self.progress_window.update_idletasks()
        x = (self.progress_window.winfo_screenwidth() // 2) - (self.progress_window.winfo_width() // 2)
        y = (self.progress_window.winfo_screenheight() // 2) - (self.progress_window.winfo_height() // 2)
        self.progress_window.geometry(f"+{x}+{y}")
        
        # Content
        content_frame = tkk.Frame(self.progress_window, padding="20")
        content_frame.pack(fill=tk.BOTH, expand=True)
        
        tkk.Label(content_frame, text="Converting RAW to FITS format...", 
                 font=('Arial', 10)).pack(pady=(0, 10))
        
        self.progressbar = tkk.Progressbar(content_frame, mode='determinate', length=400)
        self.progressbar.pack(pady=10)
        
        # Start conversion in separate thread
        threading.Thread(target=self.conv, daemon=True).start()

    def photometry_progress(self):
        """Show progress window for photometry analysis."""
        self.phot_progress_window = tk.Toplevel(self.root)
        self.phot_progress_window.title("Photometry Analysis")
        self.phot_progress_window.geometry("450x120")
        self.phot_progress_window.resizable(False, False)
        
        # Make window modal
        self.phot_progress_window.transient(self.root)
        self.phot_progress_window.grab_set()
        
        # Center window
        self.phot_progress_window.update_idletasks()
        x = (self.phot_progress_window.winfo_screenwidth() // 2) - (self.phot_progress_window.winfo_width() // 2)
        y = (self.phot_progress_window.winfo_screenheight() // 2) - (self.phot_progress_window.winfo_height() // 2)
        self.phot_progress_window.geometry(f"+{x}+{y}")
        
        # Content
        content_frame = tkk.Frame(self.phot_progress_window, padding="20")
        content_frame.pack(fill=tk.BOTH, expand=True)
        
        tkk.Label(content_frame, text="Performing aperture photometry...", 
                 font=('Arial', 10)).pack(pady=(0, 10))
        
        self.progressbar2 = tkk.Progressbar(content_frame, mode='determinate', length=400)
        self.progressbar2.pack(pady=10)
        
        # Store reference for cleanup
        self.phot_barwindow = self.phot_progress_window

    def conv(self):
        """Perform RAW to FITS conversion."""
        try:
            success = C1.raw_to_fits(self.current_raw_path, self.savepath)
            if success:
                self.call_me()  # Show completion message
        except Exception as e:
            logger.error(f"Conversion error: {str(e)}")
            messagebox.showerror('Conversion Error', f'Failed to convert: {e}')
        finally:
            # Always close progress window
            if hasattr(self, 'progress_window') and self.progress_window:
                self.progress_window.destroy()
                self.progress_window = None

    def call_me(self):
        """Show conversion complete message."""
        self._update_status("RAW to FITS conversion completed")
        if hasattr(self, 'savepath') and self.savepath:
            output_file = self.savepath if self.savepath.endswith('.fits') else f"{self.savepath}.fits"
            messagebox.showinfo("Conversion Complete", 
                f"RAW file successfully converted to FITS format!\\n\\nSaved to: {os.path.basename(output_file)}")
        else:
            messagebox.showinfo("Conversion Complete", "RAW file successfully converted to FITS format!")

    def _update_status(self, message):
        """Update status bar with message."""
        if self.statusbar:
            self.statusbar.config(text=message)
            self.root.update_idletasks()
        logger.info(message)

    def mainwindow(self):
        """Create and configure the main application window."""
        self.root = tk.Tk()
        self.root.geometry("800x600+200+100")
        self.root.title(f"{APP_CONFIG['name']} v{APP_CONFIG['version']}")
        self.root.minsize(600, 400)
        
        # Configure style
        style = tkk.Style()
        style.theme_use('clam')
        
        # Create menu bar
        main_menu = tk.Menu(self.root)
        self.root.config(menu=main_menu)

        # File Menu
        self.file_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="File", menu=self.file_menu)
        self.file_menu.add_command(label="Open Image", command=self.open, accelerator="Ctrl+O")
        self.file_menu.add_command(label="Open Raw File", command=self.raw)
        self.file_menu.add_command(label="Load Fits", command=self.fit)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Show Image", command=self.show_image, accelerator="Ctrl+I", state=tk.DISABLED)
        self.file_menu.add_command(label="Show Fits", command=self.fitsopen, state=tk.DISABLED)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command=self._safe_exit, accelerator="Ctrl+Q")

        # Process Menu
        self.process_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="Process", menu=self.process_menu)
        self.process_menu.add_command(label="Convert to FITS", command=self.conversion, state=tk.DISABLED)
        self.process_menu.add_separator()
        self.process_menu.add_command(label="Stack Images", command=self.stack)

        # Analysis Menu
        self.analysis_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="Analysis", menu=self.analysis_menu)
        self.analysis_menu.add_command(label="Aperture Photometry", command=self.photometry, state=tk.DISABLED)

        # Settings Menu
        self.settings_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="Settings", menu=self.settings_menu)
        self.settings_menu.add_command(label="WSL Astrometry Settings", command=self.open_wsl_settings)

        # Help Menu
        help_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="User Guide", command=self.help, accelerator="F1")
        help_menu.add_separator()
        help_menu.add_command(label="About", command=self.about)

        # Create main content area
        self._create_main_content()

        # Create status bar
        self.statusbar = tkk.Label(self.root, text="Ready", relief=tk.SUNKEN, 
                                  anchor=tk.W, padding=5)
        self.statusbar.pack(side=tk.BOTTOM, fill=tk.X)

        # Keyboard shortcuts
        self.root.bind('<Control-o>', lambda e: self.open())
        self.root.bind('<Control-i>', lambda e: self.show_image())
        self.root.bind('<Control-q>', lambda e: self._safe_exit())
        self.root.bind('<F1>', lambda e: self.help())
        
        # Handle window closing
        self.root.protocol("WM_DELETE_WINDOW", self._safe_exit)
        
        # Set initial focus
        self.root.focus_force()
        
        logger.info("Main window created and configured")

    def _create_main_content(self):
        """Create the main content area."""
        # Main frame
        main_frame = tkk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Welcome section
        welcome_frame = tkk.LabelFrame(main_frame, text="Welcome", padding="10")
        welcome_frame.pack(fill=tk.X, pady=(0, 10))
        
        tkk.Label(welcome_frame, text=APP_CONFIG['name'], 
                 font=('Arial', 14, 'bold')).pack()
        tkk.Label(welcome_frame, text="Astronomical Image Processing and Analysis", 
                 font=('Arial', 10)).pack()
        
        # Quick actions
        actions_frame = tkk.LabelFrame(main_frame, text="Quick Actions", padding="10")
        actions_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Button grid with state management
        # Row 1: File operations
        tkk.Button(actions_frame, text="Open Image", command=self.open, width=15).grid(row=0, column=0, padx=5, pady=5, sticky='ew')
        tkk.Button(actions_frame, text="Open RAW", command=self.raw, width=15).grid(row=0, column=1, padx=5, pady=5, sticky='ew')
        tkk.Button(actions_frame, text="Open FITS", command=self.fit, width=15).grid(row=0, column=2, padx=5, pady=5, sticky='ew')
        
        # Row 2: Processing operations (initially disabled)
        self.convert_btn = tkk.Button(actions_frame, text="Convert RAW", command=self.conversion, width=15, state=tk.DISABLED)
        self.convert_btn.grid(row=1, column=0, padx=5, pady=5, sticky='ew')
        
        self.stack_btn = tkk.Button(actions_frame, text="Stack Images", command=self.stack, width=15)
        self.stack_btn.grid(row=1, column=1, padx=5, pady=5, sticky='ew')
        
        self.photometry_btn = tkk.Button(actions_frame, text="Photometry", command=self.photometry, width=15, state=tk.DISABLED)
        self.photometry_btn.grid(row=1, column=2, padx=5, pady=5, sticky='ew')
        
        # Configure grid weights
        for i in range(3):
            actions_frame.columnconfigure(i, weight=1)
        
        # Information panel
        info_frame = tkk.LabelFrame(main_frame, text="Session Information", padding="5")
        info_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create text widget for session info
        self.info_text = tk.Text(info_frame, height=8, state=tk.DISABLED, wrap=tk.WORD,
                                font=('Consolas', 9))
        info_scrollbar = tkk.Scrollbar(info_frame, orient=tk.VERTICAL, command=self.info_text.yview)
        self.info_text.configure(yscrollcommand=info_scrollbar.set)
        
        self.info_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        info_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Add welcome message
        self._add_info_message("Application started successfully")
        self._add_info_message(f"Version {APP_CONFIG['version']} ready for astronomical processing")
        
        # Update initial UI state
        self._update_ui_state()

    def _add_info_message(self, message):
        """Add a message to the information panel."""
        timestamp = datetime.now().strftime("%H:%M:%S")
        formatted_message = f"[{timestamp}] {message}\\n"
        
        self.info_text.config(state=tk.NORMAL)
        self.info_text.insert(tk.END, formatted_message)
        self.info_text.see(tk.END)
        self.info_text.config(state=tk.DISABLED)

    def _update_ui_state(self):
        """Update UI state based on loaded files."""
        # Update menu states
        if self.file_menu:
            # Show Image menu item (enable for both regular images and RAW files)
            if self.current_image_path or self.current_raw_path:
                self.file_menu.entryconfig("Show Image", state=tk.NORMAL)
            else:
                self.file_menu.entryconfig("Show Image", state=tk.DISABLED)
            
            # Show FITS menu item
            if self.current_fits_path:
                self.file_menu.entryconfig("Show Fits", state=tk.NORMAL)
            else:
                self.file_menu.entryconfig("Show Fits", state=tk.DISABLED)
        
        # Update process menu states
        if self.process_menu:
            # Convert to FITS menu item
            if self.current_raw_path:
                self.process_menu.entryconfig("Convert to FITS", state=tk.NORMAL)
            else:
                self.process_menu.entryconfig("Convert to FITS", state=tk.DISABLED)
        
        # Update analysis menu states
        if self.analysis_menu:
            # Aperture Photometry menu item
            if self.current_fits_path:
                self.analysis_menu.entryconfig("Aperture Photometry", state=tk.NORMAL)
            else:
                self.analysis_menu.entryconfig("Aperture Photometry", state=tk.DISABLED)
        
        # Update button states
        if hasattr(self, 'convert_btn') and self.convert_btn:
            if self.current_raw_path:
                self.convert_btn.config(state=tk.NORMAL)
            else:
                self.convert_btn.config(state=tk.DISABLED)
        
        if hasattr(self, 'photometry_btn') and self.photometry_btn:
            if self.current_fits_path:
                self.photometry_btn.config(state=tk.NORMAL)
            else:
                self.photometry_btn.config(state=tk.DISABLED)
    
    def _safe_exit(self):
        """Safely exit the application."""
        if messagebox.askyesno("Exit Application", 
                              "Are you sure you want to exit?"):
            logger.info("Application closing by user request")
            self.root.destroy()

    def update_status(self):
        """Update status to ready."""
        self._update_status("Ready")


class Photometry:
    """Enhanced photometry class with better error handling."""
    
    def photometry(self, filepath, fwhm, threshold):
        """Perform aperture photometry with enhanced error handling."""
        try:
            logger.info(f"Starting photometry: {filepath}, FWHM={fwhm}, threshold={threshold}")
            
            # Load FITS file
            hdulist = fits.open(filepath, ignore_missing_end=True)
            hdu = hdulist[0]
            
            G1.progressbar2['value'] = 20
            G1.progressbar2.update()
            
            image = hdu.data.astype(float)
            if image is None or image.size == 0:
                raise ValueError("No image data found in FITS file")
            
            image -= np.median(image)
            bkg_sigma = mad_std(image)
            
            G1.progressbar2['value'] = 40
            G1.progressbar2.update()
            
            # Perform star detection
            daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * bkg_sigma)
            sources = daofind(image)
            
            if sources is None or len(sources) == 0:
                messagebox.showwarning("No Stars Found", 
                    "No stellar sources detected. Try adjusting FWHM or threshold parameters.")
                return False
            
            # Format source table
            for col in sources.colnames:
                sources[col].info.format = '%.8g'
            
            logger.info(f"Found {len(sources)} stellar sources")
            print(f"Detected {len(sources)} sources:")
            print(sources)
            
            G1.progressbar2['value'] = 60
            G1.progressbar2.update()

            # Perform aperture photometry
            positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
            apertures = CircularAperture(positions, r=17.)
            phot_table = aperture_photometry(image, apertures)
            
            for col in phot_table.colnames:
                phot_table[col].info.format = '%.8g'

            G1.progressbar2['value'] = 80
            G1.progressbar2.update()

            G1.progressbar2['value'] = 90
            G1.progressbar2.update()
            
            # Clean up progress window first
            if hasattr(G1, 'phot_barwindow') and G1.phot_barwindow:
                G1.phot_barwindow.destroy()
                G1.phot_barwindow = None
            G1.update_status()
            
            # Display results first
            self._display_photometry_results(image, apertures, len(sources))
            
            # Then save results (user can cancel without affecting display)
            self._save_photometry_results(sources, phot_table)
            
            messagebox.showinfo("Analysis Complete", 
                f"Photometry analysis completed successfully!\\n"
                f"Found {len(sources)} stellar sources")
            
            logger.info("Photometry analysis completed successfully")
            return True
            
        except Exception as e:
            logger.error(f"Photometry error: {str(e)}")
            messagebox.showerror('Photometry Error', f'Analysis failed: {e}')
            return False
        finally:
            # Always clean up progress window
            if hasattr(G1, 'phot_barwindow') and G1.phot_barwindow:
                try:
                    G1.phot_barwindow.destroy()
                    G1.phot_barwindow = None
                except:
                    pass
            G1.update_status()

    def _display_photometry_results(self, image, apertures, num_sources):
        """Display photometry results visualization."""
        try:
            plt.figure(figsize=(12, 10))
            vmin, vmax = np.percentile(image, [1, 99])
            plt.imshow(image, cmap='gray_r', origin='lower', vmin=vmin, vmax=vmax)
            apertures.plot(color='red', lw=1.5, alpha=0.8)
            plt.title(f'Aperture Photometry Results ({num_sources} sources)')
            plt.xlabel('X (pixels)')
            plt.ylabel('Y (pixels)')
            plt.colorbar(label='Counts')
            plt.tight_layout()
            plt.show()
        except Exception as e:
            logger.error(f"Error displaying results: {str(e)}")

    def _save_photometry_results(self, sources, phot_table):
        """Save combined photometry results to a single CSV file."""
        try:
            # Ask user if they want to save results
            if not messagebox.askyesno("Save Results", "Would you like to save the photometry results to a CSV file?"):
                return
            
            # Get save location
            save_file = filedialog.asksaveasfilename(
                title="Save Photometry Results",
                defaultextension=".csv",
                filetypes=FILE_TYPES['csv']
            )
            
            if not save_file:
                return
            
            # Combine sources and photometry data by matching IDs
            import pandas as pd
            
            # Convert astropy tables to pandas DataFrames
            sources_df = sources.to_pandas()
            phot_df = phot_table.to_pandas()
            
            # Merge the dataframes on the id column (both should have the same number of rows)
            if len(sources_df) == len(phot_df):
                # Add index to ensure proper matching
                sources_df['source_id'] = range(len(sources_df))
                phot_df['source_id'] = range(len(phot_df))
                
                # Merge on source_id
                combined_df = pd.merge(sources_df, phot_df, on='source_id', how='inner')
                
                # Remove the temporary source_id column
                combined_df = combined_df.drop('source_id', axis=1)
                
                # Save to CSV
                combined_df.to_csv(save_file, index=False)
                
                logger.info(f"Combined photometry results saved: {save_file}")
                messagebox.showinfo("Results Saved", f"Photometry results saved successfully!\\n\\nFile: {os.path.basename(save_file)}\\nSources: {len(combined_df)}")
            else:
                logger.error("Sources and photometry tables have different lengths")
                messagebox.showerror('Data Error', 'Sources and photometry data have mismatched lengths.')
                
        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")
            messagebox.showerror('Save Error', f'Failed to save results: {e}')


class Conversion:
    """Enhanced conversion class with better error handling."""
    
    def raw_to_fits(self, raw_path, save_location):
        """Convert RAW to FITS with enhanced error handling and metadata preservation."""
        try:
            logger.info(f"Starting RAW conversion: {raw_path} -> {save_location}")
            
            # Read RAW file
            with rawpy.imread(raw_path) as raw:
                rgb = raw.postprocess(no_auto_bright=True, use_auto_wb=False, gamma=None)

            rgb_array = np.array(rgb)
            G1.progressbar['value'] = 25
            G1.progressbar.update()

            # Convert to grayscale for astronomical analysis
            gray = np.mean(rgb_array, axis=2)

            G1.progressbar['value'] = 50
            G1.progressbar.update()

            # Read EXIF metadata
            metadata = {}
            try:
                with open(raw_path, 'rb') as f:
                    tags = exifread.process_file(f)
                
                # Extract relevant metadata
                if 'EXIF ExposureTime' in tags:
                    metadata['EXPTIME'] = float(tags['EXIF ExposureTime'].values[0])
                if 'EXIF FNumber' in tags:
                    metadata['APERTURE'] = float(tags['EXIF FNumber'].values[0])
                if 'EXIF ISOSpeedRatings' in tags:
                    metadata['ISO'] = int(tags['EXIF ISOSpeedRatings'].values[0])
                if 'EXIF FocalLength' in tags:
                    metadata['FOC_LEN'] = float(tags['EXIF FocalLength'].values[0])
                if 'Image Make' in tags and 'Image Model' in tags:
                    metadata['CAMERA'] = f"{tags['Image Make'].values[0]} {tags['Image Model'].values[0]}"
                if 'EXIF DateTime' in tags:
                    metadata['DATE-OBS'] = str(tags['EXIF DateTime'].values[0])
                    
                logger.info(f"Extracted {len(metadata)} metadata fields")
            except Exception as e:
                logger.warning(f"Could not read EXIF data: {str(e)}")

            G1.progressbar['value'] = 75
            G1.progressbar.update()

            # Create FITS file with metadata
            hdu = fits.PrimaryHDU(data=gray)
            
            # Add standard headers
            hdu.header.set('ORIGIN', 'Digital Exoplanet Hunting Observatory')
            hdu.header.set('INSTRUME', 'Camera')
            hdu.header.set('OBJECT', 'Unknown')
            
            # Add extracted metadata
            for key, value in metadata.items():
                try:
                    hdu.header.set(key, value)
                except Exception as e:
                    logger.warning(f"Could not add header {key}: {str(e)}")

            # Save FITS file
            output_path = save_location if save_location.endswith('.fits') else f"{save_location}.fits"
            hdu.writeto(output_path, overwrite=True)

            G1.progressbar['value'] = 100
            G1.progressbar.update()
            
            # Update status
            G1.update_status()
            
            logger.info(f"RAW conversion completed: {output_path}")
            return True

        except Exception as e:
            logger.error(f"RAW conversion error: {str(e)}")
            messagebox.showerror('Conversion Error', f'Failed to convert RAW to FITS: {e}')
            
            # Update status on error
            G1.update_status()
            return False


class Fits:
    """Enhanced FITS handling class."""
    
    def load_fits(self, filepath):
        """Load and display FITS file with enhanced visualization."""
        try:
            logger.info(f"Loading FITS file: {filepath}")
            
            hdulist = fits.open(filepath)
            data = hdulist[0].data
            header = hdulist[0].header
            hdulist.close()
            
            if data is None:
                raise ValueError("No image data found in FITS file")
            
            print(f"FITS data shape: {data.shape}")
            print(f"Data type: {data.dtype}")
            print(f"Min value: {np.min(data):.2f}")
            print(f"Max value: {np.max(data):.2f}")
            print(f"Mean value: {np.mean(data):.2f}")
            
            max_point = np.unravel_index(np.argmax(data), data.shape)
            print(f"Maximum at position: {max_point}")
            
            # Enhanced visualization
            plt.figure(figsize=(12, 10))
            
            # Use better contrast scaling
            vmin, vmax = np.percentile(data, [1, 99])
            plt.imshow(data, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
            plt.xlabel('X (pixels)')
            plt.ylabel('Y (pixels)')
            plt.title(f'FITS Image: {os.path.basename(filepath)}')
            plt.colorbar(label='Counts')
            
            # Add metadata information if available
            if header:
                info_text = []
                if 'EXPTIME' in header:
                    info_text.append(f"Exposure: {header['EXPTIME']}s")
                if 'CAMERA' in header:
                    info_text.append(f"Camera: {header['CAMERA']}")
                if 'ISO' in header:
                    info_text.append(f"ISO: {header['ISO']}")
                
                if info_text:
                    plt.figtext(0.02, 0.02, '\\n'.join(info_text), fontsize=10,
                              bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
            
            plt.tight_layout()
            plt.show()
            
            return data
            
        except Exception as e:
            logger.error(f"Error loading FITS: {str(e)}")
            messagebox.showerror('FITS Error', f'Failed to load FITS file: {e}')
            return None

    def open_wsl_settings(self):
        """Open WSL Astrometry Settings dialog."""
        try:
            from photometry import Photometry
            
            # Create settings dialog
            settings_window = tk.Toplevel(self.root)
            settings_window.title("WSL Astrometry Settings")
            settings_window.geometry("600x500")
            settings_window.resizable(False, False)
            
            # Center the window
            settings_window.update_idletasks()
            x = (settings_window.winfo_screenwidth() // 2) - (settings_window.winfo_width() // 2)
            y = (settings_window.winfo_screenheight() // 2) - (settings_window.winfo_height() // 2)
            settings_window.geometry(f"+{x}+{y}")
            
            # Create a temporary photometry instance to access settings
            temp_photometry = Photometry()
            current_settings = temp_photometry.get_wsl_settings()
            
            # Main frame
            main_frame = tkk.Frame(settings_window, padding="20")
            main_frame.pack(fill=tk.BOTH, expand=True)
            
            # Title
            title_label = tk.Label(main_frame, text="WSL Astrometry Configuration", 
                                  font=('Arial', 14, 'bold'))
            title_label.pack(pady=(0, 20))
            
            # Settings variables
            workdir_var = tk.StringVar(value=current_settings.get('wsl_workdir', '$HOME/astrometry-work'))
            config_var = tk.StringVar(value=current_settings.get('wsl_config_file', '/home/burhan/astrometry-4211-4216.cfg'))
            scale_low_var = tk.DoubleVar(value=current_settings.get('scale_low', 4.5))
            scale_high_var = tk.DoubleVar(value=current_settings.get('scale_high', 5.5))
            
            # Work directory setting
            work_frame = tkk.LabelFrame(main_frame, text="WSL Work Directory", padding="10")
            work_frame.pack(fill=tk.X, pady=(0, 15))
            
            tk.Label(work_frame, text="Directory path (supports $HOME and ~/):").pack(anchor=tk.W)
            workdir_entry = tk.Entry(work_frame, textvariable=workdir_var, width=60)
            workdir_entry.pack(fill=tk.X, pady=(5, 0))
            
            # Config file setting
            config_frame = tkk.LabelFrame(main_frame, text="Astrometry Config File", padding="10")
            config_frame.pack(fill=tk.X, pady=(0, 15))
            
            tk.Label(config_frame, text="Path to astrometry.cfg file:").pack(anchor=tk.W)
            config_entry = tk.Entry(config_frame, textvariable=config_var, width=60)
            config_entry.pack(fill=tk.X, pady=(5, 0))
            
            # Scale settings
            scale_frame = tkk.LabelFrame(main_frame, text="Field of View Scale (degrees)", padding="10")
            scale_frame.pack(fill=tk.X, pady=(0, 20))
            
            scale_inner = tk.Frame(scale_frame)
            scale_inner.pack(fill=tk.X)
            
            tk.Label(scale_inner, text="Scale Low:").grid(row=0, column=0, sticky=tk.W, padx=(0, 10))
            scale_low_entry = tk.Entry(scale_inner, textvariable=scale_low_var, width=10)
            scale_low_entry.grid(row=0, column=1, padx=(0, 20))
            
            tk.Label(scale_inner, text="Scale High:").grid(row=0, column=2, sticky=tk.W, padx=(0, 10))
            scale_high_entry = tk.Entry(scale_inner, textvariable=scale_high_var, width=10)
            scale_high_entry.grid(row=0, column=3)
            
            # Status label
            status_label = tk.Label(main_frame, text="", fg="blue")
            status_label.pack(pady=(10, 0))
            
            # Button frame
            button_frame = tk.Frame(main_frame)
            button_frame.pack(fill=tk.X, pady=(20, 0))
            
            def save_settings():
                """Save the settings."""
                try:
                    new_settings = {
                        'wsl_workdir': workdir_var.get().strip(),
                        'wsl_config_file': config_var.get().strip(),
                        'scale_low': scale_low_var.get(),
                        'scale_high': scale_high_var.get(),
                        'scale_units': 'degwidth'
                    }
                    
                    success = temp_photometry.update_wsl_settings(new_settings)
                    if success:
                        status_label.config(text="✓ Settings saved successfully!", fg="green")
                        settings_window.after(2000, settings_window.destroy)
                    else:
                        status_label.config(text="✗ Failed to save settings.", fg="red")
                        
                except Exception as e:
                    logger.error(f"Error saving WSL settings: {e}")
                    status_label.config(text=f"✗ Error: {str(e)}", fg="red")
            
            def test_settings():
                """Test the current settings."""
                try:
                    status_label.config(text="Testing WSL configuration...", fg="blue")
                    settings_window.update()
                    
                    # Create temporary photometry with new settings
                    test_settings_dict = {
                        'wsl_workdir': workdir_var.get().strip(),
                        'wsl_config_file': config_var.get().strip(),
                        'scale_low': scale_low_var.get(),
                        'scale_high': scale_high_var.get(),
                        'scale_units': 'degwidth'
                    }
                    
                    test_photometry = Photometry(use_local_astrometry=True)
                    test_photometry.update_wsl_settings(test_settings_dict)
                    
                    # Run preflight check
                    if test_photometry._preflight_check_wsl():
                        status_label.config(text="✓ WSL configuration test passed!", fg="green")
                    else:
                        status_label.config(text="✗ WSL configuration test failed.", fg="red")
                        
                except Exception as e:
                    logger.error(f"Error testing WSL settings: {e}")
                    status_label.config(text=f"✗ Test error: {str(e)}", fg="red")
            
            # Buttons
            tk.Button(button_frame, text="Test Settings", command=test_settings).pack(side=tk.LEFT, padx=(0, 10))
            tk.Button(button_frame, text="Save", command=save_settings).pack(side=tk.LEFT, padx=(0, 10))
            tk.Button(button_frame, text="Cancel", command=settings_window.destroy).pack(side=tk.RIGHT)
            
        except Exception as e:
            logger.error(f"Error opening WSL settings dialog: {e}")
            messagebox.showerror("Settings Error", f"Failed to open settings dialog: {e}")

    def mean_fits(self, filepaths):
        """Stack multiple FITS files with enhanced error handling."""
        try:
            if not filepaths:
                raise ValueError("No files provided for stacking")
            
            logger.info(f"Stacking {len(filepaths)} FITS files")
            
            # Load first image
            hdulist = fits.open(filepaths[0])
            data = hdulist[0].data.astype(float)
            first_shape = data.shape
            hdulist.close()
            
            print(f"Stacking {len(filepaths)} images of shape {first_shape}")
            
            # Add remaining images
            for i, filepath in enumerate(filepaths[1:], 1):
                try:
                    hdulist = fits.open(filepath)
                    current_data = hdulist[0].data
                    hdulist.close()
                    
                    if current_data.shape != first_shape:
                        logger.error(f"Shape mismatch in file {i+1}: {current_data.shape} vs {first_shape}")
                        messagebox.showerror('Shape Mismatch', 
                            f'Image {i+1} has different dimensions.\\n'
                            f'Expected: {first_shape}\\nFound: {current_data.shape}')
                        return None
                    
                    data += current_data.astype(float)
                    
                except Exception as e:
                    logger.error(f"Error processing file {filepath}: {str(e)}")
                    continue

            # Calculate mean
            mean_data = data / len(filepaths)
            
            logger.info(f"Stacking completed. Mean data range: {np.min(mean_data):.2f} to {np.max(mean_data):.2f}")
            
            return mean_data
            
        except Exception as e:
            logger.error(f"Stacking error: {str(e)}")
            messagebox.showerror('Stacking Error', f'Failed to stack images: {e}')
            return None


def main():
    """Main application entry point with dependency checking."""
    try:
        # Log startup
        logger.info("="*60)
        logger.info(f"Starting {APP_CONFIG['name']} v{APP_CONFIG['version']}")
        logger.info("="*60)
        
        # Create instances
        global P1, C1, G1, F1
        P1 = Photometry()
        C1 = Conversion()
        G1 = GUI()
        F1 = Fits()
        
        # Start GUI
        G1.mainwindow()
        
    except Exception as e:
        logger.critical(f"Critical error in main: {str(e)}", exc_info=True)
        try:
            messagebox.showerror('Critical Error', f'Application failed to start: {e}')
        except:
            print(f"Critical error: {e}")


if __name__ == '__main__':
    main()