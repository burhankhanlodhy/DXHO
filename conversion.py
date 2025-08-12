"""
RAW to FITS conversion module for astronomical image processing.
"""

import numpy as np
from astropy.io import fits
import rawpy
import exifread
from tkinter import messagebox
import logging

logger = logging.getLogger(__name__)


class Conversion:
    """
    Class for converting RAW astronomical images to FITS format.
    """
    
    def __init__(self, gui_reference=None):
        """
        Initialize the Conversion class.
        
        Args:
            gui_reference: Reference to the GUI instance for progress updates
        """
        self.gui_ref = gui_reference
    
    def raw_to_fits(self, raw_path, save_location):
        """
        Convert a RAW image file to FITS format with metadata preservation.
        
        Args:
            raw_path (str): Path to the RAW image file
            save_location (str): Path where the FITS file should be saved
            
        Returns:
            bool: True if conversion successful, False otherwise
        """
        try:
            logger.info(f"Starting RAW to FITS conversion: {raw_path}")
            
            # Read the RAW file
            with rawpy.imread(raw_path) as raw:
                rgb = raw.postprocess(
                    no_auto_bright=True, 
                    use_auto_wb=False, 
                    gamma=None
                )
            
            rgb_array = np.array(rgb)
            
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar'):
                self.gui_ref.progressbar['value'] = 25
                self.gui_ref.progressbar.update()
            
            # Convert to grayscale (astronomical images are typically monochrome)
            gray = np.mean(rgb_array, axis=2)
            
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar'):
                self.gui_ref.progressbar['value'] = 50
                self.gui_ref.progressbar.update()
            
            # Read EXIF metadata
            metadata = self._read_exif_metadata(raw_path)
            
            if self.gui_ref and hasattr(self.gui_ref, 'progressbar'):
                self.gui_ref.progressbar['value'] = 75
                self.gui_ref.progressbar.update()
            
            # Create FITS file with metadata
            hdu = fits.PrimaryHDU(data=gray)
            self._add_metadata_to_header(hdu.header, metadata)
            
            # Save FITS file
            output_path = save_location if save_location.endswith('.fits') else f"{save_location}.fits"
            hdu.writeto(output_path, overwrite=True)
            
            if self.gui_ref:
                if hasattr(self.gui_ref, 'progressbar'):
                    self.gui_ref.progressbar['value'] = 100
                    self.gui_ref.progressbar.update()
                if hasattr(self.gui_ref, 'progress_window') and self.gui_ref.progress_window:
                    self.gui_ref.progress_window.destroy()
                    self.gui_ref.progress_window = None
                if hasattr(self.gui_ref, 'statusbar'):
                    self.gui_ref.statusbar['text'] = "Ready"
                    self.gui_ref.statusbar.update()
            
            logger.info(f"Conversion completed successfully: {output_path}")
            messagebox.showinfo('Conversion Complete', f'Successfully converted to FITS format!\\n{output_path}')
            return True
            
        except Exception as e:
            logger.error(f"Error in RAW to FITS conversion: {str(e)}")
            messagebox.showerror('Conversion Error', f'Failed to convert RAW to FITS: {e}')
            
            # Clean up progress bar if conversion fails
            if self.gui_ref:
                if hasattr(self.gui_ref, 'progress_window') and self.gui_ref.progress_window:
                    self.gui_ref.progress_window.destroy()
                    self.gui_ref.progress_window = None
                if hasattr(self.gui_ref, 'statusbar'):
                    self.gui_ref.statusbar['text'] = "Ready"
                    self.gui_ref.statusbar.update()
            
            return False
    
    def _read_exif_metadata(self, raw_path):
        """
        Read EXIF metadata from RAW file.
        
        Args:
            raw_path (str): Path to the RAW file
            
        Returns:
            dict: Dictionary containing extracted metadata
        """
        metadata = {}
        
        try:
            with open(raw_path, 'rb') as f:
                tags = exifread.process_file(f)
            
            # Extract relevant EXIF tags
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
            logger.warning(f"Failed to read EXIF metadata: {str(e)}")
        
        return metadata
    
    def _add_metadata_to_header(self, header, metadata):
        """
        Add metadata to FITS header.
        
        Args:
            header: FITS header object
            metadata (dict): Dictionary containing metadata
        """
        # Add standard FITS keywords
        header.set('ORIGIN', 'Digital Exoplanet Hunting Observatory', 'File origin')
        header.set('INSTRUME', 'Camera', 'Instrument used')
        header.set('OBJECT', 'Unknown', 'Target object')
        
        # Add extracted metadata
        for key, value in metadata.items():
            try:
                if key == 'EXPTIME':
                    header.set(key, value, 'Exposure time (seconds)')
                elif key == 'APERTURE':
                    header.set(key, value, 'F-number of lens')
                elif key == 'ISO':
                    header.set(key, value, 'ISO sensitivity')
                elif key == 'FOC_LEN':
                    header.set(key, value, 'Focal length (mm)')
                elif key == 'CAMERA':
                    header.set(key, str(value), 'Camera make and model')
                elif key == 'DATE-OBS':
                    header.set(key, str(value), 'Observation date')
                else:
                    header.set(key, value)
                    
            except Exception as e:
                logger.warning(f"Failed to add metadata {key}: {str(e)}")
        
        logger.info("Metadata added to FITS header")