"""
FITS file handling module for loading and processing astronomical images.
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from tkinter import messagebox
import logging

logger = logging.getLogger(__name__)


class FitsHandler:
    """
    Class for handling FITS file operations including loading and stacking.
    """
    
    def load_fits(self, filepath):
        """
        Load a FITS file and display it.
        
        Args:
            filepath (str): Path to the FITS file
            
        Returns:
            numpy.ndarray: Image data array or None if failed
        """
        try:
            logger.info(f"Loading FITS file: {filepath}")
            
            with fits.open(filepath) as hdulist:
                data = hdulist[0].data
                header = hdulist[0].header
            
            if data is None:
                messagebox.showerror('Load Error', 'No image data found in FITS file')
                return None
            
            # Display basic information
            logger.info(f"Image dimensions: {data.shape}")
            logger.info(f"Data type: {data.dtype}")
            
            # Find maximum brightness point
            max_point = np.unravel_index(np.argmax(data), data.shape)
            max_value = data[max_point]
            
            logger.info(f"Maximum value: {max_value} at position {max_point}")
            
            # Display the image
            self._display_fits_image(data, filepath, header)
            
            return data
            
        except Exception as e:
            logger.error(f"Error loading FITS file: {str(e)}")
            messagebox.showerror('Load Error', f'Failed to load FITS file: {e}')
            return None
    
    def mean_fits(self, filepaths):
        """
        Create a mean-combined image from multiple FITS files (stacking).
        
        Args:
            filepaths (list): List of FITS file paths
            
        Returns:
            numpy.ndarray: Mean-combined image data or None if failed
        """
        try:
            if not filepaths:
                messagebox.showwarning('No Files', 'No FITS files selected for stacking')
                return None
            
            logger.info(f"Stacking {len(filepaths)} FITS files")
            
            # Load first image to get dimensions
            with fits.open(filepaths[0]) as hdulist:
                data = hdulist[0].data.astype(float)
                first_shape = data.shape
            
            # Verify all images have the same dimensions
            for i, filepath in enumerate(filepaths[1:], 1):
                with fits.open(filepath) as hdulist:
                    current_data = hdulist[0].data
                    if current_data.shape != first_shape:
                        messagebox.showerror(
                            'Dimension Mismatch', 
                            f'Image {i+1} has different dimensions ({current_data.shape}) '
                            f'than the first image ({first_shape})'
                        )
                        return None
                    data += current_data.astype(float)
            
            # Calculate mean
            mean_data = data / len(filepaths)
            
            logger.info(f"Stacking completed. Final image shape: {mean_data.shape}")
            
            # Display the stacked image
            self._display_stacked_image(mean_data, len(filepaths))
            
            return mean_data
            
        except Exception as e:
            logger.error(f"Error in FITS stacking: {str(e)}")
            messagebox.showerror('Stacking Error', f'Failed to stack FITS files: {e}')
            return None
    
    def _display_fits_image(self, data, filepath, header=None):
        """
        Display a FITS image with proper scaling and information.
        
        Args:
            data (numpy.ndarray): Image data
            filepath (str): Original file path
            header: FITS header (optional)
        """
        try:
            plt.figure(figsize=(10, 8))
            
            # Use appropriate scaling for astronomical images
            vmin, vmax = np.percentile(data, [1, 99])
            
            plt.imshow(data, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
            plt.xlabel('X (pixels)')
            plt.ylabel('Y (pixels)')
            plt.title(f'FITS Image: {filepath.split("/")[-1]}')
            plt.colorbar(label='Counts')
            
            # Add header information if available
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
                              bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
            
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            logger.error(f"Error displaying FITS image: {str(e)}")
    
    def _display_stacked_image(self, data, num_images):
        """
        Display a stacked FITS image.
        
        Args:
            data (numpy.ndarray): Stacked image data
            num_images (int): Number of images that were stacked
        """
        try:
            plt.figure(figsize=(12, 10))
            
            # Use appropriate scaling
            vmin, vmax = np.percentile(data, [1, 99])
            
            plt.imshow(data, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
            plt.xlabel('X (pixels)')
            plt.ylabel('Y (pixels)')
            plt.title(f'Stacked Image ({num_images} frames)')
            plt.colorbar(label='Mean Counts')
            
            # Add statistics
            mean_val = np.mean(data)
            std_val = np.std(data)
            info_text = f"Mean: {mean_val:.2f}\\nStd: {std_val:.2f}\\nFrames: {num_images}"
            
            plt.figtext(0.02, 0.02, info_text, fontsize=10,
                      bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
            
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            logger.error(f"Error displaying stacked image: {str(e)}")
    
    def save_fits(self, data, filepath, header=None):
        """
        Save image data as a FITS file.
        
        Args:
            data (numpy.ndarray): Image data to save
            filepath (str): Output file path
            header: Optional FITS header
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            hdu = fits.PrimaryHDU(data=data, header=header)
            output_path = filepath if filepath.endswith('.fits') else f"{filepath}.fits"
            hdu.writeto(output_path, overwrite=True)
            
            logger.info(f"FITS file saved: {output_path}")
            messagebox.showinfo('Save Complete', f'FITS file saved successfully!\\n{output_path}')
            return True
            
        except Exception as e:
            logger.error(f"Error saving FITS file: {str(e)}")
            messagebox.showerror('Save Error', f'Failed to save FITS file: {e}')
            return False