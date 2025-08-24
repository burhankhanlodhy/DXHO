"""
RAW/JPG/PNG to FITS conversion module for astronomical image processing.
"""

import numpy as np
from astropy.io import fits
import rawpy
import exifread
from PIL import Image
from PIL.ExifTags import TAGS
from tkinter import messagebox
import logging
import os

logger = logging.getLogger(__name__)


class Conversion:
    """
    Class for converting RAW, JPG, and PNG astronomical images to FITS format.
    """

    # Supported file extensions
    SUPPORTED_EXTENSIONS = {'.cr2', '.cr3', '.nef', '.arw', '.dng', '.raf', '.orf',
                          '.rw2', '.pef', '.srw', '.jpg', '.jpeg', '.png'}

    def __init__(self, gui_reference=None):
        """
        Initialize the Conversion class.

        Args:
            gui_reference: Reference to the GUI instance for progress updates
        """
        self.gui_ref = gui_reference

    def convert_to_fits(self, image_path, save_location):
        """
        Convert an image file (RAW, JPG, or PNG) to FITS format with metadata preservation.

        Args:
            image_path (str): Path to the image file
            save_location (str): Path where the FITS file should be saved

        Returns:
            bool: True if conversion successful, False otherwise
        """
        try:
            # Determine file type
            file_extension = os.path.splitext(image_path)[1].lower()

            if file_extension not in self.SUPPORTED_EXTENSIONS:
                raise ValueError(f"Unsupported file format: {file_extension}")

            logger.info(f"Starting image to FITS conversion: {image_path}")

            # Process based on file type
            if file_extension in {'.jpg', '.jpeg', '.png'}:
                return self._convert_standard_image_to_fits(image_path, save_location)
            else:
                return self._convert_raw_to_fits(image_path, save_location)

        except Exception as e:
            logger.error(f"Error in image to FITS conversion: {str(e)}")
            messagebox.showerror('Conversion Error', f'Failed to convert image to FITS: {e}')
            self._cleanup_progress()
            return False

    def _convert_raw_to_fits(self, raw_path, save_location):
        """
        Convert a RAW image file to FITS format with metadata preservation.

        Args:
            raw_path (str): Path to the RAW image file
            save_location (str): Path where the FITS file should be saved

        Returns:
            bool: True if conversion successful, False otherwise
        """
        # Read the RAW file
        with rawpy.imread(raw_path) as raw:
            rgb = raw.postprocess(
                no_auto_bright=True,
                use_auto_wb=False,
                gamma=None
            )

        rgb_array = np.array(rgb)
        self._update_progress(25)

        # Convert to grayscale (astronomical images are typically monochrome)
        gray = np.mean(rgb_array, axis=2)
        self._update_progress(50)

        # Read EXIF metadata
        metadata = self._read_exif_metadata(raw_path)
        self._update_progress(75)

        # Create and save FITS file
        return self._create_fits_file(gray, metadata, save_location)

    def _convert_standard_image_to_fits(self, image_path, save_location):
        """
        Convert a JPG or PNG image file to FITS format with metadata preservation.

        Args:
            image_path (str): Path to the image file
            save_location (str): Path where the FITS file should be saved

        Returns:
            bool: True if conversion successful, False otherwise
        """
        # Read the image file
        with Image.open(image_path) as img:
            # Convert to RGB if necessary (handles PNG with transparency)
            if img.mode in ('RGBA', 'LA'):
                # Create white background for transparent images
                background = Image.new('RGB', img.size, (255, 255, 255))
                if img.mode == 'RGBA':
                    background.paste(img, mask=img.split()[-1])
                else:  # LA mode
                    background.paste(img, mask=img.split()[-1])
                img = background
            elif img.mode != 'RGB':
                img = img.convert('RGB')

            img_array = np.array(img)

        self._update_progress(25)

        # Convert to grayscale for astronomical processing
        if len(img_array.shape) == 3:
            gray = np.mean(img_array, axis=2)
        else:
            gray = img_array  # Already grayscale

        self._update_progress(50)

        # Read metadata (EXIF for JPG, limited for PNG)
        metadata = self._read_standard_image_metadata(image_path)
        self._update_progress(75)

        # Create and save FITS file
        return self._create_fits_file(gray, metadata, save_location)

    def _create_fits_file(self, image_data, metadata, save_location):
        """
        Create FITS file from image data and metadata.

        Args:
            image_data (np.array): Image data array
            metadata (dict): Metadata dictionary
            save_location (str): Path where the FITS file should be saved

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Create FITS file with metadata
            hdu = fits.PrimaryHDU(data=image_data.astype(np.float32))
            self._add_metadata_to_header(hdu.header, metadata)

            # Save FITS file
            output_path = save_location if save_location.endswith('.fits') else f"{save_location}.fits"
            hdu.writeto(output_path, overwrite=True)

            self._update_progress(100)
            self._cleanup_progress()

            logger.info(f"Conversion completed successfully: {output_path}")
            messagebox.showinfo('Conversion Complete', f'Successfully converted to FITS format!\n{output_path}')
            return True

        except Exception as e:
            logger.error(f"Error creating FITS file: {str(e)}")
            raise

    def _read_exif_metadata(self, raw_path):
        """
        Read EXIF metadata from RAW file using exifread.

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
                metadata['CAMERA'] = f"{tags['Image Make']} {tags['Image Model']}"

            if 'EXIF DateTime' in tags:
                metadata['DATE-OBS'] = str(tags['EXIF DateTime'])

            logger.info(f"Extracted {len(metadata)} metadata fields from RAW")

        except Exception as e:
            logger.warning(f"Failed to read EXIF metadata from RAW: {str(e)}")

        return metadata

    def _read_standard_image_metadata(self, image_path):
        """
        Read metadata from JPG or PNG files using PIL.

        Args:
            image_path (str): Path to the image file

        Returns:
            dict: Dictionary containing extracted metadata
        """
        metadata = {}

        try:
            with Image.open(image_path) as img:
                # Get basic image info
                metadata['WIDTH'] = img.size[0]
                metadata['HEIGHT'] = img.size[1]
                metadata['FORMAT'] = img.format

                # Try to get EXIF data (mainly for JPG)
                exif_data = img._getexif()

                if exif_data:
                    for tag_id, value in exif_data.items():
                        tag = TAGS.get(tag_id, tag_id)

                        try:
                            if tag == 'ExposureTime':
                                if hasattr(value, 'numerator') and hasattr(value, 'denominator'):
                                    metadata['EXPTIME'] = float(value.numerator) / float(value.denominator)
                                else:
                                    metadata['EXPTIME'] = float(value)

                            elif tag == 'FNumber':
                                if hasattr(value, 'numerator') and hasattr(value, 'denominator'):
                                    metadata['APERTURE'] = float(value.numerator) / float(value.denominator)
                                else:
                                    metadata['APERTURE'] = float(value)

                            elif tag == 'ISOSpeedRatings':
                                metadata['ISO'] = int(value)

                            elif tag == 'FocalLength':
                                if hasattr(value, 'numerator') and hasattr(value, 'denominator'):
                                    metadata['FOC_LEN'] = float(value.numerator) / float(value.denominator)
                                else:
                                    metadata['FOC_LEN'] = float(value)

                            elif tag == 'Make':
                                metadata['MAKE'] = str(value)

                            elif tag == 'Model':
                                metadata['MODEL'] = str(value)

                            elif tag == 'DateTime':
                                metadata['DATE-OBS'] = str(value)

                        except (ValueError, AttributeError, TypeError) as e:
                            logger.debug(f"Could not process EXIF tag {tag}: {e}")
                            continue

                    # Combine make and model if both exist
                    if 'MAKE' in metadata and 'MODEL' in metadata:
                        metadata['CAMERA'] = f"{metadata['MAKE']} {metadata['MODEL']}"
                        del metadata['MAKE'], metadata['MODEL']

            logger.info(f"Extracted {len(metadata)} metadata fields from standard image")

        except Exception as e:
            logger.warning(f"Failed to read metadata from standard image: {str(e)}")

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
                elif key == 'WIDTH':
                    header.set('NAXIS1', value, 'Image width')
                elif key == 'HEIGHT':
                    header.set('NAXIS2', value, 'Image height')
                elif key == 'FORMAT':
                    header.set('SRCFMT', str(value), 'Source image format')
                else:
                    # Add other metadata with generic comment
                    header.set(key, value, 'Image metadata')

            except Exception as e:
                logger.warning(f"Failed to add metadata {key}: {str(e)}")

        logger.info("Metadata added to FITS header")

    def _update_progress(self, value):
        """Update progress bar if GUI reference exists."""
        if self.gui_ref and hasattr(self.gui_ref, 'progressbar'):
            self.gui_ref.progressbar['value'] = value
            self.gui_ref.progressbar.update()

    def _cleanup_progress(self):
        """Clean up progress bar and status."""
        if self.gui_ref:
            if hasattr(self.gui_ref, 'progress_window') and self.gui_ref.progress_window:
                self.gui_ref.progress_window.destroy()
                self.gui_ref.progress_window = None
            if hasattr(self.gui_ref, 'statusbar'):
                self.gui_ref.statusbar['text'] = "Ready"
                self.gui_ref.statusbar.update()

    # Backward compatibility - keep the original method name
    def raw_to_fits(self, raw_path, save_location):
        """
        Legacy method for RAW to FITS conversion.
        Redirects to the new convert_to_fits method.
        """
        return self.convert_to_fits(raw_path, save_location)

    @classmethod
    def is_supported_format(cls, file_path):
        """
        Check if the file format is supported for conversion.

        Args:
            file_path (str): Path to the file

        Returns:
            bool: True if format is supported, False otherwise
        """
        file_extension = os.path.splitext(file_path)[1].lower()
        return file_extension in cls.SUPPORTED_EXTENSIONS

    @classmethod
    def get_supported_extensions(cls):
        """
        Get list of supported file extensions.

        Returns:
            set: Set of supported file extensions
        """
        return cls.SUPPORTED_EXTENSIONS.copy()