"""
Image to FITS conversion module for astronomical image processing.
Supports RAW, JPG, PNG, TIFF, and BMP formats.
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
    Class for converting various image formats to FITS format for astronomical processing.
    Supports RAW formats (CR2, CR3, NEF, ARW, DNG, RAF, ORF, RW2, PEF, SRW) 
    and standard formats (JPG, PNG, TIFF, BMP).
    """

    # Supported file extensions
    SUPPORTED_EXTENSIONS = {'.cr2', '.cr3', '.nef', '.arw', '.dng', '.raf', '.orf',
                          '.rw2', '.pef', '.srw', '.jpg', '.jpeg', '.png', '.tiff', '.tif', '.bmp'}

    def __init__(self, gui_reference=None):
        """
        Initialize the Conversion class.

        Args:
            gui_reference: Reference to the GUI instance for progress updates
        """
        self.gui_ref = gui_reference

    def convert_to_fits(self, image_path, save_location):
        """
        Convert an image file to FITS format with metadata preservation.
        Supports RAW formats and standard image formats (JPG, PNG, TIFF, BMP).

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
            if file_extension in {'.jpg', '.jpeg', '.png', '.tiff', '.tif', '.bmp'}:
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
        Convert standard image formats (JPG, PNG, TIFF, BMP) to FITS format with metadata preservation.

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

        # Read metadata (EXIF for JPG/TIFF, limited for PNG/BMP)
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
            # Ensure image data is in the correct format for astrometry.net
            # Convert to float32 and ensure proper orientation
            if len(image_data.shape) == 2:
                # Ensure data is properly scaled (0-65535 for 16-bit equivalent)
                if image_data.max() <= 1.0:  # If normalized to 0-1
                    processed_data = (image_data * 65535).astype(np.uint16)
                elif image_data.max() <= 255:  # If 8-bit
                    processed_data = (image_data * 257).astype(np.uint16)  # Scale to 16-bit
                else:
                    processed_data = image_data.astype(np.uint16)
            else:
                raise ValueError("Image data must be 2D for FITS conversion")

            # Create FITS file with proper headers
            hdu = fits.PrimaryHDU(data=processed_data)

            # Set essential FITS headers first
            header = hdu.header
            header['SIMPLE'] = True
            header['BITPIX'] = 16  # 16-bit unsigned integer
            header['NAXIS'] = 2
            header['NAXIS1'] = processed_data.shape[1]  # Width
            header['NAXIS2'] = processed_data.shape[0]  # Height
            header['BZERO'] = 32768  # For unsigned 16-bit
            header['BSCALE'] = 1

            # Add custom metadata
            self._add_metadata_to_header(header, metadata)

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
        Read metadata from standard image files (JPG, PNG, TIFF, BMP) using PIL.

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

                # Try to get EXIF data (mainly for JPG and TIFF)
                exif_data = None
                try:
                    exif_data = img._getexif()
                except AttributeError:
                    # Some formats like BMP don't support _getexif()
                    logger.debug(f"No EXIF data available for {os.path.splitext(image_path)[1]} format")

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

            logger.info(f"Extracted {len(metadata)} metadata fields from {os.path.splitext(image_path)[1].upper()} image")

        except Exception as e:
            logger.warning(f"Failed to read metadata from standard image: {str(e)}")

        return metadata

    def _add_metadata_to_header(self, header, metadata):
        """
        Add metadata to FITS header in astrometry.net compatible format.

        Args:
            header: FITS header object
            metadata (dict): Dictionary containing metadata
        """
        # Remove any existing WCS keywords that might conflict
        wcs_keywords = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2',
                       'CROTA1', 'CROTA2', 'CTYPE1', 'CTYPE2', 'CD1_1', 'CD1_2',
                       'CD2_1', 'CD2_2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2']

        for keyword in wcs_keywords:
            if keyword in header:
                del header[keyword]

        # Add standard FITS keywords
        header.set('ORIGIN', 'Digital Exoplanet Hunting Observatory', 'File origin')
        header.set('INSTRUME', 'Camera', 'Instrument used')
        header.set('OBJECT', 'Unknown', 'Target object')

        # Set image scaling information for proper display
        header.set('DATAMIN', 0, 'Minimum data value')
        header.set('DATAMAX', 65535, 'Maximum data value')

        # Add extracted metadata
        for key, value in metadata.items():
            try:
                if key == 'EXPTIME':
                    header.set(key, float(value), 'Exposure time (seconds)')
                elif key == 'APERTURE':
                    header.set('APERTURE', float(value), 'F-number of lens')
                elif key == 'ISO':
                    header.set(key, int(value), 'ISO sensitivity')
                elif key == 'FOC_LEN':
                    header.set('FOCALLEN', float(value), 'Focal length (mm)')
                elif key == 'CAMERA':
                    # Limit string length for FITS compatibility
                    camera_str = str(value)[:67]  # FITS strings max 68 chars
                    header.set(key, camera_str, 'Camera make and model')
                elif key == 'DATE-OBS':
                    # Ensure proper date format
                    date_str = str(value).replace(':', '-', 2)  # Replace first two colons
                    header.set(key, date_str, 'Observation date')
                elif key == 'WIDTH':
                    # Don't add WIDTH/HEIGHT as separate keywords since NAXIS1/2 handle this
                    continue
                elif key == 'HEIGHT':
                    continue
                elif key == 'FORMAT':
                    header.set('SRCFMT', str(value)[:67], 'Source image format')
                else:
                    # Add other metadata with length limits
                    if isinstance(value, str):
                        value = value[:67]
                    header.set(key[:8], value, 'Image metadata')  # FITS keyword max 8 chars

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

    def validate_fits_file(self, fits_path):
        """
        Validate FITS file for astrometry.net compatibility.

        Args:
            fits_path (str): Path to the FITS file

        Returns:
            tuple: (is_valid, error_messages)
        """
        errors = []

        try:
            with fits.open(fits_path) as hdul:
                header = hdul[0].header
                data = hdul[0].data

                # Check essential headers
                if 'NAXIS' not in header:
                    errors.append("Missing NAXIS keyword")
                elif header['NAXIS'] != 2:
                    errors.append(f"NAXIS should be 2, found {header['NAXIS']}")

                if 'NAXIS1' not in header or 'NAXIS2' not in header:
                    errors.append("Missing NAXIS1 or NAXIS2")

                if 'BITPIX' not in header:
                    errors.append("Missing BITPIX keyword")

                # Check data shape matches header
                if data is not None:
                    expected_shape = (header.get('NAXIS2', 0), header.get('NAXIS1', 0))
                    if data.shape != expected_shape:
                        errors.append(f"Data shape {data.shape} doesn't match header {expected_shape}")

                # Check for conflicting WCS keywords
                wcs_keywords = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2']
                partial_wcs = any(kw in header for kw in wcs_keywords)
                if partial_wcs and not all(kw in header for kw in wcs_keywords[:4]):
                    errors.append("Incomplete WCS keywords found - may cause issues")

                logger.info(f"FITS validation complete. Errors found: {len(errors)}")

        except Exception as e:
            errors.append(f"Failed to read FITS file: {str(e)}")

        return len(errors) == 0, errors

    def create_astrometry_compatible_fits(self, image_path, save_location, force_16bit=True):
        """
        Create a FITS file specifically optimized for astrometry.net compatibility.

        Args:
            image_path (str): Path to the source image
            save_location (str): Path where the FITS file should be saved
            force_16bit (bool): Force 16-bit output for maximum compatibility

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            logger.info(f"Creating astrometry.net compatible FITS: {image_path}")

            # Determine file type and read image
            file_extension = os.path.splitext(image_path)[1].lower()

            if file_extension in {'.jpg', '.jpeg', '.png', '.tiff', '.tif', '.bmp'}:
                # Read with PIL
                with Image.open(image_path) as img:
                    if img.mode in ('RGBA', 'LA'):
                        # Handle transparency
                        background = Image.new('RGB', img.size, (0, 0, 0))  # Black background
                        if img.mode == 'RGBA':
                            background.paste(img, mask=img.split()[-1])
                        else:
                            background.paste(img, mask=img.split()[-1])
                        img = background
                    elif img.mode != 'RGB':
                        img = img.convert('RGB')

                    img_array = np.array(img)

                # Convert to grayscale using proper luminance weights
                if len(img_array.shape) == 3:
                    gray = np.dot(img_array[...,:3], [0.299, 0.587, 0.114])
                else:
                    gray = img_array

            else:
                # Handle RAW files
                with rawpy.imread(image_path) as raw:
                    rgb = raw.postprocess(
                        no_auto_bright=True,
                        use_auto_wb=False,
                        gamma=None,
                        output_bps=16  # Force 16-bit output
                    )
                rgb_array = np.array(rgb)
                gray = np.mean(rgb_array, axis=2)

            self._update_progress(30)

            # Ensure proper data scaling for astrometry.net
            if gray.dtype == np.uint8 or gray.max() <= 255:
                # Scale 8-bit to 16-bit
                processed_data = (gray.astype(np.float32) * 257).astype(np.uint16)
            elif gray.max() <= 1.0:
                # Scale normalized data to 16-bit
                processed_data = (gray * 65535).astype(np.uint16)
            else:
                # Assume already in proper range
                processed_data = gray.astype(np.uint16)

            self._update_progress(60)

            # Create minimal FITS header for astrometry.net
            hdu = fits.PrimaryHDU(data=processed_data)
            header = hdu.header

            # Essential headers only - minimal approach
            header['SIMPLE'] = True
            header['BITPIX'] = 16
            header['NAXIS'] = 2
            header['NAXIS1'] = processed_data.shape[1]
            header['NAXIS2'] = processed_data.shape[0]
            header['BZERO'] = 32768
            header['BSCALE'] = 1.0

            # Minimal metadata to avoid conflicts
            header['OBJECT'] = 'UNKNOWN'
            header['ORIGIN'] = 'Custom Converter'

            self._update_progress(80)

            # Save with minimal header
            output_path = save_location if save_location.endswith('.fits') else f"{save_location}.fits"
            hdu.writeto(output_path, overwrite=True)

            # Validate the created file
            is_valid, validation_errors = self.validate_fits_file(output_path)
            if not is_valid:
                logger.warning(f"FITS validation warnings: {validation_errors}")

            self._update_progress(100)
            self._cleanup_progress()

            logger.info(f"Astrometry-compatible FITS created: {output_path}")
            messagebox.showinfo('Conversion Complete',
                              f'Successfully converted to astrometry.net compatible FITS!\n{output_path}')
            return True

        except Exception as e:
            logger.error(f"Error creating astrometry-compatible FITS: {str(e)}")
            messagebox.showerror('Conversion Error',
                               f'Failed to create astrometry-compatible FITS: {e}')
            self._cleanup_progress()
            return False

    def batch_convert_to_fits(self, image_paths, output_directory):
        """
        Convert multiple images to FITS format in batch.

        Args:
            image_paths (list): List of image file paths
            output_directory (str): Directory to save converted FITS files

        Returns:
            tuple: (successful_count, failed_count, failed_files)
        """
        successful_count = 0
        failed_count = 0
        failed_files = []

        total_files = len(image_paths)
        logger.info(f"Starting batch conversion of {total_files} images")

        for i, image_path in enumerate(image_paths):
            try:
                # Create output filename
                base_name = os.path.splitext(os.path.basename(image_path))[0]
                output_path = os.path.join(output_directory, f"{base_name}.fits")

                # Check if format is supported
                if not self.is_supported_format(image_path):
                    failed_count += 1
                    failed_files.append(f"{os.path.basename(image_path)} (unsupported format)")
                    logger.warning(f"Unsupported format: {image_path}")
                    continue

                # Update progress for batch operation
                progress = int((i / total_files) * 100)
                self._update_progress(progress)

                # Convert the image
                if self.convert_to_fits(image_path, output_path):
                    successful_count += 1
                    logger.info(f"Batch conversion successful: {os.path.basename(image_path)}")
                else:
                    failed_count += 1
                    failed_files.append(os.path.basename(image_path))
                    logger.warning(f"Batch conversion failed: {os.path.basename(image_path)}")

            except Exception as e:
                failed_count += 1
                failed_files.append(os.path.basename(image_path))
                logger.error(f"Error in batch conversion of {os.path.basename(image_path)}: {str(e)}")

        self._update_progress(100)
        self._cleanup_progress()

        logger.info(f"Batch conversion complete: {successful_count} successful, {failed_count} failed")
        
        # Show completion message
        if failed_count == 0:
            messagebox.showinfo('Batch Conversion Complete', 
                              f'Successfully converted all {successful_count} images to FITS format!')
        else:
            messagebox.showwarning('Batch Conversion Complete', 
                                 f'Converted {successful_count} images successfully.\n'
                                 f'{failed_count} images failed: {", ".join(failed_files[:5])}'
                                 f'{"..." if len(failed_files) > 5 else ""}')

        return successful_count, failed_count, failed_files