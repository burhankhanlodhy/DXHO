# Digital Exoplanet Hunting Observatory

A comprehensive astronomical image processing application for stellar photometry and exoplanet detection.

## Features

- **RAW Image Processing**: Convert astronomical RAW images (CR2, NEF, ARW, DNG) to FITS format
- **Image Stacking**: Combine multiple FITS images to reduce noise and improve signal quality
- **Aperture Photometry**: Perform precise stellar photometry analysis using DAOStarFinder
- **Metadata Preservation**: Maintain EXIF data during RAW to FITS conversion
- **User-Friendly GUI**: Intuitive interface with progress tracking and help system

## Installation

1. **Install Python 3.8 or higher**

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the application**:
   ```bash
   python main.py
   ```

## Usage

### 1. Opening Files
- Use **File → Open Image** for standard image formats (PNG, TIFF, JPEG)
- Use **File → Open RAW File** for astronomical RAW images
- Use **File → Open FITS File** for FITS astronomical images

### 2. RAW to FITS Conversion
1. Load a RAW file (File → Open RAW File)
2. Select **Process → Convert RAW to FITS**
3. Choose save location
4. EXIF metadata will be preserved in the FITS header

### 3. Image Stacking
1. Select **Process → Stack FITS Images**
2. Choose multiple FITS files
3. The application creates a mean-combined image to reduce noise

### 4. Aperture Photometry
1. Load a FITS file for analysis
2. Select **Analysis → Aperture Photometry**
3. Set parameters:
   - **FWHM**: Full Width Half Maximum (typical: 3-8 pixels)
   - **Threshold**: Detection threshold in sigma units (typical: 2-5)
4. Results are saved as CSV files containing source positions and photometry

## File Structure

```
├── main.py              # Main application entry point
├── gui.py               # Main GUI interface
├── photometry.py        # Aperture photometry module
├── conversion.py        # RAW to FITS conversion
├── fits_handler.py      # FITS file operations
├── config.py            # Configuration and logging setup
├── requirements.txt     # Python dependencies
├── Class.spec           # PyInstaller build configuration
└── README.md           # This file
```

## Dependencies

- **NumPy**: Numerical computations
- **Pandas**: Data manipulation
- **Astropy**: Astronomical data handling
- **Photutils**: Photometry utilities
- **Matplotlib**: Plotting and visualization
- **rawpy**: RAW image processing
- **exifread**: EXIF metadata extraction
- **Pillow**: Image processing

## Building Executable

To create a standalone executable:

```bash
pyinstaller Class.spec
```

## Tips

1. **Image Quality**: Use FITS format for scientific analysis
2. **Noise Reduction**: Stack multiple images of the same target
3. **Photometry Parameters**: Adjust FWHM and threshold based on your image characteristics
4. **File Organization**: Keep RAW files and processed FITS files organized in separate folders

## Troubleshooting

### Common Issues

1. **Missing Dependencies**: Run `pip install -r requirements.txt`
2. **RAW File Support**: Ensure LibRaw is properly installed (included with rawpy)
3. **Memory Issues**: For large images, close other applications to free RAM
4. **File Permissions**: Ensure write permissions for output directories

### Log Files

Application logs are saved in the `logs/` directory with timestamps for debugging.

## Version History

### v2.0
- Removed authentication system
- Modular architecture with separate files
- Improved error handling and logging
- Enhanced progress tracking
- Better user interface
- Comprehensive documentation

### v1.0
- Initial release with basic functionality

## Authors

- M. Burhan Khan Lodhi
- Shomyal Khan
- Mushahid Zafar Jafri

**Institution**: University of Management and Technology, Sialkot

## License

All rights reserved by University of Management and Technology, Sialkot.

## Support

For technical support or bug reports, please contact the development team.