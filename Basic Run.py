import numpy as np
from astropy.io import fits
from photutils import DAOStarFinder
from astropy.stats import mad_std
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture
def photometry(filepath):
 hdulist  = fits.open(filepath, ignore_missing_end=True)
 hdu = hdulist[0]
 hdu.data.shape
 image = hdu.data.astype(float)
 image -= np.median(image)
 bkg_sigma = mad_std(image)
 daofind = DAOStarFinder(fwhm=35., threshold=30.*bkg_sigma)
 sources = daofind(image)    #Save Stellar Sources from DOA Star Algorithm
 for col in sources.colnames:
     sources[col].info.format = '%.8g'  # for consistent table output
 print(sources)
 # Perform Aperture Photometry
 positions = (sources['xcentroid'], sources['ycentroid'])
 apertures = CircularAperture(positions, r=7.)
 phot_table = aperture_photometry(image, apertures)
 for col in phot_table.colnames:
     phot_table[col].info.format = '%.8g'  # for consistent table output
 filedest = input("Where to Save the result:")
 filedest1 = input("Where to Save the result 2:")
 print(phot_table)
 np.savetxt(filedest, (sources), delimiter=',')  # Save Result in CSV
 np.savetxt(filedest1, (phot_table), delimiter=',')  # Save Result in CSV
 plt.imshow(image, cmap='gray_r', origin='lower')
 apertures.plot(color='blue', lw=1.5, alpha=0.5)
 plt.show()

if __name__ == '__main__':
 filepath = input("Enter Your File Path:")
 photometry(filepath)