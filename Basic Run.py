import numpy as np
from astropy.io import fits
from photutils import DAOStarFinder
from astropy.stats import mad_std
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture
from PIL import Image
from rawkit.raw import Raw
import rawpy.enhance
from cr2fits import cr2fits

class Photometry:
    def photometry(self,filepath, b, d):
        hdulist = fits.open(filepath, ignore_missing_end=True)
        hdu = hdulist[0]
        hdu.data.shape
        image = hdu.data.astype(float)
        image -= np.median(image)
        bkg_sigma = mad_std(image)
        daofind = DAOStarFinder(fwhm=b, threshold=d * bkg_sigma)
        sources = daofind(image)  # Save Stellar Sources from DOA Star Algorithm
        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output
        print(sources)
        # Perform Aperture Photometry
        positions = (sources['xcentroid'], sources['ycentroid'])
        apertures = CircularAperture(positions, r=17.)
        phot_table = aperture_photometry(image, apertures)
        for col in phot_table.colnames:
            phot_table[col].info.format = '%.8g'  # for consistent table output
        filedest1 = input("Where to Save the result:")
        filedest2 = input("Where to Save the result 2:")
        print(phot_table)
        np.savetxt(filedest1, (sources), delimiter=',')  # Save Result in CSV
        np.savetxt(filedest2, (phot_table), delimiter=',')  # Save Result in CSV
        plt.imshow(image, cmap='gray_r', origin='lower')
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        plt.show()

class Conversion:
    def Raw_Fits(self, raw_path, save_location):
        paths = [raw_path]
        bad_pixels = rawpy.enhance.find_bad_pixels(paths)

        for path in paths:
            with rawpy.imread(path) as raw:
                rawpy.enhance.repair_bad_pixels(raw, bad_pixels, method='median')
                rgb = raw.postprocess(no_auto_bright=True, use_auto_wb=False, gamma=None)
        a = np.array(rgb)
        print(a.shape)

        filename = raw_path
        raw_image = Raw(filename)
        buffered_image = np.array(raw_image.to_buffer())

        image = Image.frombytes('RGB', (raw_image.metadata.width, raw_image.metadata.height), a).convert('LA')
        xsize, ysize = image.size
        data1 = np.array(image.getdata())
        print(data1.shape)

        r = [(d[0]) for d in data1]
        g = [(d[1]) for d in data1]
        r_1 = np.array(r)
        g_1 = np.array(g)

        r_data = np.array(r_1.data)
        g_data = np.array(g_1.data)
        print(r_data.shape)

        r_data = r_data.reshape(ysize, xsize)
        g_data = g_data.reshape(ysize, xsize)

        a = cr2fits(raw_path, 0)
        b = cr2fits.read_exif(a)

        concat = r_data + g_data
        hdu = fits.PrimaryHDU(data=concat)
        hdu.header.set('OBSTIME', a.date)
        hdu.header.set('OBSTIME', a.date)
        hdu.header.set('EXPTIME', a.shutter)
        hdu.header.set('APERTUR', a.aperture)
        hdu.header.set('ISO', a.iso)
        hdu.header.set('FOCAL', a.focal)
        hdu.header.set('ORIGIN', a.original_file)
        hdu.header.set('FILTER', a.colors[a.colorInput])
        hdu.header.set('CAMERA', a.camera)
        hdu.writeto(save_location, overwrite=True)

if __name__ == '__main__':
    P1 = Photometry()
    C1 = Conversion()
    print("Menu")
    print("1. Convert Raw file into Fits")
    print("2. Perform Aperture Photometry")
    choice = input("Enter your choice: ")

    if choice == '1':
        raw_path = input("Enter Raw file path to convert: ")
        save_location = input("Enter Destination file location ")
        C1.Raw_Fits(raw_path, save_location)

    elif choice == '2':
        filepath = input("Enter Your File Path:")
        a = input("Enter Base Value:")
        b = int(a)
        c = input("Enter threshold value:")
        d = int(c)
        P1.photometry(filepath, b, d)

    else:
        print("Invalid Command")
