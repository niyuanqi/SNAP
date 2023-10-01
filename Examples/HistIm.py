from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

filename = "N2292-1.Q1.B.161023_1631-161025_0722.XCXA.064318N2550.00005.00005.FM48.BS0512.coadd.REF.fits"
filename = "N2292-1.Q1.I.161023_1745-161025_0835.XCXA.064317N2550.00005.00005.FM37.BS0512.coadd.REF.fits"
filename = "N2292-1.Q1.V.161024_1637-161025_1638.XCXA.064318N2550.00005.00005.FM42.BS0512.coadd.REF.fits"

hdu = fits.open(filename)
image = hdu[0].data
hdu.close()

#hist, bins = np.histogram(image.flatten(), bins=100, range=(-3000, 80000))
hist, bins = np.histogram(image.flatten(), bins=100)
print hist, bins

plt.bar(bins[:-1], hist, align='edge', width=bins[1:]-bins[:-1])
plt.show()
