from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from SNAP.Photometry import *

filename = "E489-1.Q0.B.170221_1012.A.045732.061920N2336.0060.nh.fits"
hdu = fits.open(filename)
image = hdu[0].data
hdu.close()

RA = [94.49705394, 94.49953134, 94.49971101, 94.49547835]
DEC = [-23.81448059, -23.8123183, -23.81819249, -23.81783911]
psftype = [2,3,2,3]
fitsky = [1,1,1,1]
names = ['sn', 'star1','star2','star3']

from SNAP.MagCalc import *

#load fits file, get relevant data
image, time, wcs = loadFits(filename, year=2016, getwcs=True, verbosity=1)

x0, y0 = np.zeros(len(RA)), np.zeros(len(DEC))
for i in range(len(x0)):
    x0[i], y0[i] = wcs.all_world2pix(RA[i], DEC[i], 0)
print x0[0], y0[0]

"""
intens, x, y = ap_multi(image, x0, y0, [1]*len(names), 0, 15)
plt.imshow(image, cmap='Greys', vmax=0.0001*np.amax(image), vmin=0)
plt.scatter(x, y, color='b', marker='.')
plt.show()

plt.imshow(image, cmap='Greys', vmax=0.0001*np.amax(image), vmin=0)
intens, x, y = ap_multi(image, x0, y0, fitsky, 35, 50)
plt.scatter(x, y, color='r', marker='.')
intens, x, y = ap_multi(image, x0, y0, fitsky, 50, 60)
plt.scatter(x, y, color='g', marker='.')
plt.scatter(x0, y0, color='b')
plt.show()
"""

RA, DEC, I, SN, M, Merr, Mlim = magnitude(image, image, wcs, 'aavso', 'E489-1.Q0.AAVSO.cat', (RA,DEC), radius=2000.0, psf=psftype, name=names, band='B', fwhm=5.0, limsnr=2.0, satmag=15.0, refmag=16.0, fitsky=fitsky, satpix=40000.0, verbosity=2)
#output position, magnitude
print time, RA, DEC, I, SN, M, Merr, Mlim
