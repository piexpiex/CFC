import sys
from astropy.io import fits
import numpy as np
#r=814
#x0=826
#y0=826
# 1 dentro y 0 fuera

fichero=sys.argv[1]
hdulist = fits.open(fichero)
data = hdulist[0].data
x=np.shape(data)[0]
y=np.shape(data)[1]

x0=int(sys.argv[2])
y0=int(sys.argv[3])
r=int(sys.argv[4])
mask=np.zeros((x,y))
for j in range(x):
	for k in range(y):
		if (j-x0)**2+(k-y0)**2<r**2:
			mask[j,k]=1

hdu = fits.PrimaryHDU(mask)
hdul = fits.HDUList([hdu])
hdul.writeto('custom_mask.fits')
