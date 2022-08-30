import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
Target = input()
o105w = fits.open("%s_F105W_drz_wht.fits"%Target)
o105s = fits.open("%s_F105W_drz_sci.fits"%Target)
o140w = fits.open("%s_F160W_drz_wht.fits"%Target)
o140s = fits.open("%s_F160W_drz_sci.fits"%Target)

w105 = o105w[0].data
s105 = o105s[0].data
w140 = o140w[0].data
s140 = o140s[0].data

s105[np.isnan(s105)]=0
s140[np.isnan(s140)]=0

det = (w105*s105 + w140*s140)/(w105 + w140)

hdu = fits.PrimaryHDU(det)
hdu.writeto('detection.fits', overwrite = True)

wei = w105 + w140

hdu2 = fits.PrimaryHDU(wei)
hdu2.writeto('map_weight.fits', overwrite = True)
