import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clip

Target = input()

### read fits
o105w = fits.open("%s_F105W_drz_wht.fits"%Target)
o105s = fits.open("%s_F105W_drz_sci.fits"%Target)
o140w = fits.open("%s_F160W_drz_wht.fits"%Target)
o140s = fits.open("%s_F160W_drz_sci.fits"%Target)

w105 = o105w[0].data
s105 = o105s[0].data
w140 = o140w[0].data
s140 = o140s[0].data

### wht process
samw105 = w105[w105>0]
m105 = np.median(samw105)
wh105 = w105/m105
wh105 = 1./wh105
wh105 = np.sqrt(wh105)

samw140 = w140[w140>0]
m140 = np.median(samw140)
wh140 = w140/m140
wh140 = 1./wh140
wh140 = np.sqrt(wh140)

### sci process
f105 = sigma_clip(s105, sigma=3, maxiters=None)
sigma105 = np.std(f105)

f140 = sigma_clip(s140, sigma=3, maxiters=None)
sigma140 = np.std(f140)

###rms making
rms105 = np.sqrt(s105 + sigma105**2)
rms105 = rms105*wh105
rms105[wh105<=0] = 1.e+18

rms140 = np.sqrt(s140 + sigma140**2)
rms140 = rms140*wh140
rms140[wh140<=0] = 1.e+18

### save fits
hdu1 = fits.PrimaryHDU(rms105)
hdu2 = fits.PrimaryHDU(rms140)
hdu1.writeto("%s_F105W_drz_sci_rms.fits"%Target, overwrite=True)
hdu2.writeto("%s_F160W_drz_sci_rms.fits"%Target, overwrite=True)
