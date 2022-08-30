import copy as cp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

target = 'XMM1229'
filters = ['F105W','F160W']
sky0 = -0.6918
sky1 = -0.7228
ske0 = (sky0 + 0.8406)
ske1 = (sky1 + 0.8330)

exs = [1515.213,1662.000,1507.100,1494.876,1483.263,1537.497,1540.210,1548.195,1554.165,1516.566,1626.726]
eys = [1471.061,1462.000,1442.387,1442.612,1443.775,1503.692,1497.878,1490.668,1484.311,1505.242,1384.462]
ers = [13.43914,10.29680,13.70895,17.66237,13.89291,9.331967,10.80362,12.49509,10.00892,6.838820,9.302715]

sci0 = fits.open('%s_%s_drz_sci.fits'%(target,filters[0]))
sci1 = fits.open('%s_%s_drz_sci.fits'%(target,filters[1]))
sci0_data = sci0[0].data
sci1_data = sci1[0].data
sci0_header = sci0[0].header
sci1_header = sci1[0].header
exptime = [sci0_header['EXPTIME'], sci1_header['EXPTIME']]
photplam = [sci0_header['PHOTPLAM'], sci1_header['PHOTPLAM']]
photflam = [sci0_header['PHOTFLAM'], sci1_header['PHOTFLAM']]

wht0 = fits.open('%s_%s_drz_wht.fits'%(target,filters[0]))
wht1 = fits.open('%s_%s_drz_wht.fits'%(target,filters[1]))
wht0_data = wht0[0].data
wht1_data = wht1[0].data

mask = fits.open('%s_seg_BCG_banned_2.0.fits'%target)
mask_data = mask[0].data

mask_max = fits.open('%s_seg_BCG_banned_6.0.fits'%target)
mask_max_data = mask_max[0].data

### sky correction
sci0_data -= sky0
sci1_data -= sky1

### masking
sci0_mask = cp.deepcopy(sci0_data)
sci1_mask = cp.deepcopy(sci1_data)
sci0_mask[mask_data!=0]=-32768
sci1_mask[mask_data!=0]=-32768

sci0_mask_max = cp.deepcopy(sci0_data)
sci1_mask_max = cp.deepcopy(sci1_data)
sci0_mask_max[mask_max_data!=0]=-32768
sci1_mask_max[mask_max_data!=0]=-32768

nx, ny = np.shape(sci0_data)
xx = np.linspace(0,nx,nx)
yy = np.linspace(0,ny,ny)
mx, my = np.meshgrid(yy, xx)

eye_mask = np.zeros((nx, ny))
for i in range(np.int64(np.size(np.array(exs)))):
    rr = np.sqrt((mx-exs[i])**2 + (my-eys[i])**2)
    bool = (rr<ers[i])
    eye_mask += bool
    
sci0_mask[eye_mask!=0]=-32768
sci1_mask[eye_mask!=0]=-32768

sci0_mask_max[eye_mask!=0]=-32768
sci1_mask_max[eye_mask!=0]=-32768

mask_data = mask_data + eye_mask

hdu = fits.PrimaryHDU(mask_data)
hdu.writeto('%s_seg_BCG_banned_2.0_add.fits'%target, overwrite=True)

mask_max_data = mask_max_data + eye_mask

hdu = fits.PrimaryHDU(mask_max_data)
hdu.writeto('%s_seg_BCG_banned_6.0_add.fits'%target, overwrite=True)

sci0[0].data = sci0_mask
sci0.writeto('%s_%s_mask2.fits'%(target,filters[0]), overwrite = True)
sci1[0].data = sci1_mask
sci1.writeto('%s_%s_mask2.fits'%(target,filters[1]), overwrite = True)

sci0[0].data = sci0_mask_max
sci0.writeto('%s_%s_mask_max.fits'%(target,filters[0]), overwrite = True)
sci1[0].data = sci1_mask_max
sci1.writeto('%s_%s_mask_max.fits'%(target,filters[1]), overwrite = True)

det_mask = (sci0_mask * wht0_data + sci1_mask * wht1_data) / (wht0_data + wht1_data)
hdu = fits.PrimaryHDU(det_mask)
hdu.writeto('detection_mask.fits', overwrite = True)
