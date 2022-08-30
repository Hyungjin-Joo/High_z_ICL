import numpy as np
from astropy.io import fits

target = 'XMM1229'
filters = ['F105W','F160W']
sky0 = -0.6918
sky1 = -0.7228
cx, cy = 1760, 1670
size = 1000

ori0 = fits.open('%s_%s_drz_sci.fits'%(target,filters[0]))
ori1 = fits.open('%s_%s_drz_sci.fits'%(target,filters[1]))
sci0 = fits.open('%s_%s_mask.fits'%(target,filters[0]))
sci1 = fits.open('%s_%s_mask.fits'%(target,filters[1]))
sci0_max = fits.open('%s_%s_mask_max.fits'%(target,filters[0]))
sci1_max = fits.open('%s_%s_mask_max.fits'%(target,filters[1]))
detm = fits.open('detection_mask.fits')
deto = fits.open('detection.fits')
mask = fits.open('%s_seg_BCG_banned_2.0_add.fits'%target)
mosk = fits.open('%s_seg_BCG_banned_2.0.fits'%target)
mask_max = fits.open('%s_seg_BCG_banned_6.0_add.fits'%target)
mosk_max = fits.open('%s_seg_BCG_banned_6.0.fits'%target)
seg = fits.open('%s_seg.fits'%target)

o0_data = ori0[0].data - sky0
o1_data = ori1[0].data - sky1
s0_data = sci0[0].data
s0_data[np.isnan(s0_data)==True]=0
s1_data = sci1[0].data
s1_data[np.isnan(s1_data)==True]=0
s0_data_max = sci0_max[0].data
s0_data_max[np.isnan(s0_data_max)==True]=0
s1_data_max = sci1_max[0].data
s1_data_max[np.isnan(s1_data_max)==True]=0
dm_data = detm[0].data
do_data = deto[0].data
ma_data = mask[0].data
mo_data = mosk[0].data
ma_max_data = mask_max[0].data
mo_max_data = mosk_max[0].data
se_data = seg[0].data

samo0 = o0_data[cx-size:cx+size,cy-size:cy+size]
samo1 = o1_data[cx-size:cx+size,cy-size:cy+size]
samp0 = s0_data[cx-size:cx+size,cy-size:cy+size]
samp1 = s1_data[cx-size:cx+size,cy-size:cy+size]
samp0_max = s0_data_max[cx-size:cx+size,cy-size:cy+size]
samp1_max = s1_data_max[cx-size:cx+size,cy-size:cy+size]
sampd = dm_data[cx-size:cx+size,cy-size:cy+size]
sampdo = do_data[cx-size:cx+size,cy-size:cy+size]
sampm = ma_data[cx-size:cx+size,cy-size:cy+size]
sampmo = mo_data[cx-size:cx+size,cy-size:cy+size]
sampmma = ma_max_data[cx-size:cx+size,cy-size:cy+size]
sampmoma = mo_max_data[cx-size:cx+size,cy-size:cy+size]
sampse = se_data[cx-size:cx+size,cy-size:cy+size]

ori0[0].data = samo0
ori0.writeto('%s_%s_drz_sci_samp.fits'%(target,filters[0]), overwrite=True)
ori1[0].data = samo1
ori1.writeto('%s_%s_drz_sci_samp.fits'%(target,filters[1]), overwrite=True)
sci0[0].data = samp0
sci0.writeto('%s_%s_mask_samp.fits'%(target,filters[0]), overwrite=True)
sci1[0].data = samp1
sci1.writeto('%s_%s_mask_samp.fits'%(target,filters[1]), overwrite=True)
sci0_max[0].data = samp0_max
sci0_max.writeto('%s_%s_mask_max_samp.fits'%(target,filters[0]), overwrite=True)
sci1_max[0].data = samp1_max
sci1_max.writeto('%s_%s_mask_max_samp.fits'%(target,filters[1]), overwrite=True)
detm[0].data = sampd
detm.writeto('detection_mask_samp.fits', overwrite=True)
deto[0].data = sampdo
deto.writeto('detection_samp.fits', overwrite=True)
mask[0].data = sampm
mask.writeto('%s_seg_BCG_banned_2.0_add_samp.fits'%target, overwrite=True)
mosk[0].data = sampmo
mosk.writeto('%s_seg_BCG_banned_2.0_samp.fits'%target, overwrite=True)
mask_max[0].data = sampmma
mask_max.writeto('%s_seg_BCG_banned_6.0_add_samp.fits'%target, overwrite=True)
mosk_max[0].data = sampmoma
mosk_max.writeto('%s_seg_BCG_banned_6.0_samp.fits'%target, overwrite=True)
seg[0].data = sampse
seg.writeto('%s_seg_samp.fits'%target, overwrite=True)