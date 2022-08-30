import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clip

### Parameter setting
target = 'XMM1229'
redshift = 0.98
kpc_over_arcs = 7.969
arcs_over_pix = 0.05
kpc_over_pix = kpc_over_arcs * arcs_over_pix
filters = ['F105W','F160W']
x_center, y_center = 911, 743
ell, pa = 0.545, (132.828-90) * np.pi / 180
sky0 = -0.6918
sky1 = -0.7228
ske0 = (sky0 + 0.8406)
ske1 = (sky1 + 0.8330)

### Read files
sci0 = fits.open('%s_%s_mask_samp.fits'%(target,filters[0]))
sci1 = fits.open('%s_%s_mask_samp.fits'%(target,filters[1]))
sci0_data = sci0[0].data
sci0_data[sci0_data==0]=np.nan
sci1_data = sci1[0].data
sci1_data[sci1_data==0]=np.nan
sci0_head = sci0[0].header
sci1_head = sci1[0].header
photflam = [sci0_head['PHOTFLAM'], sci1_head['PHOTFLAM']]
photplam = [sci0_head['PHOTPLAM'], sci1_head['PHOTPLAM']]
exptime = [sci0_head['EXPTIME'], sci1_head['EXPTIME']]
print(ske0)
print(ske1)

### Ellipse binning
nx, ny = np.shape(sci0_data)
lx = np.linspace(0,nx,nx)
ly = np.linspace(0,ny,ny)
mx, my = np.meshgrid(lx, ly)
mx -= x_center
my -= y_center
ell = np.sqrt((mx * np.cos(pa) + my * np.sin(pa))**2 + (mx * np.sin(pa) - my * np.cos(pa))**2 / (1-ell)**2)

max_kpc = 500
max_pixel = max_kpc / kpc_over_pix
binsize = 25
initial_a = max_pixel ** (1/binsize)
ratio_bin = np.linspace(1,binsize,binsize)
sma_bin = initial_a**(ratio_bin)

### counts calculation
cnt0 = np.array([])
cnt1 = np.array([])
err0 = np.array([])
err1 = np.array([])

for i in range(binsize):
    if i == 0:
        bool_inner = np.ones((nx,ny))
    else:
        bool_inner = (ell > sma_bin[i-1])
    bool_outer = (ell <= sma_bin[i])
    bool_bin = bool_inner * bool_outer

    samp = sci0_data[bool_bin==1]
    bool = np.isnan(samp)
    samp = samp[bool==0]
    cnt0 = np.append(cnt0, np.nanmedian(samp))
    err0 = np.append(err0, np.nanstd(samp)/np.sqrt(np.size(samp)))
    
    samp = sci1_data[bool_bin==1]
    bool = np.isnan(samp)
    samp = samp[bool==0]
    cnt1 = np.append(cnt1, np.nanmedian(samp))
    err1 = np.append(err1, np.nanstd(samp)/np.sqrt(np.size(samp)))
    print(i, '/', binsize,'     ', end = '\r')



sma_bin *= kpc_over_pix
cnt0 /= arcs_over_pix**2
cnt1 /= arcs_over_pix**2
err0 /= arcs_over_pix**2
err1 /= arcs_over_pix**2
ske0 /= arcs_over_pix**2
ske1 /= arcs_over_pix**2
ter0 = np.sqrt(err0**2 + ske0**2)
ter1 = np.sqrt(err1**2 + ske1**2)


### save sfb
data = np.array([sma_bin.T, cnt0.T/exptime[0], cnt1.T/exptime[1], ter0/exptime[0], ter1/exptime[1]])
header = 'SMA[kpc], sfb[counts/s/arcs2], err[counts/s/arcs2]'
np.savetxt('%s_SFB_nature.txt'%target, data.T, header = header)
