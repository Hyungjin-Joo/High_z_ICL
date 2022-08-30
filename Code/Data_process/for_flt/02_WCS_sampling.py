import numpy as np
import copy as cp
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.wcs import WCS

lim_s = input('WHT radius min? ')
lim_S = input('WHT radius max? ')
lim_min = np.int64(lim_s)
lim_max = np.int64(lim_S)

wht140 = fits.open('TEST_F105W_drz_wht.fits')
wdata = wht140[0].data
whead = wht140[0].header
wcs = WCS(whead)

nx, ny = wdata.shape
x = np.linspace(1,nx,nx)
y = np.linspace(1,ny,ny)
xv,yv = np.meshgrid(x,y)
rs = np.sqrt((xv-np.int64(ny/2))**2 + (yv-np.int64(nx/2))**2)

rr = np.linspace(lim_min, lim_max, 5)

R00 = (rs-rr[0])>0
R10 = (rs-rr[1])<=0
R20 = (rs-rr[1])>0
R30 = (rs-rr[2])<=0
R40 = (rs-rr[2])>0
R50 = (rs-rr[3])<=0
R60 = (rs-rr[3])>0
R70 = (rs-rr[4])<=0

B01 = (yv-np.int64(nx/2))>=0
B02 = (yv-np.int64(nx/2))<0
B03 = (xv-np.int64(ny/2))>=0
B04 = (xv-np.int64(ny/2))<0
B05 = (yv-np.int64(nx/2))<(xv-np.int64(ny/2))
B06 = (yv-np.int64(nx/2))>=(xv-np.int64(ny/2))
B07 = (yv-np.int64(nx/2))>=(np.int64(ny/2)-xv)
B08 = (yv-np.int64(nx/2))<(np.int64(ny/2)-xv)

aa = 1 * B01 * B05 + 2 * B06 * B03 + 3 * B04 * B07 + 4 * B01 * B08 + 5 * B02 * B06 + 6 * B05 * B04 + 7 * B03 * B08 + 8 * B07 * B02

sub_sky = R00 * R10 * aa + R20 * R30 * (aa + 8) + R40 * R50 * (aa + 16) + R60 * R70 * (aa + 24)

nx, ny = wdata.shape
xs = np.linspace(0,nx,nx)
ys = np.linspace(0,ny,ny)
mx, my = np.meshgrid(ys,xs)

coord = wcs.pixel_to_world([mx],[my])

files = np.loadtxt('plfs.lis', dtype = 'str', unpack = True)
skies = np.zeros(files.size)
for i in range(files.size):
    plf = fits.open(files[i])
    data_single = plf[1].data
    head_single = plf[1].header
    wcs_single = WCS(head_single)
    new_weight = np.zeros((data_single.shape))
    for j in range(32):
        skybin = (sub_sky.T==j+1)
        mx_bin = mx[skybin==1]
        my_bin = my[skybin==1]
        coord_bin = wcs.pixel_to_world([mx_bin],[my_bin])
        coord_single_tuple = wcs_single.world_to_pixel(coord_bin)
        coord_single = np.asarray(coord_single_tuple)
        ra_single = np.int64(coord_single[0,:])
        dec_single = np.int64(coord_single[1,:])
        bool1 = ra_single<1014
        bool2 = dec_single<1014
        bool3 = ra_single>0
        bool4 = dec_single>0
        bool = bool1 * bool2 * bool3 * bool4
        ra_single[bool == 0] = -32768
        dec_single[bool == 0] = -32768
        ra_single = ra_single[ra_single!=-32768]
        dec_single = dec_single[dec_single!=-32768]
        new_weight[ra_single[:], dec_single[:]] += j+1
    hdu = fits.PrimaryHDU(new_weight)
    split_fname = files[i].split('.')
    new_fname = '_wht.'.join(split_fname)
    hdu.writeto(new_fname, overwrite = True)
    print(i, '/', files.size -1, '            ', end = '\r')
sky_file = np.concatenate(([files], [skies], [skies]))
sky_file = sky_file.reshape(3,files.size)
np.savetxt('sky_level.txt',sky_file.T, fmt = '%s')
