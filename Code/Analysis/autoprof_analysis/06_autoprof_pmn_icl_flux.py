import copy as cp
from sys import meta_path
import numpy as np
from numpy import fft
import matplotlib.pyplot as plt
from astropy.io import fits
from pymultinest.solve import solve
from pymultinest.analyse import Analyzer
from astropy.modeling.models import Sersic1D
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d

### Parameter setting
target = 'SpARCS1049'
redshift = 1.71
kpc_over_arcs = 8.463
arcs_over_pix = 0.05
kpc_over_pix = kpc_over_arcs * arcs_over_pix
filters = ['F105W','F160W']
x_center, y_center = 983, 940
ell, pa = 0.051, (165.664-90) * np.pi / 180
sky0 = -5.928
sky1 = -6.017
ske0 = (sky0 + 6.400) * np.sqrt(32)
ske1 = (sky1 + 6.377) * np.sqrt(32)



### Read files
sci0 = fits.open('%s_%s_mask_samp.fits'%(target,filters[0]))
sci1 = fits.open('%s_%s_mask_samp.fits'%(target,filters[1]))
sci0_data = sci0[0].data
sci1_data = sci1[0].data
sci0_head = sci0[0].header
sci1_head = sci1[0].header
sma, sfb0, sfb1, sfe0, sfe1 = np.loadtxt('%s_SFB.txt'%target, unpack = True)
photflam = [sci0_head['PHOTFLAM'], sci1_head['PHOTFLAM']]
photplam = [sci0_head['PHOTPLAM'], sci1_head['PHOTPLAM']]
exptime = [sci0_head['EXPTIME'], sci1_head['EXPTIME']]
mesk = fits.open('%s_seg_mem_banned.fits'%target)
mem = mesk[0].data

### psf interpolation
psf_r_0, psf_0 = np.loadtxt(target+'_'+filters[0]+'_psf_extend.txt', unpack = True)
psf_r_1, psf_1 = np.loadtxt(target+'_'+filters[1]+'_psf_extend.txt', unpack = True)
# psf_r : kpcs
psf_f0 = interp1d(np.log10(psf_r_0), np.log10(psf_0), kind = 'linear', fill_value = 'extrapolate')
psf_f1 = interp1d(np.log10(psf_r_1), np.log10(psf_1), kind = 'linear', fill_value = 'extrapolate')

# linear binning of bins in pixels scale converted to kpc scale
bin_kpc = np.linspace(1,1500,1500) * kpc_over_pix

psf0 = 10**psf_f0(np.log10(bin_kpc))
psf1 = 10**psf_f1(np.log10(bin_kpc))

### counts/arcs > magnitude
def mag(flux, filter):
    abmag_zpt = -2.5 * np.log10(photflam[filter]) - 21.10 - 5 * np.log10(photplam[filter]) + 18.692
    mag = -2.5*np.log10(flux/exptime[filter]) + abmag_zpt
    return mag

def sersic(ie, re, n):
    ser = Sersic1D(amplitude = ie, r_eff = re, n = n)
    return ser

def double_sersic(ie_a, re_a, n_a, ie_b, re_b, n_b):
    return sersic(ie_a, re_a, n_a) + sersic(ie_b, re_b, n_b)

def triple_sersic(ie_a, re_a, n_a, ie_b, re_b, n_b, ie_c, re_c, n_c):
    return sersic(ie_a, re_a, n_a) + sersic(ie_b, re_b, n_b) + sersic(ie_c, re_c, n_c)

def conv(a,b):
    ### original shape: \, padding shape: _/\_ [zero, mirror, original, zero]
    a_mirror = np.concatenate((np.zeros(np.size(a)),a[::-1],a,np.zeros(np.size(a))))
    b_mirror = np.concatenate((np.zeros(np.size(b)),b[::-1],b,np.zeros(np.size(b))))
    b_mirror /= np.sum(b) #normalize psf
    fft_a = fft.fft(a_mirror)
    fft_b = fft.fft(b_mirror, n = fft_a.size)
    fft_ab = fft_a*fft_b
    ab = fft.ifft(fft_ab) #shape : \__/
    x =(np.linspace(0,np.size(ab),np.size(ab)))
    ab = ab[0:np.size(a)] #crop \ only
    return ab.real

def conv_interp(r, value, r_intp):
    model = interp1d(np.log10(r), np.log10(value), kind = 'linear', fill_value = 'extrapolate')
    v_intp = model(np.log10(r_intp))
    return 10**v_intp


def sersic_process(bin, sma, psf, ie_a, re_a, n_a,):
    os = sersic(ie_a, re_a, n_a)
    os_m = os(bin)
    conv_m = conv(os_m, psf)
    m_sfb = conv_interp(bin, conv_m, sma)
    return m_sfb

def double_sersic_process(bin, sma, psf, ie_a, re_a, n_a, ie_b, re_b, n_b):
    ds = double_sersic(ie_a, re_a, n_a, ie_b, re_b, n_b)
    ds_m = ds(bin)
    conv_m = conv(ds_m, psf)
    m_sfb = conv_interp(bin, conv_m, sma)
    return m_sfb

def triple_sersic_process(bin, sma, psf, ie_a, re_a, n_a, ie_b, re_b, n_b, ie_c, re_c, n_c):
    ts = triple_sersic(ie_a, re_a, n_a, ie_b, re_b, n_b, ie_c, re_c, n_c)
    ts_m = ts(bin)
    conv_m = conv(ts_m, psf)
    m_sfb = conv_interp(bin, conv_m, sma)
    return m_sfb

def myprior(cube):
    cube[0] = cube[0] * 1000000 + 100
    cube[1] = cube[1] * 10 + 0.1
    cube[2] = cube[2] * 1.5 + 0.1
    cube[3] = cube[3] * 10000 + 100
    cube[4] = cube[4] * 1000 + 15
    cube[5] = cube[5] * 0.8 + 0.2
    return cube

def loglikelihood_105(cube):
    global bin_kpc, sma, sfb0, sfe0, psf0
    m_sfb = double_sersic_process(bin_kpc, sma, psf0, *cube)
    dsfb = sfb0 - m_sfb
    chi_sq = np.sum(dsfb**2 / sfe0**2)
    re_chi = - chi_sq / (np.size(sma) - 6)
    return re_chi

def loglikelihood_140(cube):
    global sma, sfb1, sfe1, psf1
    m_sfb = double_sersic_process(bin_kpc, sma, psf1, *cube)
    dsfb = sfb1 - m_sfb
    chi_sq = np.sum(dsfb**2 / sfe1**2)
    re_chi = - chi_sq / (np.size(sma) - 6)
    return re_chi

### load MultiNest
#params_1 = np.loadtxt('%s_ds/1-.txt'%filters[0], usecols = (2,3,4,5,6,7))
params_0 = np.loadtxt('%s_ts_conv/1-phys_live.points'%filters[0], usecols = (0,1,2,3,4,5,6,7,8))

#params_1 = np.loadtxt('%s_ds/1-.txt'%filters[1], usecols = (2,3,4,5,6,7))
params_1 = np.loadtxt('%s_ts_conv/1-phys_live.points'%filters[1], usecols = (0,1,2,3,4,5,6,7,8))


sfb0_uerr = sfb0 + sfe0
sfb0_derr = sfb0 - sfe0
sfb1_uerr = sfb1 + sfe1
sfb1_derr = sfb1 - sfe1
mag0 = mag(sfb0, 0)
mag0_uerr = mag0 - mag(sfb0_uerr, 0)
mag0_derr = mag(sfb0_derr, 0) - mag0
mag1 = mag(sfb1, 1)
mag1_uerr = mag1 - mag(sfb1_uerr, 1)
mag1_derr = mag(sfb1_derr, 1) - mag1

nx, ny = np.shape(sci0_data)
lx = np.linspace(0,nx,nx)
ly = np.linspace(0,ny,ny)
mx, my = np.meshgrid(lx, ly)
mx -= x_center
my -= y_center
ell = np.sqrt((mx * np.cos(pa) + my * np.sin(pa))**2 + (mx * np.sin(pa) - my * np.cos(pa))**2 / (1-ell)**2) * kpc_over_pix

ICL_amounts = []
BCG_amounts = []
BCL_amounts = []
loop = np.int64(np.shape(params_0)[0])

plt.figure()
plt.subplot(1,2,1)
plt.errorbar(sma, mag0, yerr = [mag0_uerr, mag0_derr], fmt = 'o', markersize = 5, c = 'k', label = filters[0], zorder = loop+2)

for i in range(loop):
    ICL_comp = sersic_process(bin_kpc, sma, psf0, params_0[i,6], params_0[i,7], params_0[i,8])
    ICL_map = sersic_process(bin_kpc, ell, psf0, params_0[i,6], params_0[i,7], params_0[i,8]) * arcs_over_pix**2
    ICL_amounts.append(np.sum(ICL_map))
    plt.plot(sma, mag(ICL_comp,0), c = 'r', alpha = 0.005, zorder = i)
    BCG_comp = double_sersic_process(bin_kpc, sma, psf0, params_0[i,0], params_0[i,1], params_0[i,2], params_0[i,3], params_0[i,4], params_0[i,5])
    BCG_map = double_sersic_process(bin_kpc, ell, psf0, params_0[i,0], params_0[i,1], params_0[i,2], params_0[i,3], params_0[i,4], params_0[i,5]) * arcs_over_pix**2
    BCG_amounts.append(np.sum(BCG_map))
    plt.plot(sma, mag(BCG_comp,0), c = 'g', alpha = 0.005, zorder = i)
    BCL_comp = triple_sersic_process(bin_kpc, sma, psf0, params_0[i,0], params_0[i,1], params_0[i,2], params_0[i,3], params_0[i,4], params_0[i,5], params_0[i,6], params_0[i,7], params_0[i,8])
    BCL_map = triple_sersic_process(bin_kpc, ell, psf0, params_0[i,0], params_0[i,1], params_0[i,2], params_0[i,3], params_0[i,4], params_0[i,5], params_0[i,6], params_0[i,7], params_0[i,8]) * arcs_over_pix**2
    BCL_amounts.append(np.sum(BCL_map))
    plt.plot(sma, mag(BCL_comp,0), c = 'b', alpha = 0.005, zorder = i)
    print(i, '/', loop, end = '\r')
print('                             ')
plt.axhline(mag(ske0/0.05**2,0), c = 'k', linestyle = '--')
plt.xscale('log')
plt.xlabel('SMA [kpc]')
plt.ylabel(r'$SB [mag \ arcs^{-2} ]$')
plt.legend(loc=0)
plt.ylim(31,18)

ICL_amount_0 = np.median(ICL_amounts)
medid0 = np.argsort(ICL_amounts)[len(ICL_amounts)//2]
median_params_0 = params_0[medid0,:]
ICL_map = sersic_process(bin_kpc, ell, psf0, median_params_0[6], median_params_0[7], median_params_0[8]) * arcs_over_pix**2

hdu = fits.PrimaryHDU(ICL_map)
hdu.writeto(target+'_'+filters[0]+'_ICL_mask_conv.fits', overwrite = True)
kernel_0 = gaussian_kde(ICL_amounts, bw_method = 'silverman')
min_0, max_0 = np.amin(ICL_amounts)*0.9, np.amax(ICL_amounts)*1.1
x_0 = np.linspace(min_0, max_0, 1000)
plt.subplot(3,2,2)
n, bins, patches = plt.hist(np.array(ICL_amounts), bins = 50, density = True, histtype = 'step', color = 'black')
plt.ylim(0,np.amax(n)*1.2)
plt.plot(x_0, kernel_0(x_0), c = 'r')
plt.axvline(np.median(ICL_amounts), c = 'r', linestyle = ':')

dx = x_0[2] - x_0[1]
x = cp.deepcopy(ICL_amount_0)
while True:
    x += dx
    prob =  kernel_0.integrate_box_1d(ICL_amount_0, x)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        udx = cp.deepcopy(x)
        break
x = cp.deepcopy(ICL_amount_0)
while True:
    x -= dx
    prob =  kernel_0.integrate_box_1d(x, ICL_amount_0)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        ddx = cp.deepcopy(x)
        break
plt.axvline(udx, c = 'red', linestyle = '-.')
plt.axvline(ddx, c = 'red', linestyle = '-.')
print('ICL_counts: ',ICL_amount_0,'+',udx,'-',ddx)

BCG_amount_0 = np.median(BCG_amounts)
medid0 = np.argsort(BCG_amounts)[len(BCG_amounts)//2]
median_params_0 = params_0[medid0,:]
BCG_map = double_sersic_process(bin_kpc, ell, psf0, median_params_0[0], median_params_0[1], median_params_0[2], median_params_0[3], median_params_0[4], median_params_0[5]) * arcs_over_pix**2

hdu = fits.PrimaryHDU(BCG_map)
hdu.writeto(target+'_'+filters[0]+'_BCG_mask_conv.fits', overwrite = True)
kernel_0 = gaussian_kde(BCG_amounts, bw_method = 'silverman')
min_0, max_0 = np.amin(BCG_amounts)*0.9, np.amax(BCG_amounts)*1.1
x_0 = np.linspace(min_0, max_0, 1000)
plt.subplot(3,2,4)
n, bins, patches = plt.hist(np.array(BCG_amounts), bins = 50, density = True, histtype = 'step', color = 'black')
plt.ylim(0,np.amax(n)*1.2)
plt.plot(x_0, kernel_0(x_0), c = 'r')
plt.axvline(np.median(BCG_amounts), c = 'r', linestyle = ':')

dx = x_0[2] - x_0[1]
x = cp.deepcopy(BCG_amount_0)
while True:
    x += dx
    prob =  kernel_0.integrate_box_1d(BCG_amount_0, x)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        udx = cp.deepcopy(x)
        break
x = cp.deepcopy(BCG_amount_0)
while True:
    x -= dx
    prob =  kernel_0.integrate_box_1d(x, BCG_amount_0)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        ddx = cp.deepcopy(x)
        break
plt.axvline(udx, c = 'red', linestyle = '-.')
plt.axvline(ddx, c = 'red', linestyle = '-.')
print('BCG_counts: ',BCG_amount_0,'+',udx,'-',ddx)

BCL_amount_0 = np.median(BCL_amounts)
medid0 = np.argsort(BCL_amounts)[len(BCL_amounts)//2]
median_params_0 = params_0[medid0,:]
BCL_map = triple_sersic_process(bin_kpc, ell, psf0, median_params_0[0], median_params_0[1], median_params_0[2], median_params_0[3], median_params_0[4], median_params_0[5], median_params_0[6], median_params_0[7], median_params_0[8]) * arcs_over_pix**2

hdu = fits.PrimaryHDU(BCL_map)
hdu.writeto(target+'_'+filters[0]+'_BCL_mask_conv.fits', overwrite = True)
kernel_0 = gaussian_kde(BCL_amounts, bw_method = 'silverman')
min_0, max_0 = np.amin(BCL_amounts)*0.9, np.amax(BCL_amounts)*1.1
x_0 = np.linspace(min_0, max_0, 1000)
plt.subplot(3,2,6)
n, bins, patches = plt.hist(np.array(BCL_amounts), bins = 50, density = True, histtype = 'step', color = 'black')
plt.ylim(0,np.amax(n)*1.2)
plt.plot(x_0, kernel_0(x_0), c = 'r')
plt.axvline(np.median(BCL_amounts), c = 'r', linestyle = ':')

dx = x_0[2] - x_0[1]
x = cp.deepcopy(BCL_amount_0)
while True:
    x += dx
    prob =  kernel_0.integrate_box_1d(BCL_amount_0, x)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        udx = cp.deepcopy(x)
        break
x = cp.deepcopy(BCL_amount_0)
while True:
    x -= dx
    prob =  kernel_0.integrate_box_1d(x, BCL_amount_0)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        ddx = cp.deepcopy(x)
        break
plt.axvline(udx, c = 'red', linestyle = '-.')
plt.axvline(ddx, c = 'red', linestyle = '-.')
plt.suptitle(target+', '+filters[0])
plt.savefig(target+'_'+filters[0]+'fraction_error.png')
plt.close()
print('BCL_counts: ',BCL_amount_0,'+',udx,'-',ddx)



plt.figure()
plt.subplot(331)
plt.hist(BCG_amounts)
plt.subplot(334)
plt.scatter(BCG_amounts, ICL_amounts)
plt.subplot(335)
plt.hist(ICL_amounts)
plt.subplot(337)
plt.scatter(BCG_amounts, BCL_amounts)
plt.subplot(338)
plt.scatter(ICL_amounts, BCL_amounts)
plt.subplot(339)
plt.hist(BCL_amounts)
plt.suptitle(target+'/'+filters[0])
plt.savefig(target+'_'+filters[0]+'_covariance.png')
plt.close()

amounts = np.array([BCG_amounts, ICL_amounts, BCG_amounts])
np.savetxt(target+'_'+filters[0]+'_covariance.txt', amounts.T)

############## Filter 0 done. Filter 1 start.


ICL_amounts = []
BCG_amounts = []
BCL_amounts = []
loop = np.int64(np.shape(params_1)[0])

plt.figure()
plt.subplot(1,2,1)
plt.errorbar(sma, mag0, yerr = [mag0_uerr, mag0_derr], fmt = 'o', markersize = 5, c = 'k', label = filters[1], zorder = loop+2)

for i in range(loop):
    ICL_comp = sersic_process(bin_kpc, sma, psf1, params_1[i,6], params_1[i,7], params_1[i,8])
    ICL_map = sersic_process(bin_kpc, ell, psf1, params_1[i,6], params_1[i,7], params_1[i,8]) * arcs_over_pix**2
    ICL_amounts.append(np.sum(ICL_map))
    plt.plot(sma, mag(ICL_comp,0), c = 'r', alpha = 0.005, zorder = i)
    BCG_comp = double_sersic_process(bin_kpc, sma, psf1, params_1[i,0], params_1[i,1], params_1[i,2], params_1[i,3], params_1[i,4], params_1[i,5])
    BCG_map = double_sersic_process(bin_kpc, ell, psf1, params_1[i,0], params_1[i,1], params_1[i,2], params_1[i,3], params_1[i,4], params_1[i,5]) * arcs_over_pix**2
    BCG_amounts.append(np.sum(BCG_map))
    plt.plot(sma, mag(BCG_comp,0), c = 'g', alpha = 0.005, zorder = i)
    BCL_comp = triple_sersic_process(bin_kpc, sma, psf1, params_1[i,0], params_1[i,1], params_1[i,2], params_1[i,3], params_1[i,4], params_1[i,5], params_1[i,6], params_1[i,7], params_1[i,8])
    BCL_map = triple_sersic_process(bin_kpc, ell, psf1, params_1[i,0], params_1[i,1], params_1[i,2], params_1[i,3], params_1[i,4], params_1[i,5], params_1[i,6], params_1[i,7], params_1[i,8]) * arcs_over_pix**2
    BCL_amounts.append(np.sum(BCL_map))
    plt.plot(sma, mag(BCL_comp,0), c = 'b', alpha = 0.005, zorder = i)
    print(i, '/', loop, end = '\r')
print('                             ')
plt.axhline(mag(ske0/0.05**2,0), c = 'k', linestyle = '--')
plt.xscale('log')
plt.xlabel('SMA [kpc]')
plt.ylabel(r'$SB [mag \ arcs^{-2} ]$')
plt.legend(loc=0)
plt.ylim(31,18)

ICL_amount_1 = np.median(ICL_amounts)
medid0 = np.argsort(ICL_amounts)[len(ICL_amounts)//2]
median_params_1 = params_1[medid0,:]
ICL_map = sersic_process(bin_kpc, ell, psf1, median_params_1[6], median_params_1[7], median_params_1[8]) * arcs_over_pix**2

hdu = fits.PrimaryHDU(ICL_map)
hdu.writeto(target+'_'+filters[1]+'_ICL_mask_conv.fits', overwrite = True)
kernel_1 = gaussian_kde(ICL_amounts, bw_method = 'silverman')
min_1, max_1 = np.amin(ICL_amounts)*0.9, np.amax(ICL_amounts)*1.1
x_1 = np.linspace(min_1, max_1, 1000)
plt.subplot(3,2,2)
n, bins, patches = plt.hist(np.array(ICL_amounts), bins = 50, density = True, histtype = 'step', color = 'black')
plt.ylim(0,np.amax(n)*1.2)
plt.plot(x_1, kernel_1(x_1), c = 'r')
plt.axvline(np.median(ICL_amounts), c = 'r', linestyle = ':')

dx = x_1[2] - x_1[1]
x = cp.deepcopy(ICL_amount_1)
while True:
    x += dx
    prob =  kernel_1.integrate_box_1d(ICL_amount_1, x)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        udx = cp.deepcopy(x)
        break
x = cp.deepcopy(ICL_amount_1)
while True:
    x -= dx
    prob =  kernel_1.integrate_box_1d(x, ICL_amount_1)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        ddx = cp.deepcopy(x)
        break
plt.axvline(udx, c = 'red', linestyle = '-.')
plt.axvline(ddx, c = 'red', linestyle = '-.')
print('ICL_counts: ',ICL_amount_1,'+',udx,'-',ddx)

BCG_amount_1 = np.median(BCG_amounts)
medid0 = np.argsort(BCG_amounts)[len(BCG_amounts)//2]
median_params_1 = params_1[medid0,:]
BCG_map = double_sersic_process(bin_kpc, ell, psf1, median_params_1[0], median_params_1[1], median_params_1[2], median_params_1[3], median_params_1[4], median_params_1[5]) * arcs_over_pix**2

hdu = fits.PrimaryHDU(BCG_map)
hdu.writeto(target+'_'+filters[1]+'_BCG_mask_conv.fits', overwrite = True)
kernel_1 = gaussian_kde(BCG_amounts, bw_method = 'silverman')
min_1, max_1 = np.amin(BCG_amounts)*0.9, np.amax(BCG_amounts)*1.1
x_1 = np.linspace(min_1, max_1, 1000)
plt.subplot(3,2,4)
n, bins, patches = plt.hist(np.array(BCG_amounts), bins = 50, density = True, histtype = 'step', color = 'black')
plt.ylim(0,np.amax(n)*1.2)
plt.plot(x_1, kernel_1(x_1), c = 'r')
plt.axvline(np.median(BCG_amounts), c = 'r', linestyle = ':')

dx = x_1[2] - x_1[1]
x = cp.deepcopy(BCG_amount_1)
while True:
    x += dx
    prob =  kernel_1.integrate_box_1d(BCG_amount_1, x)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        udx = cp.deepcopy(x)
        break
x = cp.deepcopy(BCG_amount_1)
while True:
    x -= dx
    prob =  kernel_1.integrate_box_1d(x, BCG_amount_1)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        ddx = cp.deepcopy(x)
        break
plt.axvline(udx, c = 'red', linestyle = '-.')
plt.axvline(ddx, c = 'red', linestyle = '-.')
print('BCG_counts: ',BCG_amount_1,'+',udx,'-',ddx)

BCL_amount_1 = np.median(BCL_amounts)
medid0 = np.argsort(BCL_amounts)[len(BCL_amounts)//2]
median_params_1 = params_1[medid0,:]
BCL_map = triple_sersic_process(bin_kpc, ell, psf1, median_params_1[0], median_params_1[1], median_params_1[2], median_params_1[3], median_params_1[4], median_params_1[5], median_params_1[6], median_params_1[7], median_params_1[8]) * arcs_over_pix**2

hdu = fits.PrimaryHDU(BCL_map)
hdu.writeto(target+'_'+filters[1]+'_BCL_mask_conv.fits', overwrite = True)
kernel_1 = gaussian_kde(BCL_amounts, bw_method = 'silverman')
min_1, max_1 = np.amin(BCL_amounts)*0.9, np.amax(BCL_amounts)*1.1
x_1 = np.linspace(min_1, max_1, 1000)
plt.subplot(3,2,6)
n, bins, patches = plt.hist(np.array(BCL_amounts), bins = 50, density = True, histtype = 'step', color = 'black')
plt.ylim(0,np.amax(n)*1.2)
plt.plot(x_1, kernel_1(x_1), c = 'r')
plt.axvline(np.median(BCL_amounts), c = 'r', linestyle = ':')

dx = x_1[2] - x_1[1]
x = cp.deepcopy(BCL_amount_1)
while True:
    x += dx
    prob =  kernel_1.integrate_box_1d(BCL_amount_1, x)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        udx = cp.deepcopy(x)
        break
x = cp.deepcopy(BCL_amount_1)
while True:
    x -= dx
    prob =  kernel_1.integrate_box_1d(x, BCL_amount_1)
    print(prob, '   ', end = '\r')
    if prob >= 0.341:
        ddx = cp.deepcopy(x)
        break
plt.axvline(udx, c = 'red', linestyle = '-.')
plt.axvline(ddx, c = 'red', linestyle = '-.')
plt.suptitle(target+', '+filters[1])
plt.savefig(target+'_'+filters[1]+'fraction_error.png')
plt.close()
print('BCL_counts: ',BCL_amount_1,'+',udx,'-',ddx)

plt.figure()
plt.subplot(331)
plt.hist(BCG_amounts)
plt.subplot(334)
plt.scatter(BCG_amounts, ICL_amounts)
plt.subplot(335)
plt.hist(ICL_amounts)
plt.subplot(337)
plt.scatter(BCG_amounts, BCL_amounts)
plt.subplot(338)
plt.scatter(ICL_amounts, BCL_amounts)
plt.subplot(339)
plt.hist(BCL_amounts)
plt.suptitle(target+'/'+filters[1])
plt.savefig(target+'_'+filters[1]+'_covariance.png')
plt.close()

amounts = np.array([BCG_amounts, ICL_amounts, BCG_amounts])
np.savetxt(target+'_'+filters[1]+'_covariance.txt', amounts.T)