from scipy.optimize import curve_fit
import numpy as np
from numpy import fft
import matplotlib.pyplot as plt
from astropy.io import fits
from pymultinest.solve import solve
from pymultinest.analyse import Analyzer
from astropy.modeling.models import Sersic1D
from scipy.interpolate import interp1d
from scipy.signal import convolve2d

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

### Read obs files
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

### pre-process data
test0 = (sfb0 - sfe0) > 0
test1 = (sfb1 - sfe1) > 0
test_all = test0 * test1
sma = sma[test_all]
sfb0 = sfb0[test_all]
sfb1 = sfb1[test_all]
sfe0 = sfe0[test_all]
sfe1 = sfe1[test_all]

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

# psf0 = psf0[sma <= psf_r_0[-1]]
# psf1 = psf1[sma <= psf_r_1[-1]]

# ## 1D psf check
# plt.figure()
# plt.subplot(211)
# plt.plot(psf_r_0, psf_0)
# plt.scatter(bin_pixels, psf0, marker='x')
# plt.xscale('log')
# plt.yscale('log')
# plt.subplot(212)
# plt.plot(psf_r_1, psf_1)
# plt.scatter(bin_pixels, psf1, marker='x')
# plt.xscale('log')
# plt.yscale('log')
# plt.show()
    
### pymultinest functions
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
    x = np.linspace(1, np.size(a_mirror), np.size(a_mirror))
    b_mirror = np.concatenate((np.zeros(np.size(b)),b[::-1],b,np.zeros(np.size(b))))
    b_mirror /= np.sum(b) #normalize psf
    fft_a = fft.fft(a_mirror)
    fft_b = fft.fft(b_mirror, n = fft_a.size)
    fft_ab = fft_a*fft_b
    ab = fft.ifft(fft_ab) #shape : \__/
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
    cube[0] = cube[0] * 2000000
    cube[1] = cube[1] * 3 + 0.1
    cube[2] = cube[2] * 2.0 + 0.01
    cube[3] = cube[3] * 500000 + 10000
    cube[4] = cube[4] * 30 + 5
    cube[5] = cube[5] * 0.6 + 0.01
    cube[6] = cube[6] * 5000
    cube[7] = cube[7] * 300 + 10
    cube[8] = cube[8] * 1.0 + 0.01
    return cube

def loglikelihood_105(cube):
    global bin_kpc, sma, sfb0, sfe0, psf0
    m_sfb = triple_sersic_process(bin_kpc, sma, psf0, *cube)
    dsfb = sfb0 - m_sfb
    chi_sq = np.sum(dsfb**2 / sfe0**2)
    re_chi = - chi_sq / (np.size(sma) - 9)
    return re_chi

def loglikelihood_140(cube):
    global sma, sfb1, sfe1, psf1
    m_sfb = triple_sersic_process(bin_kpc, sma, psf1, *cube)
    dsfb = sfb1 - m_sfb
    chi_sq = np.sum(dsfb**2 / sfe1**2)
    re_chi = - chi_sq / (np.size(sma) - 9)
    return re_chi



### conv test
# ie1, re1, n1, ie2, re2, n2 = 260830.45085510967, 13.959419228517131, 0.9708206282581932, 2991.0182422660464, 179.32496820178275, 2.342955197000516
# m_sfb_test_nconv = double_sersic(ie1, re1, n1, ie2, re2, n2)
# sfb_test_nconv = m_sfb_test_nconv(sma)
# ds_m = m_sfb_test_nconv(bin_kpc)
# conv_m = conv(ds_m, psf0)
# sfb_test_conv = conv_interp(bin_kpc, conv_m, sma)
# plt.figure()
# plt.subplot(211)
# plt.scatter(sma, sfb0, marker = 'x', c = 'k', label = 'obs')
# plt.plot(sma, sfb_test_conv, c = 'b', label = 'ok conv')
# plt.plot(sma, sfb_test_nconv, c = 'r', label = 'no conv')
# plt.xscale('log')
# plt.yscale('log')
# # plt.ylim(np.amin(sfb0)/2, np.amax(sfb0)*2)
# plt.legend(loc=0)
# plt.subplot(212)
# plt.scatter(sma, sfb_test_conv / sfb_test_nconv)

# plt.xscale('log')
# plt.show()

### multinest directiory
import os
try: os.mkdir('%s_ts_conv'%filters[0])
except OSError: pass
try: os.mkdir('%s_ts_conv'%filters[1])
except OSError: pass

# number of dimensions our problem has
parameters_w = ["ie_a","re_a","n_a","ie_b","re_b","n_b","ie_c","re_c","n_c"]
n_params_w = len(parameters_w)

#######
#Filt1#
#######

# name of the output files
prefix_w = '%s_ts_conv/1-'%filters[0]

# run MultiNest
result_w = solve(LogLikelihood=loglikelihood_105, Prior=myprior, n_live_points = 2000, n_dims=n_params_w, outputfiles_basename=prefix_w, verbose=True)
print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result_w)
print()
err105 = []
print('parameter values:')
for name, col in zip(parameters_w, result_w['samples'].transpose()):
    print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))
    err105.append(col.std())
result = Analyzer(n_params_w, outputfiles_basename=prefix_w)
result_105 = result.get_best_fit()
params_105 = result_105['parameters']
ie_a105, re_a105, n_a105, ie_b105, re_b105, n_b105, ie_c105, re_c105, n_c105 = params_105

# ######
# Filt2#
# ######

# name of the output files
prefix_w = '%s_ts_conv/1-'%filters[1]

# run MultiNest
result_w = solve(LogLikelihood=loglikelihood_140, Prior=myprior, n_live_points = 2000, n_dims=n_params_w, outputfiles_basename=prefix_w, verbose=True)
print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result_w)
print()
print('parameter values:')
err140 = []
for name, col in zip(parameters_w, result_w['samples'].transpose()):
    print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))
    err140.append(col.std())
result = Analyzer(n_params_w, outputfiles_basename=prefix_w)
result_140 = result.get_best_fit()
params_140 = result_140['parameters']
ie_a140, re_a140, n_a140, ie_b140, re_b140, n_b140, ie_c140, re_c140, n_c140 = params_140


print(params_105,params_140)
skes = np.zeros(np.size(params_105))
skes[0] = ske0
skes[1] = ske1
data = np.array([params_105, params_140, err105, err140, skes])
np.savetxt(target+'_best_params.txt', data.T)

### plotting
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

plt.figure()
plt.subplot(1,2,1)
plt.errorbar(sma, mag0, yerr = [mag0_uerr, mag0_derr], fmt = 'o', markersize = 5, c = 'k', label = 'F105W')
plt.plot(sma, mag(triple_sersic(ie_a105, re_a105, n_a105, ie_b105, re_b105, n_b105, ie_c105, re_c105, n_c105)(sma),0), linestyle = '-', c = 'r')
plt.plot(sma, mag(sersic(ie_b105, re_b105, n_b105)(sma),0), linestyle = '-', c = 'b')
plt.plot(sma, mag(sersic(ie_a105, re_a105, n_a105)(sma),0), linestyle = '-', c = 'g')
plt.plot(sma, mag(sersic(ie_c105, re_c105, n_c105)(sma),0), linestyle = '-', c = 'orange')
plt.axhline(mag(ske0/0.05**2,0), c = 'k', linestyle = '--')
plt.axvline(5, c = 'orange', linestyle = '--')
plt.axvline(80, c = 'cyan', linestyle = '-.')
plt.xscale('log')
plt.xticks(fontsize = 'x-large')
plt.yticks(fontsize = 'x-large')
plt.xlabel('SMA [kpc]', fontsize = 'x-large')
plt.ylabel(r'$SB \  [mag \ arcs^{-2} ]$', fontsize = 'x-large')
plt.legend(loc=0, fontsize = 'x-large')
plt.ylim(31,18)

plt.subplot(1,2,2)
plt.errorbar(sma, mag1, yerr = [mag1_uerr, mag1_derr], fmt = 'o', markersize = 5, c = 'k', label = 'F140W')
plt.plot(sma, mag(triple_sersic(ie_a140, re_a140, n_a140, ie_b140, re_b140, n_b140, ie_c140, re_c140, n_c140)(sma),1), linestyle = '-', c = 'r')
plt.plot(sma, mag(sersic(ie_b140, re_b140, n_b140)(sma),1), linestyle = '-', c = 'b')
plt.plot(sma, mag(sersic(ie_a140, re_a140, n_a140)(sma),1), linestyle = '-', c = 'g')
plt.plot(sma, mag(sersic(ie_c140, re_c140, n_c140)(sma),0), linestyle = '-', c = 'orange')
plt.axhline(mag(ske1/0.05**2,1), c = 'k', linestyle = '--')
plt.axvline(5, c = 'orange', linestyle = '--')
plt.axvline(80, c = 'cyan', linestyle = '-.')
plt.xscale('log')
plt.xticks(fontsize = 'x-large')
plt.yticks(fontsize = 'x-large')
plt.xlabel('SMA [kpc]', fontsize = 'x-large')
plt.legend(loc=0, fontsize = 'x-large')
plt.ylim(31,18)
plt.yticks([])
plt.show()

plt.figure()
ax1 = plt.subplot(1,2,1)
plt.errorbar(sma, mag0, yerr = [mag0_uerr, mag0_derr], fmt = 'o', markersize = 5, c = 'k', label = 'F105W')
plt.plot(sma, mag(triple_sersic_process(bin_kpc, sma, psf0, ie_a105, re_a105, n_a105, ie_b105, re_b105, n_b105, ie_c105, re_c105, n_c105),0), linestyle = '-', c = 'r')
plt.plot(sma, mag(sersic_process(bin_kpc, sma, psf0, ie_a105, re_a105, n_a105) ,0), linestyle = '-', c = 'blue')
plt.plot(sma, mag(sersic_process(bin_kpc, sma, psf0, ie_b105, re_b105, n_b105) ,0), linestyle = '-', c = 'green')
plt.plot(sma, mag(sersic_process(bin_kpc, sma, psf0, ie_c105, re_c105, n_c105) ,0), linestyle = '-', c = 'orange')
plt.axhline(mag(ske0/0.05**2,0), c = 'k', linestyle = '--')
plt.axvline(5, c = 'orange', linestyle = '--')
plt.axvline(80, c = 'cyan', linestyle = '-.')
plt.xscale('log')
plt.xticks(fontsize = 'x-large')
plt.yticks(fontsize = 'x-large')
plt.xlabel('SMA [kpc]', fontsize = 'x-large')
plt.ylabel(r'$SB \  [mag \ arcs^{-2} ]$', fontsize = 'x-large')
plt.legend(loc=0, fontsize = 'x-large')
plt.ylim(31,18)

plt.subplot(1,2,2)
plt.errorbar(sma, mag1, yerr = [mag1_uerr, mag1_derr], fmt = 'o', markersize = 5, c = 'k', label = 'F140W')
plt.plot(sma, mag(triple_sersic_process(bin_kpc, sma, psf0, ie_a140, re_a140, n_a140, ie_b140, re_b140, n_b140, ie_c140, re_c140, n_c140),1), linestyle = '-', c = 'r')
plt.plot(sma, mag(sersic_process(bin_kpc, sma, psf0, ie_a140, re_a140, n_a140) ,1), linestyle = '-', c = 'blue')
plt.plot(sma, mag(sersic_process(bin_kpc, sma, psf0, ie_b140, re_b140, n_b140) ,1), linestyle = '-', c = 'green')
plt.plot(sma, mag(sersic_process(bin_kpc, sma, psf0, ie_c140, re_c140, n_c140) ,1), linestyle = '-', c = 'orange')
plt.axhline(mag(ske1/0.05**2,1), c = 'k', linestyle = '--')
plt.axvline(5, c = 'orange', linestyle = '--')
plt.axvline(80, c = 'cyan', linestyle = '-.')
plt.xscale('log')
plt.xticks(fontsize = 'x-large')
plt.yticks(fontsize = 'x-large')
plt.xlabel('SMA [kpc]', fontsize = 'x-large')
plt.legend(loc=0, fontsize = 'x-large')
plt.ylim(31,18)
plt.yticks([])
plt.show()
