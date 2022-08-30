import copy as cp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde

### Parameter setting
target = 'MOO1142'
redshift = 1.19
kpc_over_arcs = 8.282
arcs_over_pix = 0.05
kpc_over_pix = kpc_over_arcs * arcs_over_pix
filters = ['F105W','F140W']
x_center, y_center = 583, 810
ell, pa = 0.535, (84.151-90) * np.pi / 180
sky0 = -3.1978
sky1 = -4.6000
ske0 = (sky0 + 3.5530) * np.sqrt(32)
ske1 = (sky1 + 5.280) * np.sqrt(32)

### Read files
ori0 = fits.open('%s_%s_drz_sci_samp.fits'%(target,filters[0]))
ori1 = fits.open('%s_%s_drz_sci_samp.fits'%(target,filters[1]))
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
mask = fits.open('%s_seg_BCG_banned_2.0_samp.fits'%target)
seg = mask[0].data
mask_max = fits.open('%s_seg_BCG_banned_6.0_samp.fits'%target)
seg_max = mask_max[0].data

#### BCG & spec-mem
num_bcg = 1318
num_spe = np.array([2125, 3185, 2263, 1577, 2375])
num_non = np.array([5629, 5401, 4389, 2735, 385, 419])

num, iso0, iso0_err, auto0, auto0_err, xi, yi = np.loadtxt('%s_%s_drz_sci.cat'%(target,filters[0]), usecols = (0,1,2,10,11,20,21), unpack = True)
iso1, iso1_err, auto1, auto1_err, cs = np.loadtxt('%s_%s_drz_sci.cat'%(target,filters[1]), usecols = (1,2,10,11,34), unpack = True)

def mag(flux, filter):
    abmag_zpt = -2.5 * np.log10(photflam[filter]) - 21.10 - 5 * np.log10(photplam[filter]) + 18.692
    mag = -2.5*np.log10(flux/exptime[filter]) + abmag_zpt
    return mag

def linear(x, a, b):
    return(x*a + b)

mag0_iso = mag(iso0,0)
mag1_iso = mag(iso1,1)
mag0_auto = mag(auto0,0)
mag1_auto = mag(auto1,1)

mag0_iso_u = mag0_iso - mag(iso0+iso0_err,0)
mag1_iso_u = mag1_iso - mag(iso1+iso1_err,1)
mag0_auto_u = mag0_auto - mag(auto0+auto0_err,0)
mag1_auto_u = mag1_auto - mag(auto1+auto1_err,1)

mag0_iso_d = mag(iso0-iso0_err,0) - mag0_iso
mag1_iso_d = mag(iso1-iso1_err,1) - mag1_iso
mag0_auto_d = mag(auto0-auto0_err,0) - mag0_auto
mag1_auto_d = mag(auto1-auto1_err,1) - mag1_auto

bool0iso = np.isnan(mag0_iso)
mag0_iso[bool0iso==True] = 0
bool1iso = np.isnan(mag1_iso)
mag1_iso[bool1iso==True] = 0
bool0auto = np.isnan(mag0_auto)
mag0_auto[bool0auto==True] = 0
bool1auto = np.isnan(mag1_auto)
mag1_auto[bool1auto==True] = 0

mag0_iso_spe = mag0_iso[num_spe-1]
mag1_iso_spe = mag1_iso[num_spe-1]
mag0_auto_spe = mag0_auto[num_spe-1]
mag1_auto_spe = mag1_auto[num_spe-1]

mag0_iso_bcg = mag0_iso[num_bcg-1]
mag1_iso_bcg = mag1_iso[num_bcg-1]
mag0_auto_bcg = mag0_auto[num_bcg-1]
mag1_auto_bcg = mag1_auto[num_bcg-1]

bool_cs = cs<=0.4

x_cmd = np.linspace(18,26,100)
popt, pcov = curve_fit(linear, mag0_iso_spe, mag0_auto_spe-mag1_auto_spe)
popt[0] = -0.01
popt[1] = 1
y_cmd = linear(x_cmd, *popt)

plt.figure()
plt.scatter(mag0_iso, mag0_auto - mag1_auto, c = 'k', s = 5)
plt.scatter(mag0_iso, mag0_auto - mag1_auto, c = 'k', s = 5)
plt.scatter(mag0_iso_spe, mag0_auto_spe - mag1_auto_spe, c = 'orange', s = 5)
plt.scatter(mag0_iso_bcg, mag0_auto_bcg - mag1_auto_bcg, c = 'r', s = 5)
plt.plot(x_cmd, y_cmd, linestyle = '--', c = 'r')
plt.xlim(18,26)
plt.ylim(-2,2)
plt.show()

bool = cp.deepcopy(bool_cs)
step = 0
while True:
    print(step, popt)
    num_samp = num[bool==True]
    mag0_iso_samp = mag0_iso[bool==True]
    mag0_iso_u_samp = mag0_iso_u[bool==True]
    mag0_iso_d_samp = mag0_iso_d[bool==True]
    mag1_iso_samp = mag1_iso[bool==True]
    mag1_iso_u_samp = mag1_iso_u[bool==True]
    mag1_iso_d_samp = mag1_iso_d[bool==True]
    mag0_auto_samp = mag0_auto[bool==True]
    mag0_auto_u_samp = mag0_auto_u[bool==True]
    mag0_auto_d_samp = mag0_auto_d[bool==True]
    mag1_auto_samp = mag1_auto[bool==True]
    mag1_auto_u_samp = mag1_auto_u[bool==True]
    mag1_auto_d_samp = mag1_auto_d[bool==True]
    color_samp = mag0_auto_samp - mag1_auto_samp
    criteria = linear(mag0_iso_samp, *popt)
    bool_u = criteria-color_samp >= -0.125
    bool_d = criteria-color_samp <= 0.125
    bool_l = mag0_iso_samp >= mag0_iso_bcg
    bool_r = mag0_iso_samp <= 26
    
    bool_ud = bool_u * bool_d * bool_l * bool_r
    num_plot = num_samp[bool_ud==True]
    mag0_iso_plot = mag0_iso_samp[bool_ud==True]
    mag1_iso_plot = mag1_iso_samp[bool_ud==True]
    mag0_auto_plot = mag0_auto_samp[bool_ud==True]
    mag1_auto_plot = mag1_auto_samp[bool_ud==True]

    num_plot = np.append(num_plot, num_spe, axis = None)
    mag0_iso_plot = np.append(mag0_iso_plot, mag0_iso_spe, axis = None)
    mag1_iso_plot = np.append(mag1_iso_plot, mag1_iso_spe, axis = None)
    mag0_auto_plot = np.append(mag0_auto_plot, mag0_auto_spe, axis = None)
    mag1_auto_plot = np.append(mag1_auto_plot, mag1_auto_spe, axis = None)


#    plt.figure()
#    plt.scatter(mag0_iso_samp, color_samp, c = 'k', s = 5)
#    plt.scatter(mag0_iso_plot, mag0_auto_plot - mag1_auto_plot, c = 'orange', s = 5)
#    plt.scatter(mag0_iso_spe, mag0_auto_spe - mag1_auto_spe, c = 'r', s = 5)
#    plt.plot(x_cmd, linear(x_cmd, *popt), c = 'r', linestyle = '--')
    r_popt, r_pcov = curve_fit(linear, mag0_iso_plot, mag0_auto_plot - mag1_auto_plot)
#    plt.plot(x_cmd, linear(x_cmd, *r_popt), c = 'b', linestyle = '--')
#    plt.xlim(18,26)
#    plt.ylim(-2,2)
#    plt.show()
    step += 1
#    print(np.sqrt(np.sum((popt-r_popt)**2))
    if np.sqrt(np.sum((popt-r_popt)**2)) < 0.001:
        break
    else:
        popt, pcov = r_popt, r_pcov
    
test = (mag0_auto_samp - mag1_auto_samp) - (popt[0] * mag0_iso_samp + popt[1] )

bool_plot_1_1 = mag0_iso_samp >= mag0_iso_bcg
bool_plot_1_2 = mag1_iso_samp >= mag1_iso_bcg
bool_plot_1_3 = mag0_iso_samp <= 26
bool_plot_1_4 = mag1_iso_samp <= 26
bool_plot_1_5 = test >= -2
bool_plot_1_6 = test <= 2
bool_plot_1 = bool_plot_1_1 * bool_plot_1_2  * bool_plot_1_3 * bool_plot_1_4 * bool_plot_1_5 * bool_plot_1_6

test_median = np.nanmedian(test[bool_plot_1!=0])
test_std = np.nanstd(test[bool_plot_1!=0], ddof = 1)/np.sqrt(np.size(test[bool_plot_1!=0] - 1))

kernel = gaussian_kde(test[bool_plot_1!=0])
ran = np.linspace(-2,2,10000)
pdf = kernel.pdf(ran)

def double_gauss(x, a1, m1, s1, a2, m2, s2):
    y1 = a1 * (1 / np.sqrt(2 * np.pi * s1**2)) * np.exp(-(x-m1)**2 / (2 * s1**2))
    y2 = a2 * (1 / np.sqrt(2 * np.pi * s2**2)) * np.exp(-(x-m2)**2 / (2 * s2**2))
    return y1 + y2

g_popt, g_pcov = curve_fit(double_gauss, ran, pdf, p0 = [1.0, -0.5, 0.1, 0.5, 0, 0.1])
test_median = g_popt[4]
crit = g_popt[1] + g_popt[2] 
samples = test[bool_plot_1!=0]
test_std = g_popt[5]# / np.sqrt(np.sum(samples >= crit)-1)



plt.figure()
plt.hist(test[bool_plot_1!=0], bins = 100, density = True)
plt.plot(ran,pdf)
plt.plot(ran, double_gauss(ran, *g_popt))
plt.plot(ran, g_popt[3] * (1 / np.sqrt(2 * np.pi * g_popt[5]**2)) * np.exp(-(ran-g_popt[4])**2 / (2 * g_popt[5]**2)))
plt.axvline(test_median, c = 'k')
plt.axvline(test_median + test_std, c = 'k')
plt.axvline(test_median - test_std, c = 'k')
plt.show()


print(test_median)
print(test_std)
bool_plot_2 = test + np.sqrt(mag0_auto_u_samp**2 + mag1_auto_u_samp**2) >= - test_std
bool_plot_3 = test - np.sqrt(mag0_auto_d_samp**2 + mag1_auto_d_samp**2) <= + test_std

bool_non = np.zeros(bool_plot_1.size)
for i in range(num_non.size):
    bool_non+= num_samp==num_non[i]

bool_spe = np.zeros(num_samp.size)
for i in range(num_spe.size):
    bool_spe+= num_samp==num_spe[i]

bool_plot = bool_plot_1 * bool_plot_2 * bool_plot_3 - bool_non + bool_spe
print(np.sum(bool_plot))

mag0_iso_plot = mag0_iso_samp[bool_plot!=0]
mag1_iso_plot = mag1_iso_samp[bool_plot!=0]
mag0_auto_plot = mag0_auto_samp[bool_plot!=0]
mag1_auto_plot = mag1_auto_samp[bool_plot!=0]
mag0_auto_u_plot = mag0_auto_u_samp[bool_plot!=0]
mag0_auto_d_plot = mag0_auto_d_samp[bool_plot!=0]
mag0_auto_err_plot = np.sqrt(mag0_auto_u_plot**2 + mag0_auto_d_plot**2)
mag1_auto_u_plot = mag1_auto_u_samp[bool_plot!=0]
mag1_auto_d_plot = mag1_auto_d_samp[bool_plot!=0]
mag1_auto_err_plot = np.sqrt(mag1_auto_u_plot**2 + mag1_auto_d_plot**2)
num_plot = num_samp[bool_plot!=0]


plt.figure()
# plt.subplot(1,2,1)
plt.scatter(mag0_iso_samp, color_samp, c = 'k', s = 5)
plt.errorbar(mag0_iso_plot, mag0_auto_plot - mag1_auto_plot, yerr = np.sqrt(mag0_auto_err_plot**2 + mag1_auto_err_plot**2), c = 'orange', fmt = 'x')
# plt.scatter(mag0_iso_plot, mag0_auto_plot - mag1_auto_plot, c = 'orange', s = 5)
plt.scatter(mag0_iso_spe, mag0_auto_spe - mag1_auto_spe, c = 'r', s = 5)
plt.plot(x_cmd, linear(x_cmd, *popt), c = 'r', linestyle = '--')
plt.plot(x_cmd, linear(x_cmd, *popt) + test_std, c = 'b', linestyle = '--')
plt.plot(x_cmd, linear(x_cmd, *popt) - test_std, c = 'b', linestyle = '--')
plt.xlim(18,26)
plt.ylim(-2,2)
plt.xlabel('F105W (ISO)')
plt.ylabel('F105W - F140W (Auto)')
# plt.subplot(1,2,2)
# plt.hist(test, bins = 100)
plt.show()

for i in range(np.size(num_plot)):
    seg[seg==num_plot[i]] = 0

for i in range(np.size(num_plot)):
    seg_max[seg_max==num_plot[i]] = 0


hdu = fits.PrimaryHDU(seg)
hdu.writeto('%s_seg_mem_banned.fits'%target, overwrite = True)

hdu = fits.PrimaryHDU(seg_max)
hdu.writeto('%s_seg_6.0_mem_banned.fits'%target, overwrite = True)