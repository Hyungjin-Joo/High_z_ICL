import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

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
sci1_data = sci1[0].data
sci0_head = sci0[0].header
sci1_head = sci1[0].header
sma, sfb0, sfb1, sfe0, sfe1 = np.loadtxt('%s_SFB.txt'%target, unpack = True)
photflam = [sci0_head['PHOTFLAM'], sci1_head['PHOTFLAM']]
photplam = [sci0_head['PHOTPLAM'], sci1_head['PHOTPLAM']]
exptime = [sci0_head['EXPTIME'], sci1_head['EXPTIME']]

### counts/arcs > magnitude
def mag(flux, filter):
    abmag_zpt = -2.5 * np.log10(photflam[filter]) - 21.10 - 5 * np.log10(photplam[filter]) + 18.692
    mag = -2.5*np.log10(flux/exptime[filter]) + abmag_zpt
    return mag

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


### plotting
plt.figure()
#plt.errorbar(sma, mag0, yerr = [mag0_uerr, mag0_derr], fmt = 'o', c = 'b', label = 'F105W')
plt.errorbar(sma, mag1, yerr = [mag1_uerr, mag1_derr], fmt = 'o', c = 'k', label = 'F140W')
#plt.axhline(mag(ske0/0.05**2,0), c = 'b', linestyle = '--')
plt.axhline(mag(ske1/0.05**2,1), c = 'k', linestyle = '--')
plt.xscale('log')
plt.xlabel('SMA [kpc]')
plt.ylabel(r'$SB [mag \ arcs^{-2} ]$')
#plt.legend(loc=0)
plt.ylim(31,18)
plt.title(target)
plt.show()

color_uerr = np.sqrt(mag0_uerr**2 + mag1_uerr**2)
color_derr = np.sqrt(mag0_derr**2 + mag1_derr**2)
color = mag0 - mag1

plt.figure()
plt.errorbar(sma, color, yerr = [color_derr, color_uerr], fmt = 'o', c = 'k')
plt.xscale('log')
plt.xlabel('SMA [kpc]')
plt.ylabel('F105W - F140W')
plt.title(target)
#plt.show()

plt.figure()
plt.plot(sma, mag0, c = 'b', label = 'F105W')
plt.fill_between(sma, mag0+mag0_derr, mag0-mag0_uerr, color = 'b', alpha = 0.5)
plt.plot(sma, mag1, c = 'r', label = 'F140W')
plt.fill_between(sma, mag1+mag1_derr, mag1-mag1_uerr, color = 'r', alpha = 0.5)
plt.axhline(mag(ske0/0.05**2,0), c = 'b', linestyle = '--')
plt.axhline(mag(ske1/0.05**2,1), c = 'r', linestyle = '--')
plt.xscale('log')
plt.xlabel('SMA [kpc]')
plt.ylabel(r'$SB [mag \ arcs^{-2} ]$')
plt.legend(loc=0)
plt.ylim(31,20)
plt.title(target)
plt.savefig('%s_sfb.png'%target)
plt.close()

color_uerr = np.sqrt(mag0_uerr**2 + mag1_uerr**2)
color_derr = np.sqrt(mag0_derr**2 + mag1_derr**2)
color = mag0 - mag1

plt.figure(dpi = 500)
plt.plot(sma, color, c = 'k')
plt.fill_between(sma, color-color_derr, color+color_uerr, color = 'k', alpha = 0.5)
plt.xscale('log')
plt.xlabel('SMA [kpc]')
plt.ylabel('F105W - F140W')
plt.title(target)
plt.xlim(0.35,213)
plt.savefig('%s_color.png'%target)
plt.close()

plt.figure()
plt.subplot(1,2,1)
plt.scatter(sma[1:], np.diff(np.log10(sfb0)) / (np.diff(np.log10(sma))), c= 'k', s=10)
plt.xscale('log')
plt.subplot(1,2,2)
plt.scatter(sma[1:], np.diff(np.log10(sfb1)) / (np.diff(np.log10(sma))), c= 'k', s=10)
plt.xscale('log')
#plt.show()

data = np.array([sma.T, mag0.T, mag0_uerr.T, mag0_derr.T, mag1.T, mag1_uerr.T, mag1_derr.T])
header = 'sma[kpc], sfb[mag/arcs2], sfb_ue, sfb_de'
np.savetxt('%s_SFB_mag.txt'%target, data.T, header = header)


plt.figure()
plt.plot(sma, sfe0)
plt.plot(sma, sfe1)
plt.xscale('log')
plt.yscale('log')
#plt.show()

data = np.array([sma.T, color.T, color_uerr.T, color_derr.T])
header = 'sma[kpc], sfb[mag/arcs2], sfb_ue, sfb_de'
np.savetxt('%s_color.txt'%target, data.T, header = header)

bool_nan =  1  - np.isnan(color_uerr * color_derr)
zeros = np.zeros(np.size(sma))
if filters[1] == 'F160W':
    zeros += 1
objid = np.linspace(0, np.size(sma)-1, np.size(sma))

sfb0_2 = sfb0 / exptime[0]
sfb1_2 = sfb1 / exptime[1]
sfe0_2 = sfe0 / exptime[0]
sfe1_2 = sfe1 / exptime[1]

data = np.array([objid.T, zeros.T, mag0, mag1, sfb0_2.T, sfb1_2.T, sfe0_2.T, sfe1_2.T])
header = 'objid filterset mag0 mag1 cnt0, cnt1, err0, err1'
np.savetxt('%s_photometry_prospector.dat'%target, data.T, header = header)