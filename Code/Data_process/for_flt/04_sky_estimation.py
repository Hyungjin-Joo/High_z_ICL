import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
import copy as cp
from scipy import stats
from scipy.optimize import curve_fit

files = input('File name?')
fname, type = files.split('.')

mask10 = fits.open('%s_check_1.0.fits'%fname)
mask11 = fits.open('%s_check_1.1.fits'%fname)
mask12 = fits.open('%s_check_1.2.fits'%fname)
mask13 = fits.open('%s_check_1.3.fits'%fname)
mask14 = fits.open('%s_check_1.4.fits'%fname)
mask15 = fits.open('%s_check_1.5.fits'%fname)
mask16 = fits.open('%s_check_1.6.fits'%fname)
mask17 = fits.open('%s_check_1.7.fits'%fname)
mask18 = fits.open('%s_check_1.8.fits'%fname)
mask19 = fits.open('%s_check_1.9.fits'%fname)
mask20 = fits.open('%s_check_2.0.fits'%fname)
mask21 = fits.open('%s_check_2.1.fits'%fname)
mask22 = fits.open('%s_check_2.2.fits'%fname)
mask23 = fits.open('%s_check_2.3.fits'%fname)
mask24 = fits.open('%s_check_2.4.fits'%fname)
mask25 = fits.open('%s_check_2.5.fits'%fname)
mask26 = fits.open('%s_check_2.6.fits'%fname)
mask27 = fits.open('%s_check_2.7.fits'%fname)
mask28 = fits.open('%s_check_2.8.fits'%fname)
mask29 = fits.open('%s_check_2.9.fits'%fname)
mask30 = fits.open('%s_check_3.0.fits'%fname)
mask31 = fits.open('%s_check_3.1.fits'%fname)
mask32 = fits.open('%s_check_3.2.fits'%fname)
mask33 = fits.open('%s_check_3.3.fits'%fname)
mask34 = fits.open('%s_check_3.4.fits'%fname)
mask35 = fits.open('%s_check_3.5.fits'%fname)
mask36 = fits.open('%s_check_3.6.fits'%fname)
mask37 = fits.open('%s_check_3.7.fits'%fname)
mask38 = fits.open('%s_check_3.8.fits'%fname)
mask39 = fits.open('%s_check_3.9.fits'%fname)
mask40 = fits.open('%s_check_4.0.fits'%fname)
mask41 = fits.open('%s_check_4.1.fits'%fname)
mask42 = fits.open('%s_check_4.2.fits'%fname)
mask43 = fits.open('%s_check_4.3.fits'%fname)
mask44 = fits.open('%s_check_4.4.fits'%fname)
mask45 = fits.open('%s_check_4.5.fits'%fname)
mask46 = fits.open('%s_check_4.6.fits'%fname)
mask47 = fits.open('%s_check_4.7.fits'%fname)
mask48 = fits.open('%s_check_4.8.fits'%fname)
mask49 = fits.open('%s_check_4.9.fits'%fname)
mask50 = fits.open('%s_check_5.0.fits'%fname)
mask51 = fits.open('%s_check_5.1.fits'%fname)
mask52 = fits.open('%s_check_5.2.fits'%fname)
mask53 = fits.open('%s_check_5.3.fits'%fname)
mask54 = fits.open('%s_check_5.4.fits'%fname)
mask55 = fits.open('%s_check_5.5.fits'%fname)
mask56 = fits.open('%s_check_5.6.fits'%fname)
mask57 = fits.open('%s_check_5.7.fits'%fname)
mask58 = fits.open('%s_check_5.8.fits'%fname)
mask59 = fits.open('%s_check_5.9.fits'%fname)
mask60 = fits.open('%s_check_6.0.fits'%fname)
mask61 = fits.open('%s_check_6.1.fits'%fname)
mask62 = fits.open('%s_check_6.2.fits'%fname)
mask63 = fits.open('%s_check_6.3.fits'%fname)
mask64 = fits.open('%s_check_6.4.fits'%fname)
mask65 = fits.open('%s_check_6.5.fits'%fname)
mask66 = fits.open('%s_check_6.6.fits'%fname)
mask67 = fits.open('%s_check_6.7.fits'%fname)
mask68 = fits.open('%s_check_6.8.fits'%fname)
mask69 = fits.open('%s_check_6.9.fits'%fname)
mask70 = fits.open('%s_check_7.0.fits'%fname)
mask71 = fits.open('%s_check_7.1.fits'%fname)
mask72 = fits.open('%s_check_7.2.fits'%fname)
mask73 = fits.open('%s_check_7.3.fits'%fname)
mask74 = fits.open('%s_check_7.4.fits'%fname)
mask75 = fits.open('%s_check_7.5.fits'%fname)
mask76 = fits.open('%s_check_7.6.fits'%fname)
mask77 = fits.open('%s_check_7.7.fits'%fname)
mask78 = fits.open('%s_check_7.8.fits'%fname)
mask79 = fits.open('%s_check_7.9.fits'%fname)
mask80 = fits.open('%s_check_8.0.fits'%fname)
mask81 = fits.open('%s_check_8.1.fits'%fname)
mask82 = fits.open('%s_check_8.2.fits'%fname)
mask83 = fits.open('%s_check_8.3.fits'%fname)
mask84 = fits.open('%s_check_8.4.fits'%fname)
mask85 = fits.open('%s_check_8.5.fits'%fname)
mask86 = fits.open('%s_check_8.6.fits'%fname)
mask87 = fits.open('%s_check_8.7.fits'%fname)
mask88 = fits.open('%s_check_8.8.fits'%fname)
mask89 = fits.open('%s_check_8.9.fits'%fname)
mask90 = fits.open('%s_check_9.0.fits'%fname)
mask91 = fits.open('%s_check_9.1.fits'%fname)
mask92 = fits.open('%s_check_9.2.fits'%fname)
mask93 = fits.open('%s_check_9.3.fits'%fname)
mask94 = fits.open('%s_check_9.4.fits'%fname)
mask95 = fits.open('%s_check_9.5.fits'%fname)
mask96 = fits.open('%s_check_9.6.fits'%fname)
mask97 = fits.open('%s_check_9.7.fits'%fname)
mask98 = fits.open('%s_check_9.8.fits'%fname)
mask99 = fits.open('%s_check_9.9.fits'%fname)
mask100 = fits.open('%s_check_10.0.fits'%fname)
mask101 = fits.open('%s_check_10.1.fits'%fname)
mask102 = fits.open('%s_check_10.2.fits'%fname)
mask103 = fits.open('%s_check_10.3.fits'%fname)
mask104 = fits.open('%s_check_10.4.fits'%fname)
mask105 = fits.open('%s_check_10.5.fits'%fname)
mask106 = fits.open('%s_check_10.6.fits'%fname)
mask107 = fits.open('%s_check_10.7.fits'%fname)
mask108 = fits.open('%s_check_10.8.fits'%fname)
mask109 = fits.open('%s_check_10.9.fits'%fname)
mask110 = fits.open('%s_check_11.0.fits'%fname)
mask111 = fits.open('%s_check_11.1.fits'%fname)
mask112 = fits.open('%s_check_11.2.fits'%fname)
mask113 = fits.open('%s_check_11.3.fits'%fname)
mask114 = fits.open('%s_check_11.4.fits'%fname)
mask115 = fits.open('%s_check_11.5.fits'%fname)
mask116 = fits.open('%s_check_11.6.fits'%fname)
mask117 = fits.open('%s_check_11.7.fits'%fname)
mask118 = fits.open('%s_check_11.8.fits'%fname)
mask119 = fits.open('%s_check_11.9.fits'%fname)
mask120 = fits.open('%s_check_12.0.fits'%fname)
m10 = mask10[1].data
m11 = mask11[1].data
m12 = mask12[1].data
m13 = mask13[1].data
m14 = mask14[1].data
m15 = mask15[1].data
m16 = mask16[1].data
m17 = mask17[1].data
m18 = mask18[1].data
m19 = mask19[1].data
m20 = mask20[1].data
m21 = mask21[1].data
m22 = mask22[1].data
m23 = mask23[1].data
m24 = mask24[1].data
m25 = mask25[1].data
m26 = mask26[1].data
m27 = mask27[1].data
m28 = mask28[1].data
m29 = mask29[1].data
m30 = mask30[1].data
m31 = mask31[1].data
m32 = mask32[1].data
m33 = mask33[1].data
m34 = mask34[1].data
m35 = mask35[1].data
m36 = mask36[1].data
m37 = mask37[1].data
m38 = mask38[1].data
m39 = mask39[1].data
m40 = mask40[1].data
m41 = mask41[1].data
m42 = mask42[1].data
m43 = mask43[1].data
m44 = mask44[1].data
m45 = mask45[1].data
m46 = mask46[1].data
m47 = mask47[1].data
m48 = mask48[1].data
m49 = mask49[1].data
m50 = mask50[1].data
m51 = mask51[1].data
m52 = mask52[1].data
m53 = mask53[1].data
m54 = mask54[1].data
m55 = mask55[1].data
m56 = mask56[1].data
m57 = mask57[1].data
m58 = mask58[1].data
m59 = mask59[1].data
m60 = mask60[1].data
m61 = mask61[1].data
m62 = mask62[1].data
m63 = mask63[1].data
m64 = mask64[1].data
m65 = mask65[1].data
m66 = mask66[1].data
m67 = mask67[1].data
m68 = mask68[1].data
m69 = mask69[1].data
m70 = mask70[1].data
m71 = mask71[1].data
m72 = mask72[1].data
m73 = mask73[1].data
m74 = mask74[1].data
m75 = mask75[1].data
m76 = mask76[1].data
m77 = mask77[1].data
m78 = mask78[1].data
m79 = mask79[1].data
m80 = mask80[1].data
m81 = mask81[1].data
m82 = mask82[1].data
m83 = mask83[1].data
m84 = mask84[1].data
m85 = mask85[1].data
m86 = mask86[1].data
m87 = mask87[1].data
m88 = mask88[1].data
m89 = mask89[1].data
m90 = mask90[1].data
m91 = mask91[1].data
m92 = mask92[1].data
m93 = mask93[1].data
m94 = mask94[1].data
m95 = mask95[1].data
m96 = mask96[1].data
m97 = mask97[1].data
m98 = mask98[1].data
m99 = mask99[1].data
m100 = mask100[1].data
m101 = mask101[1].data
m102 = mask102[1].data
m103 = mask103[1].data
m104 = mask104[1].data
m105 = mask105[1].data
m106 = mask106[1].data
m107 = mask107[1].data
m108 = mask108[1].data
m109 = mask109[1].data
m110 = mask110[1].data
m111 = mask111[1].data
m112 = mask112[1].data
m113 = mask113[1].data
m114 = mask114[1].data
m115 = mask115[1].data
m116 = mask116[1].data
m117 = mask117[1].data
m118 = mask118[1].data
m119 = mask119[1].data
m120 = mask120[1].data

masks = np.array([m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56,m57,m58,m59,m60,m61,m62,m63,m64,m65,m66,m67,m68,m69,m70,m71,m72,m73,m74,m75,m76,m77,m78,m79,m80,m81,m82,m83,m84,m85,m86,m87,m88,m89,m90,m91,m92,m93,m94,m95,m96,m97,m98,m99,m100,m101,m102,m103,m104,m105,m106,m107,m108,m109,m110,m111,m112,m113,m114,m115,m116,m117,m118,m119,m120])

ver = np.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12.0])

sci105 = fits.open(files)
s105 = sci105[1].data
wht105 = fits.open('%s_wht.%s'%(fname,type))
w105 = wht105[0].data
s105[w105==0] = -32768

ny, nx = np.shape(s105)
x = np.linspace(1,nx,nx)
y = np.linspace(1,ny,ny)
xv,yv = np.meshgrid(x,y)

mean105st = np.zeros(ver.size)
ste105st = np.zeros(ver.size)
for i in range(ver.size):
    s105t = cp.deepcopy(s105)
    bool_mask = (masks[i,:,:]!=0)
    subsky = np.zeros(32)
    for j in range(32):
        bool_subsky = (w105==j+1)
        bool = bool_mask * bool_subsky
        subsky[j] = np.nanmedian(s105t[bool!=0])
    mean105st[i] = np.nanmedian(subsky)
    ste105st[i] = np.nanmedian(subsky)
    print(i+1,'/', ver.size,'ing...             ', end ='\r')

means105 = np.zeros(ver.size-1)
ste105 = np.zeros(ver.size-1)
for i in range(ver.size - 1):
    dmask = (masks[i+1,:,:] - masks[i,:,:])!=0
    sample105 = cp.deepcopy(s105)
    subsky = np.zeros(32)
    for j in range(32):
        bool_subsky = (w105==j+1)
        bool = dmask * bool_subsky
        subsky[j] = np.nanmedian(sample105[bool!=0])
    means105[i] = np.nanmedian(subsky)
    ste105[i] = np.nanmedian(subsky)
    print(i, '/', ver.size,'ing...             ', end ='\r')

plt.figure()
plt.plot(ver, mean105st, c = 'green', label = 'Mask')
plt.plot(ver[1:], means105, c = 'purple', label = 'dMask')
plt.legend(loc=0)
plt.xlabel(r'$p_e$')
plt.ylabel(r'$count\/arcs^2$')
plt.title(r'%s'%fname)
plt.savefig(r'%s.png'%fname)
plt.close()

sky = np.nanmedian(means105)
fnames = np.loadtxt('sky_level.txt', usecols = (0), dtype = 'str', unpack = True)
skies = np.loadtxt('sky_level.txt', usecols = (1), dtype = 'float', unpack = True)
index = np.where(fnames == files)
skies[index] = format(sky,'.5f')

sky_file = np.concatenate(([fnames], [skies], [skies]))
sky_file = sky_file.reshape(3,fnames.size)
np.savetxt('sky_level.txt',sky_file.T, fmt = '%s')
