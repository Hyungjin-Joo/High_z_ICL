import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.stats import sigma_clip
import copy as cp
from scipy import stats
from scipy.optimize import curve_fit

Target = input('Target ')
step = 1
mask10 = fits.open('%s_seg_1.0.fits'%Target)
mask11 = fits.open('%s_seg_1.1.fits'%Target)
mask12 = fits.open('%s_seg_1.2.fits'%Target)
mask13 = fits.open('%s_seg_1.3.fits'%Target)
mask14 = fits.open('%s_seg_1.4.fits'%Target)
mask15 = fits.open('%s_seg_1.5.fits'%Target)
mask16 = fits.open('%s_seg_1.6.fits'%Target)
mask17 = fits.open('%s_seg_1.7.fits'%Target)
mask18 = fits.open('%s_seg_1.8.fits'%Target)
mask19 = fits.open('%s_seg_1.9.fits'%Target)
mask20 = fits.open('%s_seg_2.0.fits'%Target)
mask21 = fits.open('%s_seg_2.1.fits'%Target)
mask22 = fits.open('%s_seg_2.2.fits'%Target)
mask23 = fits.open('%s_seg_2.3.fits'%Target)
mask24 = fits.open('%s_seg_2.4.fits'%Target)
mask25 = fits.open('%s_seg_2.5.fits'%Target)
mask26 = fits.open('%s_seg_2.6.fits'%Target)
mask27 = fits.open('%s_seg_2.7.fits'%Target)
mask28 = fits.open('%s_seg_2.8.fits'%Target)
mask29 = fits.open('%s_seg_2.9.fits'%Target)
mask30 = fits.open('%s_seg_3.0.fits'%Target)
mask31 = fits.open('%s_seg_3.1.fits'%Target)
mask32 = fits.open('%s_seg_3.2.fits'%Target)
mask33 = fits.open('%s_seg_3.3.fits'%Target)
mask34 = fits.open('%s_seg_3.4.fits'%Target)
mask35 = fits.open('%s_seg_3.5.fits'%Target)
mask36 = fits.open('%s_seg_3.6.fits'%Target)
mask37 = fits.open('%s_seg_3.7.fits'%Target)
mask38 = fits.open('%s_seg_3.8.fits'%Target)
mask39 = fits.open('%s_seg_3.9.fits'%Target)
mask40 = fits.open('%s_seg_4.0.fits'%Target)
mask41 = fits.open('%s_seg_4.1.fits'%Target)
mask42 = fits.open('%s_seg_4.2.fits'%Target)
mask43 = fits.open('%s_seg_4.3.fits'%Target)
mask44 = fits.open('%s_seg_4.4.fits'%Target)
mask45 = fits.open('%s_seg_4.5.fits'%Target)
mask46 = fits.open('%s_seg_4.6.fits'%Target)
mask47 = fits.open('%s_seg_4.7.fits'%Target)
mask48 = fits.open('%s_seg_4.8.fits'%Target)
mask49 = fits.open('%s_seg_4.9.fits'%Target)
mask50 = fits.open('%s_seg_5.0.fits'%Target)
mask51 = fits.open('%s_seg_5.1.fits'%Target)
mask52 = fits.open('%s_seg_5.2.fits'%Target)
mask53 = fits.open('%s_seg_5.3.fits'%Target)
mask54 = fits.open('%s_seg_5.4.fits'%Target)
mask55 = fits.open('%s_seg_5.5.fits'%Target)
mask56 = fits.open('%s_seg_5.6.fits'%Target)
mask57 = fits.open('%s_seg_5.7.fits'%Target)
mask58 = fits.open('%s_seg_5.8.fits'%Target)
mask59 = fits.open('%s_seg_5.9.fits'%Target)
mask60 = fits.open('%s_seg_6.0.fits'%Target)
mask61 = fits.open('%s_seg_6.1.fits'%Target)
mask62 = fits.open('%s_seg_6.2.fits'%Target)
mask63 = fits.open('%s_seg_6.3.fits'%Target)
mask64 = fits.open('%s_seg_6.4.fits'%Target)
mask65 = fits.open('%s_seg_6.5.fits'%Target)
mask66 = fits.open('%s_seg_6.6.fits'%Target)
mask67 = fits.open('%s_seg_6.7.fits'%Target)
mask68 = fits.open('%s_seg_6.8.fits'%Target)
mask69 = fits.open('%s_seg_6.9.fits'%Target)
mask70 = fits.open('%s_seg_7.0.fits'%Target)
mask71 = fits.open('%s_seg_7.1.fits'%Target)
mask72 = fits.open('%s_seg_7.2.fits'%Target)
mask73 = fits.open('%s_seg_7.3.fits'%Target)
mask74 = fits.open('%s_seg_7.4.fits'%Target)
mask75 = fits.open('%s_seg_7.5.fits'%Target)
mask76 = fits.open('%s_seg_7.6.fits'%Target)
mask77 = fits.open('%s_seg_7.7.fits'%Target)
mask78 = fits.open('%s_seg_7.8.fits'%Target)
mask79 = fits.open('%s_seg_7.9.fits'%Target)
mask80 = fits.open('%s_seg_8.0.fits'%Target)
mask81 = fits.open('%s_seg_8.1.fits'%Target)
mask82 = fits.open('%s_seg_8.2.fits'%Target)
mask83 = fits.open('%s_seg_8.3.fits'%Target)
mask84 = fits.open('%s_seg_8.4.fits'%Target)
mask85 = fits.open('%s_seg_8.5.fits'%Target)
mask86 = fits.open('%s_seg_8.6.fits'%Target)
mask87 = fits.open('%s_seg_8.7.fits'%Target)
mask88 = fits.open('%s_seg_8.8.fits'%Target)
mask89 = fits.open('%s_seg_8.9.fits'%Target)
mask90 = fits.open('%s_seg_9.0.fits'%Target)
mask91 = fits.open('%s_seg_9.1.fits'%Target)
mask92 = fits.open('%s_seg_9.2.fits'%Target)
mask93 = fits.open('%s_seg_9.3.fits'%Target)
mask94 = fits.open('%s_seg_9.4.fits'%Target)
mask95 = fits.open('%s_seg_9.5.fits'%Target)
mask96 = fits.open('%s_seg_9.6.fits'%Target)
mask97 = fits.open('%s_seg_9.7.fits'%Target)
mask98 = fits.open('%s_seg_9.8.fits'%Target)
mask99 = fits.open('%s_seg_9.9.fits'%Target)
mask100 = fits.open('%s_seg_10.0.fits'%Target)
mask101 = fits.open('%s_seg_10.1.fits'%Target)
mask102 = fits.open('%s_seg_10.2.fits'%Target)
mask103 = fits.open('%s_seg_10.3.fits'%Target)
mask104 = fits.open('%s_seg_10.4.fits'%Target)
mask105 = fits.open('%s_seg_10.5.fits'%Target)
mask106 = fits.open('%s_seg_10.6.fits'%Target)
mask107 = fits.open('%s_seg_10.7.fits'%Target)
mask108 = fits.open('%s_seg_10.8.fits'%Target)
mask109 = fits.open('%s_seg_10.9.fits'%Target)
mask110 = fits.open('%s_seg_11.0.fits'%Target)
mask110 = fits.open('%s_seg_11.0.fits'%Target)
mask111 = fits.open('%s_seg_11.1.fits'%Target)
mask112 = fits.open('%s_seg_11.2.fits'%Target)
mask113 = fits.open('%s_seg_11.3.fits'%Target)
mask114 = fits.open('%s_seg_11.4.fits'%Target)
mask115 = fits.open('%s_seg_11.5.fits'%Target)
mask116 = fits.open('%s_seg_11.6.fits'%Target)
mask117 = fits.open('%s_seg_11.7.fits'%Target)
mask118 = fits.open('%s_seg_11.8.fits'%Target)
mask119 = fits.open('%s_seg_11.9.fits'%Target)
mask120 = fits.open('%s_seg_12.0.fits'%Target)
mask121 = fits.open('%s_seg_12.1.fits'%Target)
mask122 = fits.open('%s_seg_12.2.fits'%Target)
mask123 = fits.open('%s_seg_12.3.fits'%Target)
mask124 = fits.open('%s_seg_12.4.fits'%Target)
mask125 = fits.open('%s_seg_12.5.fits'%Target)
mask126 = fits.open('%s_seg_12.6.fits'%Target)
mask127 = fits.open('%s_seg_12.7.fits'%Target)
mask128 = fits.open('%s_seg_12.8.fits'%Target)
mask129 = fits.open('%s_seg_12.9.fits'%Target)
mask130 = fits.open('%s_seg_13.0.fits'%Target)

m10 = mask10[0].data
m11 = mask11[0].data
m12 = mask12[0].data
m13 = mask13[0].data
m14 = mask14[0].data
m15 = mask15[0].data
m16 = mask16[0].data
m17 = mask17[0].data
m18 = mask18[0].data
m19 = mask19[0].data
m20 = mask20[0].data
m21 = mask21[0].data
m22 = mask22[0].data
m23 = mask23[0].data
m24 = mask24[0].data
m25 = mask25[0].data
m26 = mask26[0].data
m27 = mask27[0].data
m28 = mask28[0].data
m29 = mask29[0].data
m30 = mask30[0].data
m31 = mask31[0].data
m32 = mask32[0].data
m33 = mask33[0].data
m34 = mask34[0].data
m35 = mask35[0].data
m36 = mask36[0].data
m37 = mask37[0].data
m38 = mask38[0].data
m39 = mask39[0].data
m40 = mask40[0].data
m41 = mask41[0].data
m42 = mask42[0].data
m43 = mask43[0].data
m44 = mask44[0].data
m45 = mask45[0].data
m46 = mask46[0].data
m47 = mask47[0].data
m48 = mask48[0].data
m49 = mask49[0].data
m50 = mask50[0].data
m51 = mask51[0].data
m52 = mask52[0].data
m53 = mask53[0].data
m54 = mask54[0].data
m55 = mask55[0].data
m56 = mask56[0].data
m57 = mask57[0].data
m58 = mask58[0].data
m59 = mask59[0].data
m60 = mask60[0].data
m61 = mask61[0].data
m62 = mask62[0].data
m63 = mask63[0].data
m64 = mask64[0].data
m65 = mask65[0].data
m66 = mask66[0].data
m67 = mask67[0].data
m68 = mask68[0].data
m69 = mask69[0].data
m70 = mask70[0].data
m71 = mask71[0].data
m72 = mask72[0].data
m73 = mask73[0].data
m74 = mask74[0].data
m75 = mask75[0].data
m76 = mask76[0].data
m77 = mask77[0].data
m78 = mask78[0].data
m79 = mask79[0].data
m80 = mask80[0].data
m81 = mask81[0].data
m82 = mask82[0].data
m83 = mask83[0].data
m84 = mask84[0].data
m85 = mask85[0].data
m86 = mask86[0].data
m87 = mask87[0].data
m88 = mask88[0].data
m89 = mask89[0].data
m90 = mask90[0].data
m91 = mask91[0].data
m92 = mask92[0].data
m93 = mask93[0].data
m94 = mask94[0].data
m95 = mask95[0].data
m96 = mask96[0].data
m97 = mask97[0].data
m98 = mask98[0].data
m99 = mask99[0].data
m100 = mask100[0].data
m101 = mask101[0].data
m102 = mask102[0].data
m103 = mask103[0].data
m104 = mask104[0].data
m105 = mask105[0].data
m106 = mask106[0].data
m107 = mask107[0].data
m108 = mask108[0].data
m109 = mask109[0].data
m110 = mask110[0].data
m111 = mask111[0].data
m112 = mask112[0].data
m113 = mask113[0].data
m114 = mask114[0].data
m115 = mask115[0].data
m116 = mask116[0].data
m117 = mask117[0].data
m118 = mask118[0].data
m119 = mask119[0].data
m120 = mask120[0].data
m121 = mask121[0].data
m122 = mask122[0].data
m123 = mask123[0].data
m124 = mask124[0].data
m125 = mask125[0].data
m126 = mask126[0].data
m127 = mask127[0].data
m128 = mask128[0].data
m129 = mask129[0].data
m130 = mask130[0].data

masks = np.array([m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56,m57,m58,m59,m60,m61,m62,m63,m64,m65,m66,m67,m68,m69,m70,m71,m72,m73,m74,m75,m76,m77,m78,m79,m80,m81,m82,m83,m84,m85,m86,m87,m88,m89,m90,m91,m92,m93,m94,m95,m96,m97,m98,m99,m100,m101,m102,m103,m104,m105,m106,m107,m108,m109,m110,m111,m112,m113,m114,m115,m116,m117,m118,m119,m120,m121,m122,m123,m124,m125,m126,m127,m128,m129,m130])

ver = np.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12.0,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13.0])


sci105 = fits.open('%s_F105W_drz_sci.fits'%Target)
sci140 = fits.open('%s_F160W_drz_sci.fits'%Target)
s105 = sci105[0].data
s140 = sci140[0].data

s105[np.isnan(s105)]=-32768
s140[np.isnan(s140)]=-32768

ny, nx = np.shape(s140)
x = np.linspace(1,nx,nx)
y = np.linspace(1,ny,ny)
xv,yv = np.meshgrid(x,y)
rs = np.sqrt((xv-np.int64(ny/2))**2 + (yv-np.int64(nx/2))**2)

r0 = 1100
r1 = 1200
r2 = 1300
r3 = 1400
r4 = 1500

R00 = (rs-r0)>0
R10 = (rs-r1)<=0
R20 = (rs-r1)>0
R30 = (rs-r2)<=0
R40 = (rs-r2)>0
R50 = (rs-r3)<=0
R60 = (rs-r3)>0
R70 = (rs-r4)<=0

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

mean105 = np.zeros((ver.size,32))
mean140 = np.zeros((ver.size,32))


for i in range(ver.size):
    s105t = cp.deepcopy(s105)
    s105t[masks[i,:,:]!=0]=-32768
    s140t = cp.deepcopy(s140)
    s140t[masks[i,:,:]!=0]=-32768
    for j in range(32):
        sky_candidate_105 = s105t[sub_sky==j+1]
        sky_candidate_105 = sky_candidate_105[sky_candidate_105!=-32768]
        mean105[i,j] = np.median(sky_candidate_105)
        sky_candidate_140 = s140t[sub_sky==j+1]
        sky_candidate_140 = sky_candidate_140[sky_candidate_140!=-32768]
        mean140[i,j] = np.median(sky_candidate_140)
    print(i+1,'/', ver.size,'ing...          ', end = '\r')

colors = np.array(['red','orange','gold','green','skyblue','blue','navy','purple'])

size105 = np.sqrt(32 - np.sum(np.isnan(mean105), axis = 1))
size140 = np.sqrt(32 - np.sum(np.isnan(mean140), axis = 1))
plt.figure()
plt.subplot(1,2,1)
for i in range(32):
    plt.plot(ver, mean105[:,i], color = colors[i%8], alpha = 0.5)
real_sammeans105 = np.nanmedian(mean105, axis = 1)
plt.plot(ver, real_sammeans105, c = 'k')
plt.plot(ver, real_sammeans105 + np.nanstd(mean105, axis =1)/size105, linestyle = '--', c = 'k')
plt.plot(ver, real_sammeans105 - np.nanstd(mean105, axis =1)/size105, linestyle = '--', c = 'k')
plt.title('F105W')
plt.xlabel('p_e')
plt.ylabel('counts/arcs2')
plt.subplot(1,2,2)
for i in range(32):
    plt.plot(ver, mean140[:,i], color = colors[i%8], alpha = 0.5)
real_sammeans140 = np.nanmedian(mean140, axis = 1)
plt.plot(ver, real_sammeans140, c = 'k')
plt.plot(ver, real_sammeans140 + np.nanstd(mean140, axis =1)/size140, linestyle = '--', c = 'k')
plt.plot(ver, real_sammeans140 - np.nanstd(mean140, axis =1)/size140, linestyle = '--', c = 'k')
plt.title('F160W')
plt.xlabel('p_e')
plt.ylabel('counts/arcs2')
plt.savefig('%s_sky_vs_pe.png'%Target)
plt.close()

mean105_d = np.zeros((ver.size-step,32))
mean140_d = np.zeros((ver.size-step,32))

for i in range(ver.size - step):
    dmask = masks[i+step,:,:] - masks[i,:,:]
    sample105 = cp.deepcopy(s105)
    sample105[dmask==0]=-32768
    sample140 = cp.deepcopy(s140)
    sample140[dmask==0]=-32768
    for j in range(32):
        sky_candidate_105 = sample105[sub_sky==j+1]
        sky_candidate_105 = sky_candidate_105[sky_candidate_105!=-32768]
        mean105_d[i,j] = np.median(sky_candidate_105)
        sky_candidate_140 = sample140[sub_sky==j+1]
        sky_candidate_140 = sky_candidate_140[sky_candidate_140!=-32768]
        mean140_d[i,j] = np.median(sky_candidate_140)
    print(i+1,'/', ver.size,'ing...          ', end = '\r')

size105_d = np.sqrt(32 - np.sum(np.isnan(mean105_d), axis = 1))
size140_d = np.sqrt(32 - np.sum(np.isnan(mean140_d), axis = 1))

plt.figure()
plt.subplot(1,2,1)
for i in range(32):
    plt.plot(ver[step:], mean105_d[:,i], color = colors[i%8], alpha = 0.5)
real_sammeans105_d = np.nanmedian(mean105_d, axis = 1)
plt.plot(ver[step:], real_sammeans105_d, c = 'k')
plt.plot(ver[step:], real_sammeans105_d + np.nanstd(mean105_d, axis =1)/size105_d, linestyle = '--', c = 'k')
plt.plot(ver[step:], real_sammeans105_d - np.nanstd(mean105_d, axis =1)/size105_d, linestyle = '--', c = 'k')
plt.title('F105W')
plt.xlabel('p_e')
plt.ylabel('counts/arcs2')
plt.subplot(1,2,2)
for i in range(32):
    plt.plot(ver[step:], mean140_d[:,i], color = colors[i%8], alpha = 0.5)
real_sammeans140_d = np.nanmedian(mean140_d, axis = 1)
plt.plot(ver[step:], real_sammeans140_d, c = 'k')
plt.plot(ver[step:], real_sammeans140_d + np.nanstd(mean140_d, axis =1)/size140_d, linestyle = '--', c = 'k')
plt.plot(ver[step:], real_sammeans140_d - np.nanstd(mean140_d, axis =1)/size140_d, linestyle = '--', c = 'k')
plt.title('F160W')
plt.xlabel('p_e')
plt.ylabel('counts/arcs2')
plt.savefig('%s_dsky_vs_pe_%s.png'%(Target,step))
plt.close()


plt.figure()
plt.subplot(1,2,1)
ax = plt.gca()
plt.plot(ver, real_sammeans105, c = 'r')
ax.fill_between(ver,real_sammeans105 + np.nanstd(mean105, axis =1)/size105, real_sammeans105 - np.nanstd(mean105, axis =1)/size105, alpha = 0.5, color = 'r')
plt.plot(ver[step:], real_sammeans105_d, c = 'b')
ax.fill_between(ver[step:],real_sammeans105_d + np.nanstd(mean105_d, axis =1)/size105_d, real_sammeans105_d - np.nanstd(mean105_d, axis =1)/size105_d, alpha = 0.5, color = 'b')
plt.title('F105W')
plt.xlabel('p_e')
plt.ylabel('counts/arcs2')
plt.subplot(1,2,2)
ax = plt.gca()
plt.plot(ver, real_sammeans140, c = 'r')
ax.fill_between(ver,real_sammeans140 + np.nanstd(mean140, axis =1)/size140, real_sammeans140 - np.nanstd(mean140, axis =1)/size140, alpha = 0.5, color = 'r')
plt.plot(ver[step:], real_sammeans140_d, c = 'b')
ax.fill_between(ver[step:],real_sammeans140_d + np.nanstd(mean140_d, axis =1)/size140_d, real_sammeans140_d - np.nanstd(mean140_d, axis =1)/size140_d, alpha = 0.5, color = 'b')
plt.title('F160W')
plt.xlabel('p_e')
plt.ylabel('counts/arcs2')
plt.savefig('%s_dsky_vs_sky_%i.png'%(Target,step))
plt.close()
