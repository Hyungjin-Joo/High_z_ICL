import copy as cp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clip

def myfitting(fname):
    sci = fits.open(fname)
    sdata = sci[1].data
    
    split_fname = fname.split('.')
    new_fname = '_plf.'.join(split_fname)
    
    c_sdata = cp.deepcopy(sdata)
    mask_data = sigma_clip(c_sdata, sigma = 3, maxiters = None, cenfunc = np.mean, masked=True, copy=False)
    n_sdata = cp.deepcopy(c_sdata)
    n_sdata[mask_data.mask == True] = 1e+8
    
    x, y = c_sdata.shape
    mx, my = np.mgrid[:x, :y]
    
    xx = mx[n_sdata!=1e+8]
    yy = my[n_sdata!=1e+8]
    zz = sdata[n_sdata!=1e+8]
    asize = xx.size
    oo = np.ones(asize)
    
    ATA = np.array([[np.sum(xx*xx),np.sum(xx*yy),np.sum(xx)], [np.sum(xx*yy),np.sum(yy*yy),np.sum(yy)], [np.sum(xx),np.sum(yy),np.sum(oo)]])
    if np.trace(ATA) == 0:
        print('Warning! Tracer = 0!')
        return(sdata)
    ATAi = np.linalg.inv(ATA)
    
    B = np.array([np.sum(xx*zz),np.sum(yy*zz),np.sum(zz)])
    
    a = np.sum(ATAi[0]*B)
    b = np.sum(ATAi[1]*B)
    c = np.sum(ATAi[2]*B)
    
    xxx = mx
    yyy = my
    new = (a*xxx + b*yyy + c)
    
    c_sdata -= new

    sdata[1].data = c_sdata
    sci.writeto(new_fname, overwrite = True)
    
    sdata = c_sdata
    return(sdata)

fnames = np.loadtxt('flts.lis', dtype='str', unpack=True)

for i in range(fnames.size):
    print(i, '/', fnames.size, ',', fnames[i], end="\r")
    aaa = myfitting(fnames[i])
