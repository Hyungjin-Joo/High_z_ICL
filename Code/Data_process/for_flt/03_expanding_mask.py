import ctypes as c
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import multiprocessing as mp
import time as t

# Functions
def multi(func, loops, ncore=8):
    '''
    Assign tasks to multiprocessing cores
    Parameters
    ----------
    func :  your func should have "only one" input argument.
    loops:  loop array
    ncore:  it must be lower than (maximum-1)
    '''
    pool = mp.Pool(ncore)
    pool.map(func,loops)
    pool.close()
    pool.join

def sum(x):
    global y, white, dseg
    print(x)
    for j in range(y):
        dseg[i,j] += white[i][j]

def cal(x):
    global k, y, rp, dseg, white
    for j in range(y):
        if dseg[x,j] !=0:
            eval = dseg[x,j]
            rt = rp[eval-1] - k
            if rt > 0:
                if x-1 > 0:
                    if dseg[x-1,j]==0:
                        if white_arr_2d[x-1,j] <= eval:
                            white_arr_2d[x-1,j] = eval
                if x+1 < 1014:
                    if dseg[x+1,j]==0:
                        if white_arr_2d[x+1,j] <= eval:
                            white_arr_2d[x+1,j] = eval
                if j-1 > 0:
                    if dseg[x,j-1]==0:
                        if white_arr_2d[x,j-1] <= eval:
                            white_arr_2d[x,j-1] = eval
                if j+1 < 1014:
                    if dseg[x][j+1]==0:
                        if white_arr_2d[x,j+1] <= eval:
                            white_arr_2d[x,j+1] = eval

files = input('File name?')
fname, type = files.split('.')
P_rs = np.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12.0])
for l in range(P_rs.size):
    P_r = P_rs[l]
    if P_r == 1.0:
        seg = fits.open('%s_check.fits'%fname)
    else:
        seg = fits.open('%s_check_%1.1f.fits'%(fname,P_r-0.1))
    r = np.loadtxt('%s.cat'%fname, usecols = (15), unpack = True)

    r[r > 20] = 20
    r[r < 0] = 0
    dseg = seg[1].data

    rp = np.int64(r*P_r)
    if P_r !=1.0:
        drp = np.int64(r*(P_r-0.1))
        rp -= drp
    x, y = dseg.shape
    for k in range(np.amax(rp)):
        print('k=',k,'/',np.amax(rp),'                    ', end = '\r')
        manager = mp.Manager()
        white = mp.Array(c.c_double, x*y)
        white_arr_1d = np.frombuffer(white.get_obj())
        white_arr_2d = white_arr_1d.reshape(x,y)
        multi(cal, np.arange(np.int64(x)), ncore = 4)
        dseg += np.int64(white_arr_2d)
    seg[1].data = dseg
    seg.writeto('%s_check_%1.1f.fits'%(fname,P_r), overwrite = True)




