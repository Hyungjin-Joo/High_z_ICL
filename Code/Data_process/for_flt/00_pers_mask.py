import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

file_name = np.loadtxt('box_mask.dat', usecols = (0), dtype = 'str', unpack = True)
x1, y1, x2, y2, x3, y3, x4, y4 = np.loadtxt('box_mask.dat', usecols = (1,2,3,4,5,6,7,8), unpack = True)

if file_name.size>1:
    for i in range(file_name.size):
        print(file_name[i])
        sci = fits.open('%s'%file_name[i])
        sdata = sci[1].data
        xx,yy = np.shape(sdata)
        lx = np.linspace(0,xx,xx)
        ly = np.linspace(0,yy,yy)
        cx, cy = np.meshgrid(lx,ly)
        white = np.ones((xx,yy))
        xs = np.array([x1[i],x2[i],x3[i],x4[i]])
        ys = np.array([y1[i],y2[i],y3[i],y4[i]])
        index = np.array([0,1,2,3])
        for j in range(4):
            for k in range(4):
                if j < k:
                    slope = (ys[k]-ys[j])/(xs[k]-xs[j])
                    ycept = ys[k] - slope * xs[k]
                    ind_j = ~np.isin(index, j)
                    ind_k = ~np.isin(index, k)
                    ind = ind_j*ind_k
                    rem = index[ind]
                    yn_k = slope * xs[rem[0]] + ycept
                    yn_j = slope * xs[rem[1]] + ycept
                    if (yn_k-ys[rem[0]])*(yn_j-ys[rem[1]])>=0:
                        if ys[rem[0]] - yn_k >=0:
                            bool = cy>=(slope * cx + ycept)
                            white *= bool
                        if ys[rem[0]] - yn_k <=0:
                            bool = cy<=(slope * cx + ycept)
                            white *= bool

        sdata = sci[3].data
        sdata[white==1] = 4
        sci[3].data = sdata
        sci.writeto('%s'%file_name[i], overwrite = True)

else:
    print(file_name)
    sci = fits.open('%s'%file_name)
    sdata = sci[1].data
    xx,yy = np.shape(sdata)
    lx = np.linspace(0,xx,xx)
    ly = np.linspace(0,yy,yy)
    cx, cy = np.meshgrid(lx,ly)
    white = np.ones((xx,yy))
    xs = np.array([x1,x2,x3,x4])
    ys = np.array([y1,y2,y3,y4])
    index = np.array([0,1,2,3])
    for j in range(4):
        for k in range(4):
            if j < k:
                slope = (ys[k]-ys[j])/(xs[k]-xs[j])
                ycept = ys[k] - slope * xs[k]
                ind_j = ~np.isin(index, j)
                ind_k = ~np.isin(index, k)
                ind = ind_j*ind_k
                rem = index[ind]
                yn_k = slope * xs[rem[0]] + ycept
                yn_j = slope * xs[rem[1]] + ycept
                if (yn_k-ys[rem[0]])*(yn_j-ys[rem[1]])>=0:
                    if ys[rem[0]] - yn_k >=0:
                        bool = cy>=(slope * cx + ycept)
                        white *= bool
                    if ys[rem[0]] - yn_k <=0:
                        bool = cy<=(slope * cx + ycept)
                        white *= bool

    sdata = sci[3].data
    sdata[white==1] = 4
    sci[3].data = sdata
    sci.writeto('%s'%file_name, overwrite = True)
