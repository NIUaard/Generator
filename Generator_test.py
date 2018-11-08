'''
 python generator based on Generator_cathode_unif.m 
'''

import numpy as np
import numpy.random as rd 
import matplotlib.pyplot as plt
import generatorTool as gt 
# set initial parameters

f2d=1
 
# X=  gt.Unif_2d_cart(1000)
# X=  gt.Unif_2d_rad(100000)
# X=  gt.Gauss_1d(10000)
# X=  gt.Gauss_1d_cut(100000,1)

#X=  gt.Gauss_2d_cart(1000,0.5)
X=  gt.Gauss_2d_cart_cut(100000,0.2,2.0)
print(np.shape(X))


if f2d==1: 
   plt.subplot (2,2,1)
   plt.plot (X[0,:], X[1,:],'.')
#   plt.ylim(-0.5,0.5)
#   plt.xlim(-0.5,0.5)
   plt.subplot (2,2,2)
   plt.hist(X[0,:], 101, normed=1, lw=3, edgecolor='red', facecolor="None")
   plt.subplot (2,2,3)
   plt.hist(X[1,:], 101, normed=1, lw=3, edgecolor='red', facecolor="None")
else:
   plt.subplot (2,1,1)
   plt.plot (X,'.')
   plt.subplot (2,1,2)
   plt.hist(X, 101, normed=1, lw=3, edgecolor='red', facecolor="None")
   


plt.show()

 

