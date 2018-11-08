'''
 generates a radial uniform distribution with Gaussian temporal profile
 and isotropic momentum distribution
 save the distribution as an impact intput file (cathod distribution)
'''

import numpy as np
import numpy.random as rd 
import matplotlib.pyplot as plt
import generatorTool as gt 


# set initial parameters

N       = 1000   # number of macroparticles
plot_it = 1 
# all space/time unit should be SI:
cathode_spot_radius = 4e-3
sigma_t_laser       = 3.4e-12 
laser_sigma_cut     = 3

# energies are in eV's
Excess_kinetic_eV   = 0.75
#Excess_kinetic_eV   = 48e6


X, Y= gt.tran_rad_unif (N, cathode_spot_radius)
T   = gt.gauss_1d_cut  (N, sigma_t_laser, laser_sigma_cut) 
# cold = no thermal emittance
PX, PY, PZ = gt.momt_cold (N, Excess_kinetic_eV) 


print(len(X), len(Y), len(T), len(PX), len(PY), len(PZ))


gt.dump_ImpactT_cathode(X,Y,T,PX,PY,PZ)


if (plot_it==1):
   plt.subplot (2,2,1)
   plt.plot (X,Y,'.')

   plt.subplot (2,2,2)
   plt.hist(T, 101, normed=1, lw=3, edgecolor='red', facecolor="None")

   plt.subplot (2,2,3)
#   plt.plot (PX,PY,'.')
   plt.plot (X,PX,'.')
 
   plt.subplot (2,2,4)
   plt.plot (T,PZ,'.')

   plt.show()

 

