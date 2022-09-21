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

Nb       = 1000   # number of macroparticles per beamlets
plot_it = 1 
# all space/time unit should be SI:
cathode_spot_radius = 1e-3
beamlet_pitch       = 4e-3
beamlet_size        = 5     # beamlet_size x beamlet_size array
sigma_t_laser       = 3.4e-12 
laser_sigma_cut     = 3

# energies are in eV's
Excess_kinetic_eV   = 0.75
#Excess_kinetic_eV   = 48e6

ind=0

for i in range(beamlet_size):
   for j in range(beamlet_size):
       print(("beamlet #", ind))
       xcurrent=(-(beamlet_size-1)/2.+i)*beamlet_pitch 
       ycurrent=(-(beamlet_size-1)/2.+j)*beamlet_pitch 
       print((xcurrent, ycurrent))
       Xtmp, Ytmp= gt.tran_rad_unif (Nb, cathode_spot_radius)
       Xtmp=Xtmp+xcurrent
       Ytmp=Ytmp+ycurrent
       if ind==0:
           X=Xtmp
           Y=Ytmp
       else:
           X=np.append(X,Xtmp) 
           Y=np.append(Y,Ytmp) 
       ind+=1

N=len(X)
T   = gt.gauss_1d_cut  (N, sigma_t_laser, laser_sigma_cut) 
# cold = no thermal emittance
PX, PY, PZ = gt.momt_cold (N, Excess_kinetic_eV) 


print((len(X), len(Y), len(T), len(PX), len(PY), len(PZ)))


gt.dump_ImpactT_cathode(X,Y,T,PX,PY,PZ)


if (plot_it==1):
   plt.subplot (2,2,1)
   plt.plot (X,Y,'.')

   plt.subplot (2,2,2)
   plt.hist(T, 101, density=True, lw=3, edgecolor='red', facecolor="None")

   plt.subplot (2,2,3)
#   plt.plot (PX,PY,'.')
   plt.plot (X,PX,'.')
 
   plt.subplot (2,2,4)
   plt.plot (T,PZ,'.')

   plt.show()

 

