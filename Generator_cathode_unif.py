'''
 python generator based on Generator_cathode_unif.m 
'''

import numpy as np
import numpy.random as rd 
import matplotlib.pyplot as plt
import generatorTool as gt 
# set initial parameters
 
q 	=  1. # charge is 1 nc 
stat 	=  1.
flag 	= -1
c 	=  299792458.
m_ec2	=  0.5109e6

#--------------------------------------------------------------------------
Xrms	=  500e-6
Yrms	=  Xrms
offsetX =  0.
offsetY =  0.
angle 	= -0*np.pi/180.
#--------------------------------------------------------------------------

# input from the user

Ekin  = 0.5	# excess kinetic energy in eV (for thermal emittance)
sig_t = 3.0	# laser pulse length in ps (laser assumed to be  Gaussian)
N = 10000	# number of macroparticles to be generated

#---------------------------------------------------------------------------
#           main program
#

P_ANGLE_SOB = rd.uniform(0,1,(7,N))

# this how you would do a uniform distribution
u = P_ANGLE_SOB[3,:]
v = P_ANGLE_SOB[4,:]
X = np.sqrt(u)*2*Xrms*np.cos(v*2*np.pi) # radius is 2 times projection rms value
Y = np.sqrt(u)*2*Xrms*np.sin(v*2*np.pi)
w1 = P_ANGLE_SOB[5,:] 
w2 = P_ANGLE_SOB[6,:]
Z = sig_t*np.sqrt(-2*np.log(w1))*np.cos(2.*np.pi*w2)  # gaussian distribution 


# generate momentum
PX, PY, PZ= gt.Cold (N, Ekin)

# Z is in ps so we convert to second (*1e-12) while passing to dump_ImpactT
gt.dump_ImpactT_c(X,Y,Z*1e-12,PX,PY,PZ)
plt.plot (X,Y,'o')
plt.figure()
plt.plot (PX, PZ,'o')
plt.show()

 

