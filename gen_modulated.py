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
cms 	=  299792458.
m_ec2	=  0.5109e6

#--------------------------------------------------------------------------
Xrms	=  1e-6
Yrms	=  Xrms
offsetX =  0.
offsetY =  0.
angle 	= -0*np.pi/180.
#--------------------------------------------------------------------------

# input from the user

Ekin  = 50e6	# excess kinetic energy in eV (for thermal emittance)
N = 100000	# number of macroparticles to be generated

#---------------------------------------------------------------------------
#           main program
#

P_ANGLE_SOB = rd.uniform(0,1,(7,N))

# this how you would do a uniform distribution
u = P_ANGLE_SOB[3,:]
v = P_ANGLE_SOB[4,:]
X = np.sqrt(u)*2*Xrms*np.cos(v*2*np.pi) # radius is 2 times projection rms value
Y = np.sqrt(u)*2*Xrms*np.sin(v*2*np.pi)

w   = 1e-4
lamb= 1e-5 
modu= 0.4
n   = 6.0
f= lambda z: np.exp(-2*(z/w)**n)*(1.+modu*np.cos(np.pi*2*z/lamb))

z=np.linspace(-2*w,2*w,1000)
F=f(z)
Norm=sum(F)*(z[1]-z[0])
F=F/Norm

u1 = (rd.rand(900000)*2.-1)*2*w  # uniform random samples scaled out
u2 =  rd.rand(900000)*2*max(F)    # uniform random samples
idx=np.where(u2<=f(u1)/Norm)[0] # rejection criterion
Z = u1[idx]



# generate momentum
PX, PY, PZ= gt.momt_cold(N, Ekin)

# Z is in ps so we convert to second (*1e-12) while passing to dump_ImpactT
gt.dump_Elegant(X,Y,Z/cms,PX,PY,PZ)
plt.plot (X,Y,'o')
plt.figure()
plt.plot (PX, PZ,'o')
plt.show()

 

