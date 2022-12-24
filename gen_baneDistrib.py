import os
import numpy as np
from scipy import constants
from scipy.stats import halfnorm, norm
import matplotlib.pyplot as plt 
import generatorTool as gt 


def doorsteppdfsmooth2(x, lo, up, k, zdoor, sigma_head, sigma_tail, blength):
    mynorm = (6 * blength + k**2 * np.cos(k * zdoor) * (blength - zdoor)**3 + 3 * k * np.sin(k * zdoor) * (blength - zdoor)**2)/6
    tail_value = 1 - k*zdoor*np.sin(k*zdoor) + \
              k**2*zdoor**2/2*np.cos(k*zdoor) + \
              blength * (k*np.sin(k*zdoor) - k**2 * zdoor*np.cos(k*zdoor)) + \
              k**2/2*np.cos(k*zdoor)*blength**2
    
    head_value = 1
    
    condition = [x < 0, (x >= 0) & (x <= zdoor), (x > zdoor) & (x <= blength), x > blength]
    choice = [head_value * np.exp(-(x)**2 / 2 /sigma_head**2),
              1, 1 - k*zdoor*np.sin(k*zdoor) + \
              k**2*zdoor**2/2*np.cos(k*zdoor) + \
              x * (k*np.sin(k*zdoor) - k**2 * zdoor*np.cos(k*zdoor)) + \
              k**2/2*np.cos(k*zdoor)*x**2, tail_value * np.exp(-(x - blength)**2 / 2 /sigma_tail**2)]
    
    mynorm = mynorm + np.sqrt(np.pi / 2) * sigma_tail * tail_value + np.sqrt(np.pi / 2) * sigma_head * head_value

    return np.select(condition, choice)/mynorm


# example of parameters
global needed, lo, up, zdoor, sigma_head, sigma_tail, blength, freq, wavelength, k

needed=100000
lo = -0.1
up =  1.3
zdoor= 0.25
sigma_head = 0.01
sigma_tail = 0.01
blength = 1 
freq = 1.47970272e+11  # frequency of the mode 
wavelength = 299792458 / freq
k = 2 * np.pi 

alphax=0
alphay=0
betax=10
betay=10
normEmit=10e-6
gamma = 80
chirp=0 
incDpp = 1e-2 
emitGeom= normEmit/80

print ('beam size X:', np.sqrt(betax*emitGeom))
print ('beam size Y:', np.sqrt(betay*emitGeom))

def doorstepodfset(x):   
   return (doorsteppdfsmooth2(x, lo, up, k, zdoor, sigma_head, sigma_tail, blength))


# longitudinal distribution:
z    = wavelength*gt.mc_1D(needed, doorstepodfset, -0.2, 1.4, maxFunc=2)
# transverse phase space: 
x, xp, y, yp = gt.gaussian_phase_space_2dof (needed, alphax, betax, emitGeom, alphay, betay, emitGeom, Cut=4)
# energy spread
zref=np.max(z)
dgam = gamma*(chirp*(z-zref)+gt.gauss_1d_cut(needed, incDpp, Cut=4))
gam = gamma+dgam 

print(np.shape(z))

plt.subplot (2,2,1)
plt.hist(z,251)
plt.ylabel ('population')
plt.xlabel ('$\zeta$ (m)')
plt.subplot (2,2,2)
plt.plot(z,gam,'.')
plt.xlabel ('$\zeta$ (mm)')
plt.ylabel ('$\gamma$')
plt.subplot (2,2,3)
plt.plot(x, y,'.')
plt.xlabel ('$x$ (m)')
plt.ylabel ('$y$ (m)')
plt.subplot (2,2,4)
plt.plot(x, xp,'.')
plt.xlabel ('$x$ (m)')
plt.ylabel ('$x_[$ ()')
plt.tight_layout()
plt.show()
