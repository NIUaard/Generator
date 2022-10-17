'''
 python generator based on Generator_cathode_unif.m 
'''

from scipy.stats import qmc
import numpy as np
import numpy.random as rd 
import matplotlib.pyplot as plt
import generatorTool as gt 
import mpl_scatter_density
cMap = plt.cm.get_cmap('gist_earth_r').copy()
cMap.set_under('white')


cMap = plt.cm.get_cmap('gist_earth_r').copy()
cMap.set_under('white')

def myfun(x):
   return(np.sqrt(x[0,:]**2+x[1,:]**2)**nxy+np.abs(x[2,:])**nz)
 


needed=100000

rmsxy=1.93e-3 
rmst=2.24e-12  
scale=1.5

nxy_arr= np.array([0.5, 1, 1.5, 2., 4, 8])
nz= 2   #4 

PX, PY, PZ = gt.momt_thermal_iso(needed, Ekin=0.55)


fig = plt.figure()

for i in range(len(nxy_arr)):
   nxy=nxy_arr[i]
   fnametmp='partcl_superellips_'+str(nxy)+'_'+str(nz)+'.data'
   ell=gt.mc_3DBoundary(needed, myfun)
   print(">>> info", i, np.shape(ell))
   x=ell[0,:]/np.std(ell[0,:])*rmsxy
   y=ell[1,:]/np.std(ell[1,:])*rmsxy
   t=ell[2,:]/np.std(ell[2,:])*rmst
   print(">>> rms x, y, t::", np.std(x), np.std(y), np.std(t))
   print(">>> exponents nxy, nz", nxy, nz)

   exec ("ax1"+str(i)+" = fig.add_subplot(6, 2, "+str(2*i+1)+", projection='scatter_density')")
   exec ("ax2"+str(i)+" = fig.add_subplot(6, 2, "+str(2*i+2)+", projection='scatter_density')")

   exec ("ax1"+str(i)+".scatter_density(t, x, cmap=cMap)")
   exec ("ax1"+str(i)+".set(xlim=(-scale*t.max(),scale*t.max()),ylim=(-scale*x.max(),scale*x.max()))")

   exec ("ax2"+str(i)+".scatter_density(x,y, cmap=cMap)")
   exec ("ax2"+str(i)+".set(xlim=(-scale*x.max(),scale*x.max()),ylim=(-scale*y.max(),scale*y.max()))")
   gt.dump_ImpactT_cathode(x,y,t, PX, PY, PZ, fname=fnametmp)
   
#print ("ell:", np.shape(ell))
#plt.title ('nt='+str(nxy)+'  nz='+str(nz))
##plt.hexbin(ell[2,:],ell[0,:],cmap=cMap, gridsize=(51,51) )
#plt.subplot (2,2,1)
#plt.plot(t, x,'.')
#plt.subplot (2,2,2)
#plt.plot(t,y,'.')
#plt.subplot (2,2,3)
#plt.plot(x,y,'.')
plt.show()
