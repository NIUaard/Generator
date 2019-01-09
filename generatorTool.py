import numpy as np
import numpy.random as rd 
import matplotlib.pyplot as plt
import sobol_lib as sob 

''''
Set of generic function for generating particle distributions for tracking code
created AUG-27-2017 P. Piot, NIU
- random generators are all based on Sobol pseudo-random generator

'''
m_ec2 = 0.5109e6
cms   = 299792458.
DEB = 1
############### elementary generating functions #############################

def Unif_1d (n, skip=None):
   '''
   generate a 1D uniform distribution in [0,1]
   '''
   global DEB
   if DEB==1:
      print '>>>> Unif_1d'
      
   if skip==None:
      skip=np.random.random_integers(0,1000)
      
   out=sob.i4_sobol_generate (1,n,skip)
   return(out.flatten())

def Gauss_1d (n, skip=None):
   '''
   generate a 1D Gaussian distribution in with a sigma of 1
   '''
   global DEB
   if DEB==1:
      print '>>>> Gauss_1d'
      
   if skip==None:
      skip=np.random.random_integers(0,1000)
   U=Unif_2d_cart (n)
   out= np.sqrt(-2.0*np.log(U[0,:]))*np.sin(2*np.pi*U[1,:]);
   return(out.flatten())

def Gauss_1d_cut (n, cut, skip=None):
   '''
   generate a 1D Gaussian distribution in with a sigma of 1 and cut in sigma unit
   '''
   global DEB
   if DEB==1:
      print '>>>> Gauss_1d_cut'
      
   out=np.zeros((n))
   needed = n
   ind_i  = 0 
   ind_f  = 0 
   while needed>0:
      if DEB==1: 
          print '<-start', needed, ind_i
      gen    = Gauss_1d (needed, skip)
      keep   = np.where(np.abs(gen)<=cut)
      Ngood  = len(gen[keep])
      ind_f  = ind_i + Ngood
      if DEB==1: 
         print '.. fill', Ngood, ind_i, ind_f
      if ind_f<=n: 
         if DEB==1: 
            print 'ind_f<n'
         out[ind_i:ind_f]=gen[keep]
      if ind_f>n:
         if DEB==1: 
            print 'ind_f>n'
         ind_f=n
	 temp=gen[keep]
	 out[ind_i:ind_f]=gen[keep[0:n-ind_i]]
	 
      needed = n-ind_f
      ind_i  = ind_f
	 
   return(out)

def Unif_2d_cart (n, skip=None):
   '''
   generates a uniform distribution in 2D within [0,1] 
   '''
   global DEB
   if DEB==1:
      print '>>>> Unif_2d_cart'
      
   if skip==None:
      skip=np.random.random_integers(0,1000)
   out=sob.i4_sobol_generate (2,n,skip)
   return(out)
   
def Unif_3d_cart (n, skip=None):
   '''
   generates a uniform distribution in 2D within [0,1] 
   '''
   global DEB
   if DEB==1:
      print '>>>> Unif_2d_cart'
      
   if skip==None:
      skip=np.random.random_integers(0,1000)
   out=sob.i4_sobol_generate (3,n,skip)
   return(out)
   
def Unif_5d (n, skip=None):
   '''
   generates a uniform distribution in 5D within [0,1] 
   '''
   global DEB
   if DEB==1:
      print '>>>> Unif_2d_cart'
      
   if skip==None:
      skip=np.random.random_integers(0,1000)
   out=sob.i4_sobol_generate (4,n,skip)
   return(out)
   
def Gauss_2d_cart (n, corr, skip=None):
   '''
   generates a 2D Gaussian distribution in x, y sigmas=1 in both directions
   the distribution with correlation corr
   '''
   global DEB
   if DEB==1:
      print '>>>> Gauss_2d_cart'
      
   U=Unif_2d_cart (n, skip)
   u=U[0,:]
   v=U[1,:]
   out=np.zeros((2,n))
   out[0,:]=np.sqrt(-2*np.log(u))*(np.sqrt(1-corr**2)*np.cos(2*np.pi*v)+corr*np.sin(2*np.pi*v))
   out[1,:]=np.sqrt(-2*np.log(u))*np.sin(2*np.pi*v)
   return(out)

   
def Unif_2d_rad (n, skip=None):
   '''
   generate a uniform distribution between r in [0,1] and phi in [0, 2pi]
   '''
   global DEB
   if DEB==1:
      print '>>>> Unif_2d_rad'
      
   U=Unif_2d_cart (n, skip)
   phi=U[0,:]*2*np.pi
   rad=U[1,:]
   out=np.zeros((2,len(phi)))
   out[0,:]=np.sqrt(rad)*np.cos(phi)
   out[1,:]=np.sqrt(rad)*np.sin(phi)
   return(out)
   
def Unif_2d_ellipse (n, r1, r2, phi0, skip=None):
   '''
   generate a uniform distribution between r in [0,1] and phi in [0, 2pi]
   '''
   global DEB
   if DEB==1:
      print '>>>> Unif_2d_ellipse'
      
   U=Unif_2d_cart (n, skip)
   phi=U[0,:]*2*np.pi
   rad=U[1,:]

   out =np.zeros((2,len(phi)))
   outt=np.zeros((2,len(phi)))
   outt[0,:]=r1*np.sqrt(rad)*np.cos(phi)
   outt[1,:]=r2*np.sqrt(rad)*np.sin(phi)
   
   out[0,:] =  outt[0,:]*np.cos(phi0) + outt[1,:]*np.sin(phi0)
   out[1,:] = -outt[0,:]*np.sin(phi0) + outt[1,:]*np.cos(phi0)
   return(out)
   
def Gauss_2d_cart_cut (n, corr, cut):
   '''
   generate a 2D Gaussian distribution in with a sigma of 1 and cut in sigma unit
   the cut applies in both directions
   '''
   global DEB
   if DEB==1:
      print '>>>> Gauss_2d_cart_cut'
      
   out=np.zeros((2,n))
   needed = n
   ind_i  = 0 
   ind_f  = 0 
   while needed>0:
      if DEB==1: 
          print '<-start', needed, ind_i
      gen    = Gauss_2d_cart (needed, corr)
      
      keep   = np.where((np.abs(gen[0,:])<=cut) & (np.abs(gen[1,:])<=cut))

      Ngood  = len(keep[0])
      ind_f  = ind_i + Ngood
      if DEB==1: 
         print '.. fill', Ngood, ind_i, ind_f
      if ind_f<=n: 
         if DEB==1: 
            print 'ind_f<n', keep
         out[0,ind_i:ind_f]=gen[0,keep]
	 out[1,ind_i:ind_f]=gen[1,keep]
      if ind_f>n:
         if DEB==1: 
            print 'ind_f>n'
         ind_f=n
	 temp=gen[keep]
	 out[0, ind_i:ind_f]=gen[0, keep[0:n-ind_i]]
	 out[1, ind_i:ind_f]=gen[1, keep[0:n-ind_i]]
	 
      needed = n-ind_f
      ind_i  = ind_f
	 
   return(out)

def mc_1D (n,func,minOrd, maxOrd):
   '''
   generate a distribution via a Monte-Carlo rejection 
   with density function given by func peak normalized to on 
   we assume max(func)=1
   - [minOrd, maxOrd] are the interval over which the monte-carlo is applied. 
   - func is the function define over [minOrd, maxOrd] note that maxOrd>minOrd
   '''

   global DEB

   if DEB==1:
      print '>>>> mc_1D'
      
   out=np.zeros((n))
   needed=n
   ind_i  = 0 
   ind_f  = 0 
   
   while needed>0:
      print "needed=", needed
      gen=Unif_2d_cart(needed)
      x=minOrd + gen[0,:]*(maxOrd-minOrd)
      y= gen[1,:]*1.2
      keep   = np.where(func(x)>gen[1,:])
      Ngood  = len(keep[0])
      ind_f  = ind_i + Ngood
      if DEB==1: 
         print '.. fill', Ngood, ind_i, ind_f
      if ind_f<=n: 
         if DEB==1: 
            print 'ind_f<n', keep
         out[ind_i:ind_f]=x[keep]
      if ind_f>n:
         if DEB==1: 
            print 'ind_f>n'
         ind_f=n
	 out[ind_i:ind_f]=x[keep[0:n-ind_i]]
	 
      needed = n-ind_f
      ind_i  = ind_f
	 
   return(out)
      
   
   
   
############### real-space generating functions #############################

def tran_rad_unif (N, Rad):
   '''
   generates a radial distribution in (x,y)
   - N:   number of macroparticles  
   - Rad: hard-edge radius  
   '''
   
   T=Rad*Unif_2d_rad(N)
   return (T[0,:], T[1,:])


def gauss_1d_cut (N, Sigma, Cut):
   '''
   generates a 1-d Gaussian distribution
   - N:      number of macroparticles  
   - Sigma: rms size
   - Cut:   cut in "number of sigma" unit
   '''
   
   X=Gauss_1d_cut (N, Cut)
   return (X*Sigma)


def momt_thermal_iso (N, Ekin):
   '''
   generates an isotropic momentum distributions:
   - N:   number of macroparticles  
   - Ekin: excess in kinetic energy (Ekin~h\nu-Phi), 
   '''
   global m_ec2
   global cms
   # calculate the thermal momentum
   p = np.sqrt(Ekin**2+2*Ekin*m_ec2) 
   
   print "total momentum", p 

   U=Unif_2d_cart (N)  
   u=U[0,:]
   v=U[1,:]

   the =2.*np.pi*u 
   phi = -np.pi/2+np.pi*v
#   costhe=(u-0.5)*2.
#   sinthe=np.sqrt(1-costhe**2)
#   sinphi=-v
#   cosphi=np.sqrt(1-sinphi**2)
# calculate momentum due to thermal emission
   PX = p*np.cos(the)*np.sin(phi) 
   PY = p*np.sin(the)*np.sin(phi) 
   PZ = p*np.cos(phi) 
#   PX = p*u 
#   PY = p*v 
#   PZ = p*v 
    
   return (PX, PY, PZ)
   
   
def momt_cold (N, Ekin):
   '''
   generates an isotropic momentum distributions:
   - N:   number of macroparticles  
   - Ekin: excess in kinetic energy (Ekin~h\nu-Phi), 
   '''
   
   global m_ec2
   global cms
   # calculate the thermal momentum
   p = np.sqrt(Ekin**2+2*Ekin*m_ec2) 
   
   print "total momentum", p 

   U=Unif_2d_cart (N)  
   u=0*U[0,:]
   v=0*U[1,:]

   the =2.*np.pi*u 
   phi =v

   PX = p*np.cos(the)*np.sin(phi) 
   PY = p*np.sin(the)*np.sin(phi) 
   PZ = p*np.cos(phi) 
    
   return (PX, PY, PZ)
   

def gaussian_phase_space(n, alpha, beta, emitgeom, cut):   
   '''
     generate a Gaussian phase space (x,x')with given 
     CS (alpha, beta) and emittance (emit) parameters. 
     cut: number of sigma for cut
   '''
   sigma_x  = np.sqrt(beta*emitgeom)
   sigma_xp = emitgeom/sigma_x  # uncorrelated momentum spread
   corr     = - alpha/beta
   out=Gauss_2d_cart_cut (n, corr, cut)
   out[0,:]=out[0,:]*sigma_x
   out[1,:]=out[1,:]*sigma_xp
   return(out[0,:], out[1,:])
   

############### dumping of particle distribution in files #############################

def dump_ImpactT_cathode(X,Y,Z,PX,PY,PZ):

   '''
    expect X[m], Y[m]. Z[sec]. PX[eV/c], PY[eV/c], PZ[eV/c]
   '''

   global m_ec2
   global cms

   TotalEmissionTimeSec = np.abs(np.max(Z)-np.min(Z))

   N       = len(X)
   #  in impact-T all the particle shoud start with z<0)
   aux     = 0.00
#   cdt     = cms*IntegStep
   # THIS is just to check Ek is the proper value
   MeanPz  = np.mean(PZ)
   MeanP   = np.mean (np.sqrt(PX**2+PY**2+PZ**2))
   MeanEk  = np.sqrt(MeanP**2 + m_ec2**2) -m_ec2

   bgx   = PX/m_ec2
   bgy   = PY/m_ec2
   bgz   = PZ/m_ec2
   bg    = np.sqrt (bgx**2 + bgy**2 + bgz**2)

   gamma = np.sqrt (bg**2+1.) 
#PP need to check the two line below when thermal emittance 
# is included sigmat was changing if I use betaz   
   betaz    = bg/gamma
# this was original choice   
#   betaz   = bgz/(gamma)
   Z2      = (Z-max(Z))*cms*np.mean(betaz)+1e-16;
   zMean   = np.mean(Z2)
#   Phi     = ( mean(Z-max(Z))-1e-16)*cms/((cms)/1.3e9*1.0/360.00)
   SigmZ   = np.std(Z2)
   bgzA    = np.mean (bgz)


   print "Total Emission time [sec]  =", np.max(Z)-np.min(Z)
   print "Mean gamma          [-]    =", np.mean(gamma)
   print "Mean beta           [-]    =",np.mean(betaz)
   print "Mean Kinetic Energy [eV]   =",MeanEk
   print "Mean momentm        [eV/c] =",MeanP 
   print "Sigma_z             [m]    =",SigmZ
   
   
   fid=open ('partcl.data','w')

   fid.write (str(N+1)+'\n')

   fid.write('{:13.5e}'.format(aux)+'{:13.5e}'.format(aux)+'{:13.5e}'.format(aux)+  \
             '{:13.5e}'.format(aux)+'{:13.5e}'.format(zMean)+'{:13.5e}'.format(bgzA)+'\n')
   for i in range(N):
        fid.write('{:13.5e}'.format(X[i])+'{:13.5e}'.format(bgx[i])+'{:13.5e}'.format(Y[i])+  \
            '{:13.5e}'.format(bgy[i])+'{:13.5e}'.format(Z2[i])+'{:13.5e}'.format(bgz[i])+'\n')

   fid.close()
   
   
def dump_Elegant(X,Y,Z,PX,PY,PZ):

   '''
    expect X[m], Y[m]. Z[sec]. PX[eV/c], PY[eV/c], PZ[eV/c]
    but convert to ELEGANT coordinate systems within this function 
   '''

   global m_ec2
   global cms

   N=len(X)
   
   f = open('bunch.sdds', 'w')
   f.write ("SDDS1\n")
   f.write ("&column name=x, units=m,   type=double, &end  \n")
   f.write ("&column name=y, units=m,    type=double, &end  \n")
   f.write ("&column name=t,   units=s,type=double, &end  \n")
   f.write ("&column name=xp,type=double, &end  \n")
   f.write ("&column name=yp,type=double, &end  \n")
   f.write ("&column name=p,  units=m$be$nc type=double, &end  \n")
   f.write ("&column name=particleID, type=long,  &end \n")
   f.write ("&data mode=ascii, &end          \n")
   f.write ("! page number 1                      \n")
   f.write (str(N)+"\n")

   x=X
   y=Y
   t=Z
   xp=PX/PZ
   yp=PY/PZ
   p =PZ/m_ec2
   for i in range(N):
      f.write(str(x[i])+"\t"+str(y[i])+"\t"+str(t[i])+"\t"+str(xp[i])+"\t"+str(yp[i])+ \
                  "\t"+str(p[i])+"\t"+str(i)+"\n")

   f.close()

