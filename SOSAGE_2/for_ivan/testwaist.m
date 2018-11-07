clear;
clf;
%%%% beam total energy in MeV 
E=15
% derived quantities
gamma=E/0.5109+1
beta=sqrt(1-1/gamma^2);
bg=beta*gamma

%%%% Number of macroparticle
N=1000
cut=4
%%%% distance to waist
L=1e-2;
%%%% transverse emittance
nemitx=1e-6
nemity=1e-6
%%%% waist size
sigxL=1e-4
sigyL=1e-4

%%%% longitudinal parameters
chirp  = 0
sigz   = 1e-8
% estimate of longitudinal emittance
sigdpp = 1e-4
nemitz = bg*sigdpp*sigz



[X,Y,Z,PX,PY,PZ]=  Gen_3D_GaussianW(L, sigxL,nemitx,sigyL,nemity,chirp,sigz,nemitz,bg,N,cut,'T');


% check waist location

zz=[0:0.01:1]*2*L;

for i=1:length(zz)
   EnvX(i)=sqrt(var(X+zz(i)*PX./PZ));
end;

 plot (zz,EnvX,'o')
 
 

