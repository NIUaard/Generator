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
sigxL=1e-5
sigyL=1e-5

%%%% longitudinal parameters
chirp  = 0
sigz   = 1e-8
% estimate of longitudinal emittance
sigdpp = 1e-4
nemitz = bg*sigdpp*sigz

sigz=[0.01,0.02, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5, 6, 10];

for j=1:length(sigz)
    [X,Y,Z,PX,PY,PZ]=  Gen_3D_GaussianW(0, sigxL,nemitx,sigyL,nemity,chirp,sigz,nemitz,bg,N,cut,'T');
    X=X-L*PX./PZ;
    Y=Y-L*PY./PZ;
% here you should write a file compatible with AldebaranX
         
% the X,Y,Z,PX,PY,PZ available here should be sae and use in the MIT code.


%plot (X,PX, XI, PX);


% check waist location
%
zz=[0:0.01:1]*2*L;

for i=1:length(zz)
   EnvX(i)=sqrt(var(X+zz(i)*PX./PZ));
end;

plot (zz,EnvX,'o')
 
 

