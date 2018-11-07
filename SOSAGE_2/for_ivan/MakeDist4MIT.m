#!/usr/bin/octave -qf

arg_list=argv();

for i=1:nargin
  printf (' %s\n',arg_list{i});
  argum(i)=str2num(arg_list{i})
end;

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
%sigz   = 1e-8
sigz   = argum(1)
% estimate of longitudinal emittance
sigdpp = 1e-4
nemitz = bg*sigdpp*sigz



[X,Y,Z,PX,PY,PZ]=  Gen_3D_GaussianW(0, sigxL,nemitx,sigyL,nemity,chirp,sigz,nemitz,bg,N,cut,'T');
X=X-L*PX./PZ;
Y=Y-L*PY./PZ;

% the X,Y,Z,PX,PY,PZ available here should be save and use in the MIT code.


%plot (X,PX, XI, PX);


 
 

