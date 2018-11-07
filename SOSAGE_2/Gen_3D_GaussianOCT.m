function [X,Y,Z,PX,PY,PZ]=  Gen_3D_Gaussian(alphax,betax,nemitx,alphay,betay,nemity,chirp,sigz,nemitz,bg,N)
% function Gen_3D_Gaussian modified to generate distributions (6-D Gaussian) described by
% the parmaters in the three plane (x,y,z) and the energy for the beam,
%
% for the x plane needed are: emitx, alphax, betax,
% for the y plance needed paramters are: emity, alphay, betay,
% for the z plane needed parameters are: emitz, sigma_z, chirp
% also the energy for the beam is needed, hence bg 
% provide the emittance values in rms in SI units 
% for example 
%nemitx = 10e-6;		% rms spot size in x
%nemity = 10e-6;		% rms spot size in y
%nemitz = 3.0e-6;		% rms spot size in z
addpath ('/work/Tools/Generator_M/SOSAGE_2')


% set initial parameters   

 
q 	= 0.1; % charge is 100 pc 
stat 	= 1;
flag 	= -1;

c 	= 299792458;
m_ec2	= 0.5109e6;
%Laser_pulse = 0.616*1e-12; % in seconds , this is also used to get sigz if needed
% hence someone can use
%sigz   = Laser_pulse*c;  	% rms bunch length in z 

%--------------------------------------------------------------------------
offsetX = 0;
offsetY = 0;
%--------------------------------------------------------------------------


	% number of macroparticles to be generated

%---------------------------------------------------------------------------
%           main program
%


%--------------------------------------------------------------------------


CutG = 4
% x-plane information :

gammax = (1+alphax^2)/betax;

sigx   = sqrt(betax*nemitx/bg);				 	% rms spot size in x
sigxp  = sqrt(gammax*nemitx/bg);
			   	% divergence  in x
Cxxp   =  -alphax/(sqrt(1+alphax^2));                                 % this is what is known as r ( check note)

%---------------------------------------------------------------------------

% y-plane information :

gammay = (1+alphay^2)/betay;

sigy   = sqrt(betay*nemity/bg);				 	% rms spot size in y
sigyp  = sqrt(gammay*nemity/bg);			   	% divergence  in y
Cyyp   = -alphay/(sqrt(1+alphay^2));                               % this is r

%---------------------------------------------------------------------------

% z-plane information :


betaz  = (sigz^2*bg)/(nemitz);  	
C      = chirp*sigz^2; % this is <z z'>


alphaz =-(C*bg)/(nemitz);
gammaz = (1+alphaz^2)/betaz;

sigzp = sqrt(gammaz*nemitz/bg);
Czzp  =  -alphaz/(sqrt(1+alphaz^2));   
%---------------------------------------------------------------------------



% generate the X Y Z distribution.

%dump1 = sobseqm2(-1,N);
% X-direction
[X,XP]  = BiGaussianRPc (N, sigx, sigxp, Cxxp, CutG);
%sqrt(det(cov(X,XP)))
% Y-direction
[Y,YP]  = BiGaussianRPc (N, sigy, sigyp, Cyyp, CutG);
% four-dimension
%[X,XP,Y,YP]  = Gaussian4D (N, sigx, sigxp, Cxxp, sigy, sigyp, Cyyp);

% Z-direction
[Z,DPP] = BiGaussianRPc (N, sigz, sigzp, Czzp, CutG);

%---------------------------------------------------------------------------


% cook the XY coordinate 
%[X,Y] = Offset(X,Y,offsetX,offsetY,2);
%[X,Y]= RotateXY(X,Y,angle,N);
% scale such that  each bin =~ 40 micron
%[X,Y] = ScaleXY(X,Y,1.00);

%[X,Y]= Rescale (X,Y, Myrmsx,Myrmsy);


q=q/N; % charge per macroparticle


% Generate the momentum 
%PX= XP*bg.*(1+DPP);
%PY= YP*bg.*(1+DPP);
%PZ=    bg.*(1+DPP);

PX= XP*bg.*(1+DPP).*m_ec2/c;
PY= YP*bg.*(1+DPP).*m_ec2/c;
PZ= bg.*(1+DPP).*m_ec2/c;

%display beam parameters
Beam_Param(X,Y,Z/c,PX*m_ec2/c,PY*m_ec2/c,PZ*m_ec2/c);
%Beam_Param(X,Y,Z/c,PX,PY,PZ);
X=X-mean(X); Y=Y-mean(Y); Z=Z-mean(Z);
PX=PX-mean(PX); PY=PY-mean(PY); 
%PZ=PZ-mean(PZ);

 
BetaGamma = sqrt(PX.^2+PY.^2+PZ.^2);

Gamma = sqrt(1+BetaGamma.^2);

