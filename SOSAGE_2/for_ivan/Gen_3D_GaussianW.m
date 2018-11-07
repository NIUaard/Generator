

function [X,Y,Z,PX,PY,PZ]=  Gen_3D_GaussianW(L, sigxL,nemitx,sigyL,nemity,chirp,sigz,nemitz,bg,N,cut,force)
%function [X,Y,Z,PX,PY,PZ]=  Gen_3D_Gaussian(alphax,betax,nemitx,alphay,betay,nemity,chirp,sigz,nemitz,bg,N)
% function Gen_3D_Gaussian modified to generate distributions (6-D Gaussian) described by
% the parmaters in the three plane (x,y,z) and the energy for the beam,
%
% Generate a 3D Gaussian bunch such that it gets focused at a distance L from the location of generation
% with corresponding waists signxL and sigyL

% Force is a flag
% Force = 'T' force the rms to be correct
%

% set initial parameters   

 
q 	= 0.1; % charge is 100 pc 
stat 	=  1;
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

% required beta and alpha

% x plane:
% beta at the waist
betaxS = sigxL^2/(nemitx/bg)
% required beta and alpha at the generation point
betax  = betaxS+L^2/betaxS
alphax = L/betaxS

% y plane:
% beta at the waist
betayS = sigyL^2/(nemity/bg)
% required beta and alpha at the generation point
betay  = betayS+L^2/betayS
alphay = L/betayS



%---------------------------------------------------------------------------
%           main program
%


%--------------------------------------------------------------------------

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

% X-direction
[X,XP]  = BiGaussianRPc (N, sigx, sigxp, Cxxp,cut);
%sqrt(det(cov(X,XP)))
% Y-direction
[Y,YP]  = BiGaussianRPc (N, sigy, sigyp, Cyyp,cut);

% Z-direction
[Z,DPP] = BiGaussianRPc (N, sigz, sigzp, Czzp,cut);

%---------------------------------------------------------------------------


% cook the XY coordinate 
%[X,Y] = Offset(X,Y,offsetX,offsetY,2);
%[X,Y]= RotateXY(X,Y,angle,N);
% scale such that  each bin =~ 40 micron
%[X,Y] = ScaleXY(X,Y,1.00);

if (force=='T')
 [X,Y]= Rescale (X,Y, sigx,sigy);
end;


q=q/N; % charge per macroparticle


% Generate the momentum 

PX= XP*bg.*(1+DPP).*m_ec2/c;
PY= YP*bg.*(1+DPP).*m_ec2/c;
PZ=    bg.*(1+DPP).*m_ec2/c;

%display beam parameters
%Beam_Param(X,Y,Z/c,PX*m_ec2/c,PY*m_ec2/c,PZ*m_ec2/c);
%Beam_Param(X,Y,Z/c,PX,PY,PZ);
X=X-mean(X); Y=Y-mean(Y); Z=Z-mean(Z);
PX=PX-mean(PX); PY=PY-mean(PY); 
%PZ=PZ-mean(PZ);

 

