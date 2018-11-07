
% gen_distrib will generate a distribution for the 6D phase to be an 
% input  for Astra or for Parmela
% the code used a quasi random code the sobol sequence sobol.m was taken
% from http://www.math.uic.edu/~hanson/mcs507/cp4f04.html.


clear all; close all; 


% set initial parameters
 
q 	= 1; % charge is 1 nc 
stat 	= 1;
flag 	= -1;
c 	= 299792458;
m_ec2	= 0.5109e6;

%--------------------------------------------------------------------------
offsetX = 0;
offsetY = 0;
%--------------------------------------------------------------------------

nemity = 3e-6;		% rms spot size in yp
nemitz = 10e-6;		% rms spot size in dpp
%Ekin   = 15e6;		% kinetic energy in eV (for thermal emittance)
gam  = 30.00;
Ekin=(gam-1)*m_ec2
N      = 500;	% number of macroparticles to be generated

%---------------------------------------------------------------------------
%           main program
%
%gam = 1 + Ekin/m_ec2;
bg  = sqrt(gam^2 - 1 );

alphax = 0
betax  = 10
nemitx = 3e-6;		
alphay = 0
betay  = 10
nemity = 3e-6;		
sigz   = 1e-3
chirp  = 10  % in meter^{-1}
nemitz = 10e-6

gammax = (1+alphax^2)/betax
sigx   = sqrt(betax*nemitx/bg)
sigxp  = sqrt(gammax*nemitx/bg)
Coupxxp= -alphax/(betax*gammax)

gammay = (1+alphay^2)/betay
sigy   = sqrt(betay*nemity/bg)
sigyp  = sqrt(gammay*nemity/bg)
Coupyyp= -alphay/(betay*gammay)


betaz  = sigz^2/(nemitz/bg);  	
Corrzzp= chirp*sigz^2
alphaz =-Corrzzp/(nemitz/bg)
gammaz = (1+alphaz^2)/betaz
sigdpp = sqrt(gammaz*nemitz/bg)
Coupzzp=-alphaz/(betaz*gammaz)

% X-direction
[X,XP]  = BiGaussianRP (N, sigx, sigxp,  Coupxxp);
[Y,YP]  = BiGaussianRP (N, sigy, sigyp,  Coupyyp);
% Z-direction
[Z,DPP] = BiGaussianRP (N, sigz, sigdpp, Coupzzp);

%XP=XP;
%YP=YP;
%DPP=DPP;

% cook the XY coordinate 
[X,Y] = Offset(X,Y,offsetX,offsetY,2);
%[X,Y]= RotateXY(X,Y,angle,N);
% scale such that  each bin =~ 40 micron
[X,Y] = ScaleXY(X,Y,1.00);


q=q/N; % charge per macroparticle

% force the X and Y rms to New values
% below values for Q = 1nc
%Myrmsx = 0.50034*1e-3;
%Myrmsy = 0.50034*1e-3;
%
%[X,Y]= Rescale (X,Y, Myrmsx,Myrmsy);


PX= XP*bg.*(1+DPP)*m_ec2/c;
PY= YP*bg.*(1+DPP)*m_ec2/c;
PZ=    bg.*(1+DPP)*m_ec2/c;

%display beam parameters
Beam_Param(X,Y,Z,PX,PY,PZ);

X=X-mean(X); Y=Y-mean(Y); Z=Z-mean(Z);
PX=PX-mean(PX); PY=PY-mean(PY); 

dump_ImpactT(X,Y,Z/c,PX,PY,PZ,1e-12);
%dump_ImpactT_c(X,Y,Z*1e-12,PX,PY,PZ,1e-12);
%dump_Parmela (X*100, PX/m_ec2, Y*100, PY/m_ec2, X*100, PZ/m_ec2,q);
% plot for image
%figure(2),subplot(2,1,1),plot(X(2:N).*1000,Y(2:N).*1000,'r.','markersize',0.5); grid on
%xlabel('x mm');ylabel('y mm');
%subplot(2,1,2),imagesc (Cropped_Image);


 

