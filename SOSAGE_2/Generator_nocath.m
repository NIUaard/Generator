
% gen_distrib will generate a distribution for the 6D phase to be an 
% input  for Astra or for Parmela
% the code used a quasi random code the sobol sequence sobol.m was taken
% from http://www.math.uic.edu/~hanson/mcs507/cp4f04.html.


clear all; close all; 


% set initial parameters
 
q 	= 3.2; % charge is 1 nc 
stat 	= 1;
flag 	= -1;
c 	= 299792458;
m_ec2	= 0.5109e6;

%--------------------------------------------------------------------------
offsetX = 0;
offsetY = 0;
%--------------------------------------------------------------------------

nemitx = 10e-7;		% rms spot size in xp
nemity = 10e-7;		% rms spot size in yp
nemitz =  3e-7;		% rms spot size in dpp
Ekin   = 15e6;		% kinetic energy in eV (for thermal emittance)
N      = 71;	% number of macroparticles to be generated

%---------------------------------------------------------------------------
%           main program
%
gam = 1 + 15e6/m_ec2;
bg  = sqrt(gam^2 - 1 );

sigx   = sqrt(nemitx/bg);		% rms spot size in x
sigy   = sqrt(nemity/bg);		% rms spot size in y
sigz   = sqrt(nemitz/bg);  	% rms bunch length in z 
sigxp  = sqrt(nemitx/bg);		% rms spot size in x
sigyp  = sqrt(nemity/bg);		% rms spot size in y
sigdpp = sqrt(nemitz/bg);;  	% rms bunch length in z 


% X-direction
emitx= nemitx/(bg);
Corrxxp = 1/(sigx*sigxp)*sqrt(sigx^2*sigxp^2-emitx^2)
[X,XP] = BiGaussianRP (N, sigx, sigxp, Corrxxp);

% Y-direction
emity= nemity/(bg);
Corryyp = 1/(sigy*sigyp)*sqrt(sigy^2*sigyp^2-emity^2)
[Y,YP] = BiGaussianRP (N, sigy, sigyp, Corryyp);

% Z-direction
emitz= nemitz/(bg); % nemitz is in z-dgamma so emitz is in z-dpp
Corrzdpp = 1/(sigz*sigdpp)*sqrt(sigz^2*sigdpp^2-emitz^2)
[Z,DPP] = BiGaussianRP (N, sigz, sigdpp, Corrzdpp);



% cook the XY coordinate 
[X,Y] = Offset(X,Y,offsetX,offsetY,2);
%[X,Y]= RotateXY(X,Y,angle,N);
% scale such that  each bin =~ 40 micron
[X,Y] = ScaleXY(X,Y,1.00);


q=q/N; % charge per macroparticle

% force the X and Y rms to New values
% ForceRMS(X,Y,1e-3, 1e-3);
[X,Y] = ScaleXY(X,Y,0.5);
% below values for Q = 1nc
%Myrmsx = 0.50034*1e-3;
%Myrmsy = 0.50034*1e-3;
%
%[X,Y]= Rescale (X,Y, Myrmsx,Myrmsy);


PX= XP*bg.*(1+DPP)*m_ec2/c;
PY= YP*bg.*(1+DPP)*m_ec2/c;
PZ= bg.*(1+DPP)*m_ec2/c;

%display beam parameters
Beam_Param(X,Y,Z,PX,PY,PZ);

X=X-mean(X); Y=Y-mean(Y); Z=Z-mean(Z);
PX=PX-mean(PX); PY=PY-mean(PY); PZ=PZ-mean(PZ);

dump_ImpactT(X,Y,Z/c,PX,PY,PZ,1e-12);
%dump_ImpactT_c(X,Y,Z*1e-12,PX,PY,PZ,1e-12);
%dump_Parmela (X*100, PX/m_ec2, Y*100, PY/m_ec2, X*100, PZ/m_ec2,q);
% plot for image
%figure(2),subplot(2,1,1),plot(X(2:N).*1000,Y(2:N).*1000,'r.','markersize',0.5); grid on
%xlabel('x mm');ylabel('y mm');
%subplot(2,1,2),imagesc (Cropped_Image);


 

