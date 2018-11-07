
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

nemitx = 5e-6;		% rms spot size in xp
nemity = 5e-6;		% rms spot size in yp

N      = 200000;	% number of macroparticles to be generated

%---------------------------------------------------------------------------
%           main program
%
bg     = 80;

% SC ON
betax  =  1.033547e+00   
betay  =  2.202231e+00

% SC OFF 
%betax  = 7.130270e-01   
%betay  =  1.546395e+00

%6.671308e-01   
%betay  =  1.737559e+00 
%betax  = 5.060499e-01   
%betay  = 2.867341e+00
%betax  = 2.422917e-01   
%betay  = 1.330498e+00
alphax = 0.000000
alphay = 0.000000


sigx   = sqrt(nemitx/bg*betax);		% rms spot size in x
sigy   = sqrt(nemity/bg*betay);		% rms spot size in y
sigz   = 0.3e-3;  	% rms bunch length in z 
sigxp  = sqrt(nemitx/bg*(1+alphax^2)/betax);		% rms spot size in xp
sigyp  = sqrt(nemity/bg*(1+alphay^2)/betay);		% rms spot size in yp
sigdpp = 1e-2;;  					% rms dpp 


dump1 = sobseqm2(-1,N);
% X-direction
[X,XP] = BiGaussianRP (N, sigx, sigxp, 0);

% Y-direction
[Y,YP] = BiGaussianRP (N, sigy, sigyp, 0);

% Z-direction
[Z,DPP] = BiGaussianRP (N, sigz, sigdpp, 0);



% cook the XY coordinate 
[X,Y] = Offset(X,Y,offsetX,offsetY,2);
%[X,Y]= RotateXY(X,Y,angle,N);
% scale such that  each bin =~ 40 micron
[X,Y] = ScaleXY(X,Y,1.00);


q=q/N; % charge per macroparticle

% force the X and Y rms to New values
% ForceRMS(X,Y,1e-3, 1e-3);
%[X,Y] = ScaleXY(X,Y,1);
% below values for Q = 1nc
%Myrmsx = 0.50034*1e-3;
%Myrmsy = 0.50034*1e-3;
%
%[X,Y]= Rescale (X,Y, Myrmsx,Myrmsy);


PX= XP*bg.*(1+DPP)*m_ec2/c;
PY= YP*bg.*(1+DPP)*m_ec2/c;
PZ= bg.*(1+DPP)*m_ec2/c;

%display beam parameters
Beam_Param(X,Y,Z/c,PX,PY,PZ);

X=X-mean(X); Y=Y-mean(Y); Z=Z-mean(Z);
PX=PX-mean(PX); PY=PY-mean(PY); 
%PZ=PZ-mean(PZ);

%dump_ImpactT(X,Y,Z,PX,PY,PZ,1e-12);
dump_ImpactT_c(X,Y,Z/c,PX,PY,PZ,1e-12);
%dump_Parmela (X*100, PX/m_ec2, Y*100, PY/m_ec2, X*100, PZ/m_ec2,q);
% plot for image
%figure(2),subplot(2,1,1),plot(X(2:N).*1000,Y(2:N).*1000,'r.','markersize',0.5); grid on
%xlabel('x mm');ylabel('y mm');
%subplot(2,1,2),imagesc (Cropped_Image);


 

