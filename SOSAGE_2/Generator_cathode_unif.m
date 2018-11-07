
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
Xrms	= 60e-9;
Yrms	= Xrms;
offsetX = 0;
offsetY = 0;
angle 	= -0*pi/180;
%--------------------------------------------------------------------------

% input from the user

Ekin = 0.55;	% excess kinetic energy in eV (for thermal emittance)
sig_t = 0.1;	% laser pulse length in ps (laser assumed to be  Gaussian)
N = 500;	% number of macroparticles to be generated

%---------------------------------------------------------------------------
%           main program
%
% read and normalize to unity the VC image:
% generate sobol sequence for the 6D phase space 
% initialization
dump1 = sobseqm2(-1,N);
%P_ANGLE_SOB = sobseqm2(6,N);
P_ANGLE_SOB = rand(6,N);

% this how you would do a uniform distribution
u = P_ANGLE_SOB(3,:)'; v = P_ANGLE_SOB(4,:)';
X = sqrt(u)*2*Xrms.*cos(v*2*pi); % radius is 2 times projection rms value
Y = sqrt(u)*2*Xrms.*sin(v*2*pi);
w1 = P_ANGLE_SOB(5,:)'; w2 = P_ANGLE_SOB(6,:)';
Z = sig_t.*sqrt(-2.*log(w1)).*cos(2.*pi.*w2); % gaussian distribution 


% generate momentum
[PX, PY, PZ]= Thermal (N, Ekin);
%[PX, PY, PZ]= Cold (N, Ekin);

% generate temploral laser profile

meanPZ = mean(PZ);

q=q/N; % charge per macroparticle
disp('scale')
% force the X and Y rms to New values
% ForceRMS(X,Y,1e-3, 1e-3);
%[X,Y] = ScaleXY(X,Y,0.9);
% below values for Q = 1nc
disp('done')


%display beam parameters
Beam_Param(X,Y,0.*Z,PX,PY,PZ);

N=length(X);
[X,Y]= Rescale (X,Y, Xrms,Yrms);
%dump_ImpactT_c(X,Y,Z*1e-12,PX,PY,PZ,1e-12);
%dump_Parmela (X*100, PX/m_ec2, Y*100, PY/m_ec2, X*100, PZ/m_ec2,q);
% plot for image
%dump_ImpactT_v16(X,Y,Z*1e-12,PX,PY,PZ)
dump_ImpactT(X,Y,Z*1e-12,PX,PY,PZ)

subplot(1,1,1);
plot(X(2:N).*1000,Y(2:N).*1000,'r.','markersize',0.5); grid on
xlabel('x mm');ylabel('y mm');


 

