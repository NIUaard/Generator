
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
angle 	= -90*pi/180;
bin_size = 39e-6;
%--------------------------------------------------------------------------

% input from the user

% 
%fname = 'LFDneg10.bmp'
fname = 'create_image/uniform.png'
%fname = 'create_image/modulation_2_05.png'
Ekin = 15e6;	% excess kinetic energy in eV (for thermal emittance)
sig_t = 3.5;	% laser pulse length in ps (laser assumed to be  Gaussian)
N = 200000;	% number of macroparticles to be generated

%---------------------------------------------------------------------------
%           main program
%
% read and normalize to unity the VC image:
A = double(imread(fname));
A = A./max(max(A));

% select ROI
Cropped_Image=Crop_Image (A);


% generate sobol sequence for the 6D phase space 
% initialization
dump = sobseqm4(-1,N,Cropped_Image);
tic 
%generation
XY_SOB = sobseqm4(4,N,Cropped_Image); % change name is MakeXY.c
toc
X = XY_SOB(1,:)'; Y = XY_SOB(2,:)';

% cook the XY coordinate 
[X,Y] = Offset(X,Y,offsetX,offsetY,2);
%[X,Y]= RotateXY(X,Y,angle,N);
% scale such that  each bin =~ 40 micron
[X,Y] = ScaleXY(X,Y,bin_size);

% generate momentum
[PX, PY, PZ]= Thermal (N, Ekin);

% generate temploral laser profile
Z= Gaussian1 (N,sig_t);

meanPZ = mean(PZ);

q=q/N; % charge per macroparticle

% force the X and Y rms to New values
% ForceRMS(X,Y,1e-3, 1e-3);
[X,Y] = ScaleXY(X,Y,0.5);
% below values for Q = 1nc
Myrmsx = 0.50034*1e-3;
Myrmsy = 0.50034*1e-3;

[X,Y]= Rescale (X,Y, Myrmsx,Myrmsy);


%display beam parameters
Beam_Param(X,Y,0.*Z,PX,PY,PZ);

dump_ImpactT(X,Y,Z*1e-12,PX,PY,PZ,1e-12);
%dump_ImpactT_c(X,Y,Z*1e-12,PX,PY,PZ,1e-12);
%dump_Parmela (X*100, PX/m_ec2, Y*100, PY/m_ec2, X*100, PZ/m_ec2,q);
% plot for image
figure(2),subplot(2,1,1),plot(X(2:N).*1000,Y(2:N).*1000,'r.','markersize',0.5); grid on
xlabel('x mm');ylabel('y mm');
subplot(2,1,2),imagesc (Cropped_Image);


 

