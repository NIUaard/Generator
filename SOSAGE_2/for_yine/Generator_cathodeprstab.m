
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
angle 	= -45*pi/180;
bin_size = 10e-6;
%--------------------------------------------------------------------------

% input from the user

% 
%fname = 'LFDneg10.bmp'
%fname = 'vc.png'
fname1 = 'CAM-XVC-20081216-123710_524.pgm'
fname2 = 'CAM-XVC-20081216-123711_507.pgm'
fname3 = 'CAM-XVC-20081216-123712_553.pgm'
fbkgrd1=  'CAM-XVC-20081216-123747_526.pgm'
fbkgrd2=  'CAM-XVC-20081216-123746_483.pgm'
fbkgrd3=  'CAM-XVC-20081216-123749_542.pgm'

Ekin = 0.55;	% excess kinetic energy in eV (for thermal emittance)
sig_t = 1.915;	% laser pulse length in ps (laser assumed to be  Gaussian)
N = 20000;	% number of macroparticles to be generated

%---------------------------------------------------------------------------
%           main program
%
% read and normalize to unity the VC image:
%I = (double(imread(fname1))+ double(imread(fname2))+double(imread(fname3)))/3.00;
%B = (double(imread(fbkgrd1))+double(imread(fbkgrd2))+double(imread(fbkgrd3)))/3.00;
I = double (imread('LFDneg10.bmp'));
B = 0*I;
A0 = I-B; A=A0;
A = medfilt2( A0, [5 3] );
A = A./max(max(A));
% remove some background
% select ROI
Cropped_Image=Crop_Image (A);
[ind1, ind2]=find (Cropped_Image<0.05);

for i=1:length(ind1)
  Cropped_Image(ind1(i), ind2(i))=0;
end;


% generate sobol sequence for the 6D phase space 
% initialization
dump = sobseqm5(-1,N,Cropped_Image);
tic 
%generation
XY_SOB = sobseqm5(4,N,Cropped_Image); % change name is MakeXY.c
toc
X = XY_SOB(1,:)'; Y = XY_SOB(2,:)';

% cook the XY coordinate 
[X,Y] = Offset(X,Y,offsetX,offsetY,2);
% scale such that  each bin =~ 40 micron
%[X,Y] = ScaleXY(X,Y,bin_size);
%[X,Y]= RotateXY(X,Y,angle,N);

% generate momentum
[PX, PY, PZ]= Thermal (N, Ekin);
%[PX, PY, PZ]= Cold (N, Ekin);

% generate temploral laser profile
Z= Gaussian1 (N,sig_t);

meanPZ = mean(PZ);

q=q/N; % charge per macroparticle

% force the X and Y rms to New values
% ForceRMS(X,Y,1e-3, 1e-3);
[X,Y] = ScaleXY(X,Y,0.9);
% below values for Q = 1nc

Myrmsx = 1e-3;
Myrmsy = 1e-3;
%

[X,Y]= Rescale (X,Y, Myrmsx,Myrmsy);

%R=sqrt(X.^2+Y.^2);
%index=find (R>0.5*1e-3);
%X(index,:)=[];
%Y(index,:)=[];
%Z(index,:)=[];
%PX(index,:)=[];
%PY(index,:)=[];
%PZ(index,:)=[];


%display beam parameters
%Beam_Param(X,Y,0.*Z,PX,PY,PZ);

%dump_ImpactT(X,Y,Z*1e-12,PX,PY,PZ)
%,1e-12);
%dump_ImpactT_c(X,Y,Z*1e-12,PX,PY,PZ,1e-12);
%dump_Parmela (X*100, PX/m_ec2, Y*100, PY/m_ec2, X*100, PZ/m_ec2,q);
% plot for image
N=length(X);
figure(2),subplot(2,2,2),plot(X(2:N)*1e3,Y(2:N)*1e3,'r.','markersize',0.5); 
xlabel('x (mm)','fontsize',22);ylabel('y (mm)','fontsize',22);
subplot(2,2,1),imagesc (Cropped_Image);
axis xy

 

