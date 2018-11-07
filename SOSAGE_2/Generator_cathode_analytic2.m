
% gen_distrib will generate a distribution for the 6D phase to be an 
% input  for Astra or for Parmela
% the code used a quasi random code the sobol sequence sobol.m was taken
% from http://www.math.uic.edu/~hanson/mcs507/cp4f04.html.




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
Ekin = 0.55;	% excess kinetic energy in eV (for thermal emittance)
sig_t = 0.1;	% laser pulse length in ps (laser assumed to be  Gaussian)
N = 500;	% number of macroparticles to be generated

%---------------------------------------------------------------------------
%           main program
%
lambdamu=[50 100 150 200 250 300 400 500]
for j=1:length(lambdamu)
% you need to have run the Generator_cathode_analystic(1).m and did not clear 
% any variiable. This just generate other longitudinal distribuion
% generate temploral laser profile
zrang =-20:0.001:20;
modul = 0.1
%lambdamu=500
lambda = lambdamu(j)*1e-6/(299792458)*1e12
zhist =temporaldist_mubunching(zrang,9.6,0.5,lambda,modul);
zhist = zhist/max(zhist);

W=[zrang',zhist'];
sobseqlong(-1,N,W) ;
ztempo=sobseqlong(4,N,W) ;
Z=ztempo(1,:);


Z= Gaussian1 (N,sig_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanPZ = mean(PZ);

q=q/N; % charge per macroparticle

% force the X and Y rms to New values
% ForceRMS(X,Y,1e-3, 1e-3);
[X,Y] = ScaleXY(X,Y,0.5);
% below values for Q = 1nc
Myrmsx = 60e-8;
Myrmsy = 60e-8;

[X,Y]= Rescale (X,Y, Myrmsx,Myrmsy);


%display beam parameters
Beam_Param(X,Y,0.*Z,PX,PY,PZ);

dump_ImpactT(X,Y,Z*1e-12,PX,PY,PZ);
%dump_ImpactT_c(X,Y,Z*1e-12,PX,PY,PZ,1e-12);
%dump_Parmela (X*100, PX/m_ec2, Y*100, PY/m_ec2, X*100, PZ/m_ec2,q);
% plot for image
figure(2),subplot(2,1,1),plot(X(2:N).*1000,Y(2:N).*1000,'r.','markersize',0.5); grid on
xlabel('x mm');ylabel('y mm');
subplot(2,1,2),imagesc (Cropped_Image);
eval(sprintf ('system(''cp partcl.data partcl.data_m%2d_l%d'')',modul*100, lambdamu(j)))
end
 

