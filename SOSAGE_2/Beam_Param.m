% function Beam_Param take the phase space distribution and return an array
% of normailized emittance (e_x, e_y,e_z) and the beam size

function [e_x,e_y,e_z,sig_x,sig_y] = Beam_Param(x,y,z,px,py,pz)

M = [cov(x,x) cov(x,px); cov(px,x) cov(px,px)]
N = [cov(y,y) cov(y,py); cov(py,y) cov(py,py)]
O = [cov(z,z) cov(z,pz); cov(pz,z) cov(pz,pz)];
c = 299792458;
m_ec2= 0.5109e6;

% normalized emittance

% e_x = (c/m_ec2)*sqrt(M(1,1)*M(2,2) - (M(1,2)^2));
% e_y = (c/m_ec2)*sqrt(N(1,1)*N(2,2) - (N(1,2)^2));

e_x = (c/m_ec2)*sqrt(det(M));
e_y = (c/m_ec2)*sqrt(det(N));
e_z = (c/m_ec2)*sqrt(det(O));

sig_x = sqrt((M(1,1)));
sig_y = sqrt((N(1,1)));


fprintf('e_x(transverse emittance) = %e mm mrad \n',e_x*10^6)
fprintf('e_y(transverse emittance) = %e mm mrad \n',e_y*10^6)
fprintf('e_z(long emittance) = %e Kev mm \n',e_z*0.511*10^6)
fprintf('sig_x (mm) = %e\n',sig_x*1000);
fprintf('sig_y (mm) = %e\n',sig_y*1000);


% unnormalized emittance 

% M_un = cov(x,px/pz);
% N_un = cov(y,py/pz);
% 
% e_x_un = (c/m_ec2)*sqrt(det(M_un));
% e_y_un = (c/m_ec2)*sqrt(det(N_un));
