function [R P] = BiGaussianRP (N, SigmaR, SigmaP, CorrRP)
m_ec2= 0.5109e6;
c = 299792458;

% calculate the thermal momentum

dump1 = sobseqm2(-1,N);
V = sobseqm2(4,N);

u = V(3,:)'; v = V(4,:)';

R =SigmaR*(sqrt(-2*log(u)).*(sqrt(1-CorrRP^2).*cos(2*pi*v)+CorrRP.*sin(2*pi*v)));
P= SigmaP*sqrt(-2*log(u)).*sin(2*pi*v);
