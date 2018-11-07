function [R P] = BiGaussianRP (N, SigmaR, SigmaP, CorrRP, cut)
m_ec2= 0.5109e6;
c = 299792458;

% calculate the thermal momentum

dump1 = sobseqm2(-1,5*N);
V = sobseqm2(4,5*N);

u = V(3,:)'; v = V(4,:)';

R = (sqrt(-2*log(u)).*(sqrt(1-CorrRP^2).*cos(2*pi*v)+CorrRP.*sin(2*pi*v)));
P = sqrt(-2*log(u)).*sin(2*pi*v);

index=find (abs(R)<cut);
R=R(index);
P=P(index);

index=find (abs(P)<cut);
R=R(index);
P=P(index);



R=SigmaR*R(1:N);
P=SigmaP*P(1:N);
