function Z = Gaussian1 (N, Sigma)
m_ec2= 0.5109e6;
c = 299792458;

% calculate the thermal momentum

dump1 = sobseqm2(-1,N);
P_ANGLE_SOB = sobseqm2(4,N);

u = P_ANGLE_SOB(3,:)'; v = P_ANGLE_SOB(4,:)';
Z = Sigma.*sqrt(-2.*log(u)).*cos(2.*pi.*v); % gaussian distribution 

