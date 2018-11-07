function [PX, PY, PZ]= Cold (N, Ekin)
m_ec2= 0.5109e6;
c = 299792458;

% calculate the thermal momentum
p = (1/c)*sqrt(2*Ekin*m_ec2) % +Ekin^2);

dump1 = sobseqm2(-1,N);
P_ANGLE_SOB = sobseqm2(4,N);
the =2.*pi.*P_ANGLE_SOB(1,:)'; phi = -pi/2+pi.*P_ANGLE_SOB(2,:)';

% calculate momentum due to thermal emission
    PX = 0.*cos(the).*sin(phi);
    PY = 0.*sin(the).*sin(phi);
    PZ = p+0*cos(phi);

