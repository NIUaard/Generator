function Dummy = dump_ImpactT(X,Y,Z,PX,PY,PZ);
%
% expect X[m], Y[m]. Z[sec]. PX[eV/c], PY[eV/c], PZ[eV/c]
%
% History:
%  11-19-2007: removed the integstep input sincxe not need in V1.5
% 
 disp ('UPDATED TO RUN WITH VERSION 1.5 of IMPACT-T, 03/31/2007')
 
m_ec2 = 0.5109e6;
cms   = 299792458;

TotalEmissionTimeSec = abs(max(Z)-min(Z))
ps_per_deg	=  1/1.3e9/360
TranslatedTime	= abs(mean(Z)-max(Z))
phaseShift 	=  TranslatedTime/(ps_per_deg)
N	= length(X);
% in impact-T all the particle shoud start with z<0)
aux	= 0;
% THIS is just to check Ek is the proper value
MeanPz	= mean(PZ);
MeanP   = mean (sqrt(PX.^2+PY.^2+PZ.^2))
MeanEk  = sqrt(MeanP^2*cms^2 + m_ec2^2) -m_ec2
MeanEkZ = sqrt(mean(PZ)^2*cms^2 + m_ec2^2) -m_ec2

bgx	= PX*cms/m_ec2;
bgy	= PY*cms/m_ec2;
bgz     = PZ*cms/m_ec2;
bg      = sqrt (bgx.^2 + bgy.^2 + bgz.^2);

gamma   = sqrt (bg.^2+1); 
betaz	= PZ*cms./(gamma*m_ec2);
%betaz   = sqrt (1-1./gamma.^2);
betazA  = mean(betaz)
Z2	= (Z-max(Z)-1e-16)*cms.*mean(betaz)+1e-16;
%Z2	= (Z-mean(Z)-1e-16)*cms.*mean(betaz);
zMean	= mean(Z2)
Phi	= (mean(Z-max(Z))-1e-16)*cms/((cms)/1.3e9*1.0/360.00)
SigmZ	= std(Z2)
bgzA    = mean (bgz)
%

fid=fopen ('partcl.data','w');

fprintf(fid, '%d\n',N+1);
		
fprintf(fid, '%13.5e%13.5e%13.5e%13.5e%13.5e%13.5e\n', ...
        	aux,aux,aux,aux,zMean,bgzA);
for i=1:N
	fprintf (fid, '%13.5e%13.5e%13.5e%13.5e%13.5e%13.5e\n',...
    		X(i),  bgx(i), ...
		Y(i),  bgy(i), ...
		Z2(i), bgz(i));
end

fclose(fid);
