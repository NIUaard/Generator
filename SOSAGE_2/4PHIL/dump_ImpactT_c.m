function Dummy = dump_ImpactT_c(X,Y,Z,PX,PY,PZ,IntegStep);
%
% expect X[m], Y[m]. Z[sec]. PX[eV/c], PY[eV/c], PZ[eV/c]
m_ec2 = 0.5109e6;
cms   = 299792458;

TotalEmissionTimeSec = abs(max(Z)-min(Z))

N	= length(X);
% in impact-T all the particle shoud start with z<0)
aux	= 0;
cdt	= cms*IntegStep;	
% THIS is just to check Ek is the proper value
MeanPz	= mean(PZ);
MeanP   = mean (sqrt(PX.^2+PY.^2+PZ.^2))
MeanEk  = sqrt(MeanP^2*cms^2 + m_ec2^2) -m_ec2

bgx	= PX*cms/m_ec2;
bgy	= PY*cms/m_ec2;
bgz     = PZ*cms/m_ec2;
bg      = sqrt (bgx.^2 + bgy.^2 + bgz.^2);

gamma   = sqrt (bg.^2+1); 
betaz	= PZ*cms./(gamma*m_ec2);
%betaz   = sqrt (1-1./gamma.^2);
Z2	= (Z-max(Z)-1e-16)*cms;
zMean	= mean(Z2)
Phi	= ( mean(Z-max(Z))-1e-16)*cms/((cms)/1.3e9*1.0/360.00);
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
