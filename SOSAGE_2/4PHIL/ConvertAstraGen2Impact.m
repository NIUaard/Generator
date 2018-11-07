% convert an output from Astra generator into an input for 
% impact distribution

fname = 'radial2k_short.ini'
fname = 'deck.0300.001'

fid=fopen (fname,'r');

a=fscanf (fid, '%e ', [10 inf]);

cms   = 299792458;

phsp=a(:,2:length(a));
aver=a(:,1);

X = phsp(1,:) + aver(1);
Y = phsp(2,:) + aver(2);
% is generator output
% Z = 1.e-9*(phsp(7,:) + aver(7));
% if astra output
Z = 1/cms*(phsp(3,:) + aver(3));
PX= 1/cms*(phsp(4,:) + aver(4));
PY= 1/cms*(phsp(5,:) + aver(5));
PZ= 1/cms*(phsp(6,:) + aver(6));

Z=Z-min(Z)+1e-10

dump_ImpactT(X,Y,Z,PX,PY,PZ);
