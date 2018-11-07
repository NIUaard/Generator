% this script read a generator output and replace the 3 and 7 column
% with a luser-define distribution 
clear;


sigmat=30;
filename='ell_5k.ini'

x=0:1000;
% parametrization of the CsTe cathode response time
% see /nfs/pamir/pdata/piot/A0/ellispoidal_bunch_proposal

a =      0.2513 
b =   -0.004084  
c =      0.7378  
d =    -0.01241  

y=a*exp(b*x)+c*exp(d*x);


xg=-1000:1:1000;
yg= exp (-xg.^2/(2*sigmat^2));

Y=conv(yg,y);
X=0:length(Y)-1;
[maxxg, indmax]=max(Y);
subplot (2,2,1)
plot (X,Y/max(Y), xg+indmax-1, yg,'--');
title ('Convoluted DF');

fid=fopen (filename,'r');
a=fscanf(fid,'%e %e %e %e %e %e %e %e %d %d',[10 inf]);
fclose (fid);
a=a';
N=length(a)
%N=10000;

X=X-mean(X);

%
zrang = min(X):0.1:max(X);
zhist =interp1 (X, Y, zrang);
zhist=zhist/max(zhist);

subplot (2,2,2)
plot (zrang, zhist,'--');
title ('Interpolated DF');


W=[zrang',zhist'];
sobseqlong(-1,N,W) ;
ztempo=sobseqlong(4,N,W) ;
Z=ztempo(1,:);

[NNN,zzz]=hist(Z,100);
[Mymax, index]=max(NNN);

Z=Z-zzz(index);

subplot (2,2,3)
%Z=Z-mean(Z);
Z(1)=0;
hist(Z,1000);
title ('Monte-Carlo Sampled');

subplot (2,2,4)
Z=Z-mean(Z);
Z(1)=0;
Z=-Z*1e-6;
hist(Z,1000);
title ('Distribution for Astra');
a(:,7)=Z;

filenameout=[filename,'_cste']

fid=fopen (filenameout,'w');
for i=1:N;
  fprintf (fid,'% .4E % .4E % .4E % .4E % .4E % .4E % .4E % .4E % d % d \n', ...
               a(i,1), a(i,2), a(i,3), a(i,4), a(i,5), a(i,6), a(i,7), a(i,8), a(i,9), a(i,10));
end;
fclose(fid);
 
