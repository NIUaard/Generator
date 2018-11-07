z=-10:0.001:10;
pz=temporaldist_mubunching(z,6,0.5,0.3,0.3);

pz=pz/max(pz);


W=[z',pz'];
N=500000;
sobseqlong(-1,N,W) ;
Z=sobseqlong(4,N,W) ;

hist (Z(1,:),1000);
