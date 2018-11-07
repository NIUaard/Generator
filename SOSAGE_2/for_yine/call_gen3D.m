clear all;
close all
c = 299792458;

alphax = 0;
betax = 18;
nemitx= 10e-6;
alphay = 0;
betay= 18.75;
nemity= 10e-6;
chirp = -8;
sigz = 3.58e-4;
nemitz = 3e-6;
bg = 30;
N=50000;
[X,Y,Z,PX,PY,PZ]=  Gen_3D_Gaussian(alphax,betax,nemitx,alphay,betay,nemity,chirp,sigz,nemitz,bg,N);

PHSP=[X,Y,Z,PX,PY,PZ];

%PHSPS = slits(PHSP, 1e-3, 2e-3, 1, -999)

%plot (PHSP(:,1), PHSP(:,2),'.',PHSPS(:,1), PHSPS(:,2),'.r','markersize',1)
