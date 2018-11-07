clear;
maxX=640
maxY=480

meanx=320
meany=240
radius=100

Delta=1;
N=1;

for i=1:maxX
    theta=atan2(([1:maxY]-meany),(i-meanx));
    radiu=sqrt(([1:maxY]-meany).^2+(i-meanx)^2);
% UNIFORM TRIANGLE
    I(i,:)=triangle (i-meanx,[1:maxY],meany+100, meany-100,1);
% ramped triangle
%    I(i,:)=triangle (i-meanx,[1:maxY],meany+100, meany-100,1).*...
%                 (1+Delta.*sin(1*theta).*(radiu/radius).^N);
%  standard distribution used in the prstans
%    I(i,:)=spot (i,[1:maxY],meanx,meany,radius).* ...
%                 (1+Delta.*cos(1*theta).*(radiu/radius).^N);
end;
I=I/max(max(I));
imagesc (I);
axis xy
imwrite(I,'tmp.png','png')
%imwrite(I,'dipole_triang_1.png','png')
