clear;
maxX=640
maxY=480

meanx=320
meany=240
thick=3;
n=1;
Delta=39;  %  separation between beamlets in pixels
radius = 80;  % radius of each beamlets




for i=1:maxX
    I(i,:)=spir (i,[1:maxY],meanx,meany,radius,thick);
end;

imwrite(I,'spirale.png','png')

%save spirale.asc -ASCII I
