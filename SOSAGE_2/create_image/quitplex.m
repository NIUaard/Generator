clear;
maxX=640
maxY=480

meanx=320
meany=240
radius=50

n=1;
Delta=141.4;  %  separation between beamlets in pixels
radius = 50;  % radius of each beamlets




for i=1:maxX
    I(i,:)=spot (i,[1:maxY],meanx,meany,radius) + ...
        spot (i,[1:maxY],meanx+Delta,meany,radius) + ...
		spot (i,[1:maxY],meanx-Delta,meany,radius) + ...
		spot (i,[1:maxY],meanx,meany+Delta,radius) + ...
		spot (i,[1:maxY],meanx,meany-Delta,radius);
end;
imwrite(I,'quitplex.png','png')

