clear;
maxX=640
maxY=480

meanx=320
meany=240

n=1;
Delta=0.2;
radius=79.92/2

for i=1:maxX
    I(i,:)=spot (i,[1:maxY],meanx,meany,radius);
end;
imwrite(I,'cathode.png','png')


