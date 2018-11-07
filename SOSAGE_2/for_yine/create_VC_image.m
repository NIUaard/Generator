clear;
maxX=640
maxY=480

meanx=320
meany=240
radius=100

n=2;
Delta=0.5;


for i=1:maxX
    I(i,:)=... %spot (i,[1:maxY],meanx,meany,radius).* ...
              1 *   (1+Delta*cos(n*2*pi/radius*(i-meanx))).* ...
		 (1+Delta*cos(n*2*pi/radius*([1:maxY]-meany)));
end;
I=I/max(max(I));
imagesc (I);

imwrite(I,'modulation_2_05.png','png')
