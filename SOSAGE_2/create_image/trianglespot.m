clear;
maxX=640
maxY=480

meanx=320
meany=240
radius=100

n=1;
Delta=0.9;


for i=1:maxY
    ii=(i-meanx);
    I(i,:)=triangle ([1:maxX]-meanx,i,meany+100, meany-100);
end;
I=I/max(max(I));
subplot (2,2,1)
imagesc (I);
axis xy
subplot (2,2,2)
plot (sum(I,2))
subplot (2,2,3)
plot (sum(I,1))

imwrite(I,'uniformtriangle.png','png')
