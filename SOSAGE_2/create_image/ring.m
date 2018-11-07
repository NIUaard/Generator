clear all;
maxX=640;
maxY=480;

meanx=320;
meany=240;
Delta = 6;

radius1=150.0;
radius2 =140.0;
thick = radius2 - radius1;
for i=1:maxX
     %I(i,:)=spot (i,[1:maxY],meanx,meany,radius1) +...
      I(i,:) =   ring_fun(i,[1:maxY],meanx,meany,radius1,radius2) ;
    
end;
imwrite(I,'ring.png','png')

save ring.asc -ASCII I
