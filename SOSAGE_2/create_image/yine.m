clear;
a=imread('vc.png');
a(480,:)=0;
a(:,640)=0;

for i=1:480;
	M(i,:)=spot(i,[1:640], 225.87, 341.71, 48);
end;

imagesc (M.*double(a))

T=M.*double(a);

save yine.dat T -ascii
