% make picture in the center

function [X,Y] = center_PIC(X,Y,offsetX,offsetY, FLAG)

if (FLAG==1)
  disp ('Offset:: mean is used')
  a = mean(X);
  b = mean(Y);
end;
if (FLAG==2)
  disp ('Offset:: geometric center is used')
  a = (max(X)+min(X))/2;
  b = (max(Y)+min(Y))/2;
end;

    X = X - a + offsetX;
    Y = Y - b + offsetY;



