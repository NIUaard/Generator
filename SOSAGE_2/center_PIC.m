% make picture in the center

function [X,Y] = center_PIC(X,Y,offsetX,offsetY)

a = mean(X);
b = mean(Y);

    X = X - a + offsetX;
    Y = Y - b + offsetY;



