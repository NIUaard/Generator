function [X,Y] = RotateXY(X,Y,theta,N)

for i = 1:N
    
    tmp = X(i)*cos(theta) - Y(i)*sin(theta);
    Y(i) = X(i)*sin(theta) + Y(i)*cos(theta);
    X(i) = tmp;
    
end
