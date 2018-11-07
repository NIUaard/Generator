
function [X,Y] = ScaleXY(X,Y,bin_size);


    X = X.*bin_size;
    Y = Y.*bin_size;

