function [X,Y] = ForceRMS(X,Y,NewRmsX, NewRmsY);
disp ('ForceRMS is on');
    rmsX=std(X);
    rmsY=std(Y);
    X = X.*NewRmsX/rmsX;
    Y = Y.*NewRmsY/rmsY;

