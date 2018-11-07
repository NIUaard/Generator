
function [Xnew,Ynew]=Rescale (X,Y, Myrmsx, Myrmsy)


Xnew =X.*(Myrmsx./std(X));
Ynew = Y.*(Myrmsy./std(Y));