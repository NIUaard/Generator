function z=ring_fun (x,y,x0,y0,r1,r2)
  A = r1^2-((x-x0).^2+(y-y0).^2);
  B =  r2^2-((x-x0).^2+(y-y0).^2);
 
%   z1 = 1-heaviside(A);
%   z2 = heaviside(B);    
% z = z1  + z2;

z1 =  heaviside(A);
z2 = - heaviside(B);
z = z1+z2;
