function z=spir (x,y,x0,y0,radius,thick)
    r=sqrt(((x-x0).^2+(y-y0).^2));
    th = atan2((x-x0),(y-y0));
    sp = 0*r+radius*sqrt(0.5*(1.00 + th/(2*pi)));
    z=heaviside(thick-abs(r-sp));
