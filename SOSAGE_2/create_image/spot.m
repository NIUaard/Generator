function z=spot (x,y,x0,y0,r)
    z=heaviside(r^2-((x-x0).^2+(y-y0).^2));
    
