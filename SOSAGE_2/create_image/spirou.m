r0=1;
th=0:0.002:2*pi;

sp = r0*sqrt(0.5*(1+ th/(2*pi)));
 
 
 plot (sp.*cos(th),sp.*sin(th));
 
