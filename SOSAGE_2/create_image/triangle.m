function z=triangle (x,y,ytopline,yabspoint,f)
  for i=1:length(y)
    if (y(i)<ytopline)
       z(i)=1-heaviside(abs(f*x)-y(i)+yabspoint);
%       z(i)=1;
    else
       z(i)=0;
    end;
   end;
   
    
