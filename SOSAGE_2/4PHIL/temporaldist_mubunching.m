function y=temporaldistrib(z,dt,rt,lambda, m);





y=1+z*0;



for i=1:length(z)

  if (z(i)<-dt+8*rt)

      y(i)=(1./(1+exp(-(z(i)+dt)/rt)));

  end;

  if (z(i)>dt-8*rt)

      y(i)=(1./(1+exp((z(i)-dt)/rt)));

  end;

  if (abs(z(i))<(dt-8*rt))

      y(i)=y(i)+m*cos(2*pi/lambda*z(i));

  end;

end;
