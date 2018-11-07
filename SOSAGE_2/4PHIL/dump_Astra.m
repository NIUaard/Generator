function f= dump_Astra(X,Y,Z,PX,PY,PZ,clock,q,stat,flag,N,meanE);

         % dump to Astra 
fid=fopen ('Input_Astra.dat','w');

fprintf(fid,'%e %e %e %e %e %e %e %e %d %d\n',0,0,0,0,0,meanE,0,q,stat,flag);
    for i=1:N
	fprintf(fid,'%e %e %e %e %e %e %e %e %d %d\n',X(i),Y(i),0,PX(i),PY(i),PZ(i)-meanE,clock(i),q,stat,flag);
    end
    % px py and pz are in ev/c

fclose(fid);
f = 'Input_Astra.dat';
