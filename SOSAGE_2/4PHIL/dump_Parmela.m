function f = dump_Parmela(X,PX,Y,PY,Z,PZ,q);


 fid=fopen ('Input_Parmela.dat','w');
% dump to parmela 

    for i=1:N
	fprintf(fid,'%e %e %e %e %e %e %e %e\n',X(i),PX(i),Y(i),PY(i),Z(i),PZ(i),q);
    end

fclose(fid); 

