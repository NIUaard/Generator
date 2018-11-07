function Y=slit (X,slitsize,slitspacing,colu, NSlit)
% 
% slitspacing	slit sepration (center-to-center)
% slitsize	size (full width) of the slit
% colu		[depends on  the way X is defined] 
%		1: horizontal slits, 2: vertical slits, 3: longitudinal slits 
% NSlit		number of slits (actually there will be NSlit*2+1 total slits
%       if NSlit<0 then the number of slit is automatically 
%     	chosen to cover the full beam

ind=0;
if (NSlit<0)
   NslitM=abs(round(max(X(:,colu))/slitspacing))
   Nslitm=abs(round(min(X(:,colu))/slitspacing))
   NSlit = max([NslitM, Nslitm]);
end;
for i=1:length(X(:,colu))
  for j=-NSlit:1:NSlit
    if (abs(X(i,colu)-j*slitspacing)<slitsize/2.0)
       ind=ind+1;
       Y(ind,:)=X(i,:);
    end;
  end;
end
