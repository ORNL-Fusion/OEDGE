global Rwall
global Zwall

lw=length(Rwall);
Rs=circshift(Rwall,-1);
Zs=circshift(Zwall,-1);
Rs2=circshift(Rwall,1);
Zs2=circshift(Zwall,1);

shift=1e-5; %m

npt=0;
Rnew=[];
Znew=[];
ind=[];
for k=1:lw
    d=sqrt((Rs2(k)-Rwall(k))^2+(Zs2(k)-Zwall(k))^2);
    if(d>0)
        npt=npt+1;
        Rnew=[Rnew,Rwall(k)];
        Znew=[Znew,Zwall(k)];
    end
    ind=[ind,npt];
end
if(ind(1)==0)
    ind(1)=ind(end);
end
        
lw=length(Rnew);
Rs=circshift(Rnew',-1)';
Zs=circshift(Znew',-1)';
Rs2=circshift(Rnew',1)';
Zs2=circshift(Znew',1)';

for k=1:lw
    
        
   vR=Rs(k)-Rnew(k);
   vZ=Zs(k)-Znew(k);
   norm=sqrt(vR^2+vZ^2);
   vR=vR/norm;
   vZ=vZ/norm;
   vR2=Rs2(k)-Rnew(k);
   vZ2=Zs2(k)-Znew(k);
   norm=sqrt(vR2^2+vZ2^2);
   vR2=vR2/norm;
   vZ2=vZ2/norm;
   normR=vR+vR2;
   normZ=vZ+vZ2;
   norm=sqrt(normR.^2+normZ.^2);
   normR=normR./norm;
   normZ=normZ./norm;
   if(norm==0) 
       %colinear thus 90°
       normR=-vZ2;
       normZ=vR2;
   end
   pR=Rnew(k)+normR*shift;
   pZ=Znew(k)+normZ*shift;
   if(inpolygon(pR,pZ,Rwall,Zwall))
       Rnew2(k)=pR;
       Znew2(k)=pZ;
   else
       pR=Rnew(k)-normR*shift;
       pZ=Znew(k)-normZ*shift;
       Rnew2(k)=pR;
       Znew2(k)=pZ;
   end
end

eirene.Rwall=Rnew2(ind);
eirene.Zwall=Znew2(ind);

% eirene.Rwall=[eirene.Rwall,eirene.Rwall(1)];
% eirene.Zwall=[eirene.Zwall,eirene.Zwall(1)];

% figure
% plot(Rwall,Zwall,'k.-')
% hold on
% plot(eirene.Rwall,eirene.Zwall,'r+-')
% axis equal