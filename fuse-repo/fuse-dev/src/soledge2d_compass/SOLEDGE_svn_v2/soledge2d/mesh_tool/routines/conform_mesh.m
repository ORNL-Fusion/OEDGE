global Pmegazone
global zones
global megazone
global X_points

for k=1:Pmegazone.num
   for k1=1:length(Pmegazone.mz(k).list)-1
        nz1=Pmegazone.mz(k).list(k1);
        nz2=Pmegazone.mz(k).list(k1+1);
       zones.zone(nz2).gridR(1,:)=zones.zone(nz1).gridR(end,:);
       zones.zone(nz2).gridZ(1,:)=zones.zone(nz1).gridZ(end,:);
   end
end

for k=1:megazone.num
   for k1=1:length(megazone.mz(k).list)-1
        nz1=megazone.mz(k).list(k1);
        nz2=megazone.mz(k).list(k1+1);
       zones.zone(nz2).gridR(:,1)=zones.zone(nz1).gridR(:,end);
       zones.zone(nz2).gridZ(:,1)=zones.zone(nz1).gridZ(:,end);
   end
end

for k=1:zones.num
   if((zones.zone(k).pA.coord(3)==-1)&&(zones.zone(k).pA.coord(4)==1))
       zones.zone(k).gridR(1,1)=X_points.R(zones.zone(k).pA.coord(2));
       zones.zone(k).gridZ(1,1)=X_points.Z(zones.zone(k).pA.coord(2));
   end
   if((zones.zone(k).pB.coord(3)==-1)&&(zones.zone(k).pB.coord(4)==1))
       zones.zone(k).gridR(end,1)=X_points.R(zones.zone(k).pB.coord(2));
       zones.zone(k).gridZ(end,1)=X_points.Z(zones.zone(k).pB.coord(2));
   end
   if((zones.zone(k).pC.coord(3)==-1)&&(zones.zone(k).pC.coord(4)==1))
       zones.zone(k).gridR(end,end)=X_points.R(zones.zone(k).pC.coord(2));
       zones.zone(k).gridZ(end,end)=X_points.Z(zones.zone(k).pC.coord(2));
   end
   if((zones.zone(k).pD.coord(3)==-1)&&(zones.zone(k).pD.coord(4)==1))
       zones.zone(k).gridR(1,end)=X_points.R(zones.zone(k).pD.coord(2));
       zones.zone(k).gridZ(1,end)=X_points.Z(zones.zone(k).pD.coord(2));
   end
   
end