global eirene;
global zones;

vector1x=eirene.R(zones.zone(1).knotC(1,1))-eirene.R(zones.zone(1).knotA(1,1));
vector1y=eirene.Z(zones.zone(1).knotC(1,1))-eirene.Z(zones.zone(1).knotA(1,1));

vector2x=eirene.R(zones.zone(1).knotB(1,1))-eirene.R(zones.zone(1).knotA(1,1));
vector2y=eirene.Z(zones.zone(1).knotB(1,1))-eirene.Z(zones.zone(1).knotA(1,1));

vecprod=vector1x*vector2y-vector1y*vector2x;
if(vecprod>=0)
    disp('indirect orientation');
    eirene.direct=0;
else
    disp('direct orientation');
    eirene.direct=1;
end