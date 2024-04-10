function reverse_megazone(nmz)

global megazone;
global zones;

for k=1:length(megazone.mz(nmz).list)
    %return zone
    num=megazone.mz(nmz).list(k);
    temp.R=zones.zone(num).east.R;
    temp.Z=zones.zone(num).east.Z;
    zones.zone(num).east.R=zones.zone(num).west.R;
    zones.zone(num).east.Z=zones.zone(num).west.Z;
    zones.zone(num).west.R=temp.R;
    zones.zone(num).west.Z=temp.Z;
    zones.zone(num).north.R=zones.zone(num).north.R(end:-1:1);
    zones.zone(num).north.Z=zones.zone(num).north.Z(end:-1:1);
    zones.zone(num).south.R=zones.zone(num).south.R(end:-1:1);
    zones.zone(num).south.Z=zones.zone(num).south.Z(end:-1:1);
    pA=zones.zone(num).pA;
    pB=zones.zone(num).pB;
    zones.zone(num).pA=zones.zone(num).pD;
    zones.zone(num).pB=zones.zone(num).pC;
    zones.zone(num).pC=pB;
    zones.zone(num).pD=pA;
end
%return megazone
megazone.mz(nmz).list=megazone.mz(nmz).list(end:-1:1);


end