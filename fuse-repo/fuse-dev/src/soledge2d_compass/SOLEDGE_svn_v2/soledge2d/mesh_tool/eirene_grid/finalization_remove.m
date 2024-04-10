global zones
global eirene

for k=1:zones.num
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            zones.zone(k).knotA(i,j)=eirene.nknots_new(zones.zone(k).knotA(i,j));
            zones.zone(k).knotB(i,j)=eirene.nknots_new(zones.zone(k).knotB(i,j));
            zones.zone(k).knotC(i,j)=eirene.nknots_new(zones.zone(k).knotC(i,j));
            zones.zone(k).knotD(i,j)=eirene.nknots_new(zones.zone(k).knotD(i,j));
        end
    end
end

eirene.R=[Rknots,Rnewknots];
eirene.Z=[Zknots,Znewknots];