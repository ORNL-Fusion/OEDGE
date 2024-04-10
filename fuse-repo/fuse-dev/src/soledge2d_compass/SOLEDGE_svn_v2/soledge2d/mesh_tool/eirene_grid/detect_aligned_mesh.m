global zones;

for k=1:zones.num
    Nx=zones.zone(k).Nx;
    Nz=zones.zone(k).Nz;
    zones.zone(k).isaligned=zeros(Nx,Nz);
end