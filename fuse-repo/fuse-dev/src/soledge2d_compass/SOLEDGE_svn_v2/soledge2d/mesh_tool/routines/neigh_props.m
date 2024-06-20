global zones
global r2D
global z2D
global Br2D
global Bz2D
global Bphi2D

for k=1:zones.num
    zones.zone(k).NorthP.R=(zones.zone(k).gridR(end,2:end)+zones.zone(k).gridR(end,1:end-1))*0.5;
    zones.zone(k).NorthP.Z=(zones.zone(k).gridZ(end,2:end)+zones.zone(k).gridZ(end,1:end-1))*0.5;
    zones.zone(k).NorthP.Br=interp2(r2D,z2D,Br2D,zones.zone(k).NorthP.R,zones.zone(k).NorthP.Z);
    zones.zone(k).NorthP.Bz=interp2(r2D,z2D,Bz2D,zones.zone(k).NorthP.R,zones.zone(k).NorthP.Z);
    zones.zone(k).NorthP.Bphi=interp2(r2D,z2D,Bphi2D,zones.zone(k).NorthP.R,zones.zone(k).NorthP.Z);
    zones.zone(k).NorthP.B=sqrt(zones.zone(k).NorthP.Br.^2+zones.zone(k).NorthP.Bz.^2+zones.zone(k).NorthP.Bphi.^2);
    
    zones.zone(k).SouthP.R=(zones.zone(k).gridR(1,2:end)+zones.zone(k).gridR(1,1:end-1))*0.5;
    zones.zone(k).SouthP.Z=(zones.zone(k).gridZ(1,2:end)+zones.zone(k).gridZ(1,1:end-1))*0.5;
    zones.zone(k).SouthP.Br=interp2(r2D,z2D,Br2D,zones.zone(k).SouthP.R,zones.zone(k).SouthP.Z);
    zones.zone(k).SouthP.Bz=interp2(r2D,z2D,Bz2D,zones.zone(k).SouthP.R,zones.zone(k).SouthP.Z);
    zones.zone(k).SouthP.Bphi=interp2(r2D,z2D,Bphi2D,zones.zone(k).SouthP.R,zones.zone(k).SouthP.Z);
    zones.zone(k).SouthP.B=sqrt(zones.zone(k).SouthP.Br.^2+zones.zone(k).SouthP.Bz.^2+zones.zone(k).SouthP.Bphi.^2);
    
    zones.zone(k).EastP.R=(zones.zone(k).gridR(2:end,end)+zones.zone(k).gridR(1:end-1,end))*0.5;
    zones.zone(k).EastP.Z=(zones.zone(k).gridZ(2:end,end)+zones.zone(k).gridZ(1:end-1,end))*0.5;
    zones.zone(k).EastP.Br=interp2(r2D,z2D,Br2D,zones.zone(k).EastP.R,zones.zone(k).EastP.Z);
    zones.zone(k).EastP.Bz=interp2(r2D,z2D,Bz2D,zones.zone(k).EastP.R,zones.zone(k).EastP.Z);
    zones.zone(k).EastP.Bphi=interp2(r2D,z2D,Bphi2D,zones.zone(k).EastP.R,zones.zone(k).EastP.Z);
    zones.zone(k).EastP.B=sqrt(zones.zone(k).EastP.Br.^2+zones.zone(k).EastP.Bz.^2+zones.zone(k).EastP.Bphi.^2);
    
    zones.zone(k).WestP.R=(zones.zone(k).gridR(2:end,1)+zones.zone(k).gridR(1:end-1,1))*0.5;
    zones.zone(k).WestP.Z=(zones.zone(k).gridZ(2:end,1)+zones.zone(k).gridZ(1:end-1,1))*0.5;
    zones.zone(k).WestP.Br=interp2(r2D,z2D,Br2D,zones.zone(k).WestP.R,zones.zone(k).WestP.Z);
    zones.zone(k).WestP.Bz=interp2(r2D,z2D,Bz2D,zones.zone(k).WestP.R,zones.zone(k).WestP.Z);
    zones.zone(k).WestP.Bphi=interp2(r2D,z2D,Bphi2D,zones.zone(k).WestP.R,zones.zone(k).WestP.Z);
    zones.zone(k).WestP.B=sqrt(zones.zone(k).WestP.Br.^2+zones.zone(k).WestP.Bz.^2+zones.zone(k).WestP.Bphi.^2);
    
end