global zones
global Rwall
global Zwall
global r2D
global z2D
global Br2D
global Bz2D
global Bphi2D

Bphi2D=abs(Bphi2D);

for k=1:zones.num
   
    zones.zone(k).gridRc=0.25*(zones.zone(k).gridR(1:end-1,1:end-1)+...
        zones.zone(k).gridR(1:end-1,2:end)+...
        zones.zone(k).gridR(2:end,2:end)+...
        zones.zone(k).gridR(2:end,1:end-1));
    
    zones.zone(k).gridZc=0.25*(zones.zone(k).gridZ(1:end-1,1:end-1)+...
        zones.zone(k).gridZ(1:end-1,2:end)+...
        zones.zone(k).gridZ(2:end,2:end)+...
        zones.zone(k).gridZ(2:end,1:end-1));
    
    zones.zone(k).chi=zeros(size(zones.zone(k).gridZc));
    [m,p]=size(zones.zone(k).chi);
    zones.zone(k).Nx=m;
    zones.zone(k).Nz=p;
    
    for i=1:m
        for j=1:p
            if(~inpolygon(zones.zone(k).gridRc(i,j),zones.zone(k).gridZc(i,j),Rwall,Zwall))
                zones.zone(k).chi(i,j)=1;
            end
        end
    end
    
    zones.zone(k).Br=interp2(r2D,z2D,Br2D,zones.zone(k).gridRc,zones.zone(k).gridZc);
    zones.zone(k).Bz=interp2(r2D,z2D,Bz2D,zones.zone(k).gridRc,zones.zone(k).gridZc);
    zones.zone(k).Bphi=interp2(r2D,z2D,Bphi2D,zones.zone(k).gridRc,zones.zone(k).gridZc);
    
    
end

clean_chis;