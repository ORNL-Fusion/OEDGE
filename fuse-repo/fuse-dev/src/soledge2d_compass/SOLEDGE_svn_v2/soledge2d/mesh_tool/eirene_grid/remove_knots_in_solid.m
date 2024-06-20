global eirene;
global zones;
% global Rwall;
% global Zwall;

nknots_new=zeros(1,eirene.nknots);
nknot_new=1;
knots=[];
Rknots=[];
Zknots=[];
for n=1:eirene.nknots
    aligned_so_save=0;
    for k=1:zones.num
        [i,j]=find(zones.zone(k).knotA==n);
        if(zones.zone(k).isaligned(i,j)==1)
            if(inpolygon(zones.zone(k).gridRc(i,j),zones.zone(k).gridZc(i,j),eirene.Rwall,eirene.Zwall))
                aligned_so_save=1;
            end
        end
        [i,j]=find(zones.zone(k).knotB==n);
        if(zones.zone(k).isaligned(i,j)==1)
            if(inpolygon(zones.zone(k).R(i,j),zones.zone(k).Z(i,j),eirene.Rwall,eirene.Zwall))
                aligned_so_save=1;
            end
        end
        [i,j]=find(zones.zone(k).knotC==n);
        if(zones.zone(k).isaligned(i,j)==1)
            if(inpolygon(zones.zone(k).R(i,j),zones.zone(k).Z(i,j),eirene.Rwall,eirene.Zwall))
                aligned_so_save=1;
            end
        end
        [i,j]=find(zones.zone(k).knotD==n);
        if(zones.zone(k).isaligned(i,j)==1)
            if(inpolygon(zones.zone(k).R(i,j),zones.zone(k).Z(i,j),eirene.Rwall,eirene.Zwall))
                aligned_so_save=1;
            end
        end
    end            
    if((inpolygon(eirene.R(n),eirene.Z(n),eirene.Rwall,eirene.Zwall)==1)||(aligned_so_save==1))
        nknots_new(n)=nknot_new;
        knots=[knots,nknot_new];
        Rknots=[Rknots,eirene.R(n)];
        Zknots=[Zknots,eirene.Z(n)];
        nknot_new=nknot_new+1;
    end
end

eirene.nknots=nknot_new-1;
eirene.Rknots=Rknots;
eirene.Zknots=Zknots;

for k=1:zones.num
    zones.zone(k).inplasma=ones(zones.zone(k).Nx,zones.zone(k).Nz);
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(inpolygon(eirene.R(zones.zone(k).knotA(i,j)),eirene.Z(zones.zone(k).knotA(i,j)),eirene.Rwall,eirene.Zwall)==0)
                if(inpolygon(eirene.R(zones.zone(k).knotB(i,j)),eirene.Z(zones.zone(k).knotB(i,j)),eirene.Rwall,eirene.Zwall)==0)
                    if(inpolygon(eirene.R(zones.zone(k).knotC(i,j)),eirene.Z(zones.zone(k).knotC(i,j)),eirene.Rwall,eirene.Zwall)==0)
                        if(inpolygon(eirene.R(zones.zone(k).knotD(i,j)),eirene.Z(zones.zone(k).knotD(i,j)),eirene.Rwall,eirene.Zwall)==0)
                            zones.zone(k).inplasma(i,j)=0;
                        end
                    end
                end
            end
        end
    end
end

eirene.nknots_new=nknots_new;