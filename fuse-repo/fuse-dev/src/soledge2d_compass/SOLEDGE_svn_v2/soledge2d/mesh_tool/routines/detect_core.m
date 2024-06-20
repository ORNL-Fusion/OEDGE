global megazone;
global zones;

for k=1:zones.num
     zones.zone(k).MagNeighbour.north=0;
     zones.zone(k).MagNeighbour.south=0;
     zones.zone(k).MagNeighbour.east=0;
     zones.zone(k).MagNeighbour.west=0;
end

for k=1:zones.num
    if(zones.zone(k).Neighbour.north<=0)
        zones.zone(k).Neighbour.north=-2;
        zones.zone(k).MagNeighbour.north=1;
    end
    if(zones.zone(k).Neighbour.south<=0)
        zones.zone(k).Neighbour.south=-2;
        zones.zone(k).MagNeighbour.south=1;
    end
    if(zones.zone(k).Neighbour.east<=0)
        zones.zone(k).Neighbour.east=-5;
        zones.zone(k).MagNeighbour.east=1;
    end
    if(zones.zone(k).Neighbour.west<=0)
        zones.zone(k).Neighbour.west=-5;
        zones.zone(k).MagNeighbour.west=1;
    end
end

for k=1:megazone.num
    nze=megazone.mz(k).list(end);
    nz1=megazone.mz(k).list(1);
    
    if(zones.zone(nze).Neighbour.east==nz1) %periodic
        megazone.mz(k).isperiodic=1;
       for k1=1:length(megazone.mz(k).list)
           nz=megazone.mz(k).list(k1);
           if(zones.zone(nz).Neighbour.south<=0)
               zones.zone(nz).Neighbour.south=-1;
               zones.zone(nz).MagNeighbour.south=1;
           end
       end
    else
        megazone.mz(k).isperiodic=0;
    end
end

