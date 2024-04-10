global zones;
global megazone;
global r2D;
global z2D;
global flux2D;
global Rwall;
global Zwall;
global X_points;
global psicore;
global use_penalization;

mesh_file=uiputfile();

hdf5write(mesh_file,'/NZones',int64(zones.num));
for k=1:zones.num
    %Nx
    hdf5write(mesh_file,['/zone',num2str(k),'/Nx'],int64(zones.zone(k).Nx),'WriteMode','append');
    %Nz
    hdf5write(mesh_file,['/zone',num2str(k),'/Nz'],int64(zones.zone(k).Nz),'WriteMode','append');
    %x
    hdf5write(mesh_file,['/zone',num2str(k),'/x'],zones.zone(k).x,'WriteMode','append');
    %z
    hdf5write(mesh_file,['/zone',num2str(k),'/z'],zones.zone(k).z,'WriteMode','append');
    %xm
    hdf5write(mesh_file,['/zone',num2str(k),'/xm'],zones.zone(k).xb(1:end-1,1:end-1),'WriteMode','append');
    %xp
    hdf5write(mesh_file,['/zone',num2str(k),'/xp'],zones.zone(k).xb(2:end,1:end-1),'WriteMode','append');
    %zm
    hdf5write(mesh_file,['/zone',num2str(k),'/zm'],zones.zone(k).zb(1:end-1,1:end-1),'WriteMode','append');
    %zp
    hdf5write(mesh_file,['/zone',num2str(k),'/zp'],zones.zone(k).zb(1:end-1,2:end),'WriteMode','append');
    %xmin
    hdf5write(mesh_file,['/zone',num2str(k),'/xmin'],zones.zone(k).xb(1,1),'WriteMode','append');
    %xmax
    hdf5write(mesh_file,['/zone',num2str(k),'/xmax'],zones.zone(k).xb(end,1),'WriteMode','append');
    %zmin
    hdf5write(mesh_file,['/zone',num2str(k),'/zmin'],zones.zone(k).zb(1,1),'WriteMode','append');
    %zmax
    hdf5write(mesh_file,['/zone',num2str(k),'/zmax'],zones.zone(k).zb(1,end),'WriteMode','append');
    %chi
    if(use_penalization)
        hdf5write(mesh_file,['/zone',num2str(k),'/chi'],zones.zone(k).chi,'WriteMode','append');
    else
        hdf5write(mesh_file,['/zone',num2str(k),'/chi'],zeros(size(zones.zone(k).chi)),'WriteMode','append');
    end
    
    %Neighbors
    hdf5write(mesh_file,['/zone',num2str(k),'/Neighbors/North'],int64(zones.zone(k).Neighbour.north),'WriteMode','append');
    hdf5write(mesh_file,['/zone',num2str(k),'/Neighbors/South'],int64(zones.zone(k).Neighbour.south),'WriteMode','append');
    if(zones.zone(k).Neighbour.east==-5)
        if(use_penalization)
            hdf5write(mesh_file,['/zone',num2str(k),'/Neighbors/East'],int64(zones.zone(k).Neighbour.east),'WriteMode','append');
        else
            hdf5write(mesh_file,['/zone',num2str(k),'/Neighbors/East'],int64(-3),'WriteMode','append');
        end
    else
        hdf5write(mesh_file,['/zone',num2str(k),'/Neighbors/East'],int64(zones.zone(k).Neighbour.east),'WriteMode','append');
    end
    if(zones.zone(k).Neighbour.west==-5)
        if(use_penalization)
            hdf5write(mesh_file,['/zone',num2str(k),'/Neighbors/West'],int64(zones.zone(k).Neighbour.west),'WriteMode','append');
        else
            hdf5write(mesh_file,['/zone',num2str(k),'/Neighbors/West'],int64(-3),'WriteMode','append');
        end
    else
        hdf5write(mesh_file,['/zone',num2str(k),'/Neighbors/West'],int64(zones.zone(k).Neighbour.west),'WriteMode','append');
    end
    %MagNeighbors
    hdf5write(mesh_file,['/zone',num2str(k),'/MagNeighbors/North'],int64(zones.zone(k).MagNeighbour.north),'WriteMode','append');
    hdf5write(mesh_file,['/zone',num2str(k),'/MagNeighbors/South'],int64(zones.zone(k).MagNeighbour.south),'WriteMode','append');
    hdf5write(mesh_file,['/zone',num2str(k),'/MagNeighbors/East'],int64(zones.zone(k).MagNeighbour.east),'WriteMode','append');
    hdf5write(mesh_file,['/zone',num2str(k),'/MagNeighbors/West'],int64(zones.zone(k).MagNeighbour.west),'WriteMode','append');
    
    Nx=zones.zone(k).Nx;
    Nz=zones.zone(k).Nz;
    %Bphi
    Bp=zeros(Nx+2,Nz+2);
    Bp(2:end-1,2:end-1)=zones.zone(k).Bphi;
    Bp(1,2:end-1)=zones.zone(k).SouthP.Bphi';
    Bp(end,2:end-1)=zones.zone(k).NorthP.Bphi';
    Bp(2:end-1,1)=zones.zone(k).WestP.Bphi;
    Bp(2:end-1,end)=zones.zone(k).EastP.Bphi;
    hdf5write(mesh_file,['/zone',num2str(k),'/Bphi'],Bp,'WriteMode','append');
    %Br
    Bp=zeros(Nx+2,Nz+2);
    Bp(2:end-1,2:end-1)=zones.zone(k).Br;
    Bp(1,2:end-1)=zones.zone(k).SouthP.Br';
    Bp(end,2:end-1)=zones.zone(k).NorthP.Br';
    Bp(2:end-1,1)=zones.zone(k).WestP.Br;
    Bp(2:end-1,end)=zones.zone(k).EastP.Br;
    hdf5write(mesh_file,['/zone',num2str(k),'/Br'],Bp,'WriteMode','append');
    %Bz
    Bp=zeros(Nx+2,Nz+2);
    Bp(2:end-1,2:end-1)=zones.zone(k).Bz;
    Bp(1,2:end-1)=zones.zone(k).SouthP.Bz';
    Bp(end,2:end-1)=zones.zone(k).NorthP.Bz';
    Bp(2:end-1,1)=zones.zone(k).WestP.Bz;
    Bp(2:end-1,end)=zones.zone(k).EastP.Bz;
    hdf5write(mesh_file,['/zone',num2str(k),'/Bz'],Bp,'WriteMode','append');
    
    %Rgeom
    Bp=zeros(Nx+2,Nz+2);
    Bp(2:end-1,2:end-1)=zones.zone(k).gridRc;
    Bp(1,2:end-1)=zones.zone(k).SouthP.R';
    Bp(end,2:end-1)=zones.zone(k).NorthP.R';
    Bp(2:end-1,1)=zones.zone(k).WestP.R;
    Bp(2:end-1,end)=zones.zone(k).EastP.R;
    hdf5write(mesh_file,['/zone',num2str(k),'/Rgeom'],Bp,'WriteMode','append');
    %Zgeom
    Bp=zeros(Nx+2,Nz+2);
    Bp(2:end-1,2:end-1)=zones.zone(k).gridZc;
    Bp(1,2:end-1)=zones.zone(k).SouthP.Z';
    Bp(end,2:end-1)=zones.zone(k).NorthP.Z';
    Bp(2:end-1,1)=zones.zone(k).WestP.Z;
    Bp(2:end-1,end)=zones.zone(k).EastP.Z;
    hdf5write(mesh_file,['/zone',num2str(k),'/Zgeom'],Bp,'WriteMode','append');
    %Rcorner
    hdf5write(mesh_file,['/zone',num2str(k),'/Rcorner'],zones.zone(k).gridR,'WriteMode','append');
    %Zcorner
    hdf5write(mesh_file,['/zone',num2str(k),'/Zcorner'],zones.zone(k).gridZ,'WriteMode','append');
    
end

hdf5write(mesh_file,'/NMegazones',int64(megazone.num),'WriteMode','append');
for k=1:megazone.num
    l=length(megazone.mz(k).list);
    hdf5write(mesh_file,['/megazone',num2str(k),'/size'],int64(l),'WriteMode','append');
    hdf5write(mesh_file,['/megazone',num2str(k),'/configuration'],int64(megazone.mz(k).list),'WriteMode','append');
    hdf5write(mesh_file,['/megazone',num2str(k),'/isperiodic'],int64(megazone.mz(k).isperiodic),'WriteMode','append');
end

hdf5write(mesh_file,'/config/psi',flux2D,'WriteMode','append');
hdf5write(mesh_file,'/config/r',r2D,'WriteMode','append');
hdf5write(mesh_file,'/config/z',z2D,'WriteMode','append');
hdf5write(mesh_file,'/config/nsep',int64(X_points.num),'WriteMode','append');
for k=1:X_points.num
    hdf5write(mesh_file,['/config/psisep',num2str(k)],X_points.psi(k),'WriteMode','append');
end
hdf5write(mesh_file,'/config/psicore',psicore,'WriteMode','append');

hdf5write(mesh_file,'/wall/R',Rwall,'WriteMode','append');
hdf5write(mesh_file,'/wall/Z',Zwall,'WriteMode','append');
