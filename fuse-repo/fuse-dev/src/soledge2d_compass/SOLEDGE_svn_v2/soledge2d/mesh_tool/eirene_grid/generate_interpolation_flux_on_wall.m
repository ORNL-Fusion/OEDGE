global eirene
global zones
global Rwall
global Zwall
%chargement des données utiles

R=100*R_;
Z=100*Z_;

ntri_wall=eirene.wall.num;
for k=1:eirene.wall.num
    triangle(k).ntri=eirene.wall.tri(k);
    ntri=triangle(k).ntri;
    triangle(k).k=eirene.wall.tri_k(k);
    triangle(k).i=eirene.wall.tri_i(k);
    triangle(k).j=eirene.wall.tri_j(k);
    A=[triangles_(triangle(k).ntri).neigh1, triangles_(triangle(k).ntri).neigh2, triangles_(triangle(k).ntri).neigh3];
    [a,b]=find(A==0);
    triangle(k).side=b;
end

for n=1:ntri_wall
    k=triangle(n).k;
    i=triangle(n).i;
    j=triangle(n).j;
end

for k=1:zones.num
    Nx=zones.zone(k).Nx;
    Nz=zones.zone(k).Nz;
    zone(k).chi2=zeros(Nx+2,Nz+2);
    zone(k).R2=zeros(Nx+2,Nz+2);
    zone(k).Z2=zeros(Nx+2,Nz+2);
    zone(k).chi2(2:end-1,2:end-1)=zones.zone(k).chi;
    zone(k).R2(2:end-1,2:end-1)=zones.zone(k).gridRc;
    zone(k).Z2(2:end-1,2:end-1)=zones.zone(k).gridZc;
    if(zones.zone(k).Neighbour.north>0) %nord
        zone(k).chi2(Nx+2,2:Nz+1)=zones.zone(zones.zone(k).Neighbour.north).chi(1,:);
        zone(k).R2(Nx+2,2:Nz+1)=zones.zone(zones.zone(k).Neighbour.north).gridRc(1,:);
        zone(k).Z2(Nx+2,2:Nz+1)=zones.zone(zones.zone(k).Neighbour.north).gridZc(1,:);
    else
        zone(k).chi2(Nx+2,2:Nz+1)=zone(k).chi2(Nx+1,2:Nz+1);
        zone(k).R2(Nx+2,2:Nz+1)=zone(k).R2(Nx+1,2:Nz+1);
        zone(k).Z2(Nx+2,2:Nz+1)=zone(k).Z2(Nx+1,2:Nz+1);
    end
    if(zones.zone(k).Neighbour.south>0) %sud
        zone(k).chi2(1,2:Nz+1)=zones.zone(zones.zone(k).Neighbour.south).chi(end,:);
        zone(k).R2(1,2:Nz+1)=zones.zone(zones.zone(k).Neighbour.south).gridRc(end,:);
        zone(k).Z2(1,2:Nz+1)=zones.zone(zones.zone(k).Neighbour.south).gridZc(end,:);
    else
        zone(k).chi2(1,2:Nz+1)=zone(k).chi2(2,2:Nz+1);
        zone(k).R2(1,2:Nz+1)=zone(k).R2(2,2:Nz+1);
        zone(k).Z2(1,2:Nz+1)=zone(k).Z2(2,2:Nz+1);
    end
    if(zones.zone(k).Neighbour.east>0) %est
        zone(k).chi2(2:Nx+1,Nz+2)=zones.zone(zones.zone(k).Neighbour.east).chi(:,1);
        zone(k).R2(2:Nx+1,Nz+2)=zones.zone(zones.zone(k).Neighbour.east).gridRc(:,1);
        zone(k).Z2(2:Nx+1,Nz+2)=zones.zone(zones.zone(k).Neighbour.east).gridZc(:,1);
    else
        zone(k).chi2(2:Nx+1,Nz+2)=zone(k).chi2(2:Nx+1,Nz+1);
        zone(k).R2(2:Nx+1,Nz+2)=zone(k).R2(2:Nx+1,Nz+1);
        zone(k).Z2(2:Nx+1,Nz+2)=zone(k).Z2(2:Nx+1,Nz+1);
    end
    if(zones.zone(k).Neighbour.west>0) %ouest
        zone(k).chi2(2:Nx+1,1)=zones.zone(zones.zone(k).Neighbour.west).chi(:,end);
        zone(k).R2(2:Nx+1,1)=zones.zone(zones.zone(k).Neighbour.west).gridRc(:,end);
        zone(k).Z2(2:Nx+1,1)=zones.zone(zones.zone(k).Neighbour.west).gridZc(:,end);
    else
        zone(k).chi2(2:Nx+1,1)=zone(k).chi2(2:Nx+1,2);
        zone(k).R2(2:Nx+1,1)=zone(k).R2(2:Nx+1,2);
        zone(k).Z2(2:Nx+1,1)=zone(k).Z2(2:Nx+1,2);
    end
end

for k=1:zones.num
    Nx=zones.zone(k).Nx;
    Nz=zones.zone(k).Nz;
    for i=1:Nx
        for j=1:Nz-1
            zone(k).point_est(i,j).i=i;
            zone(k).point_est(i,j).j=j+1;
            zone(k).point_est(i,j).k=k;
        end
        if(zones.zone(k).Neighbour.east>0)
            zone(k).point_est(i,Nz).i=i;
            zone(k).point_est(i,Nz).j=1;
            zone(k).point_est(i,Nz).k=zones.zone(k).Neighbour.east;
        else
            zone(k).point_est(i,Nz).i=0;
            zone(k).point_est(i,Nz).j=0;
            zone(k).point_est(i,Nz).k=0;
        end
    end
    for i=1:Nx
        for j=2:Nz
            zone(k).point_west(i,j).i=i;
            zone(k).point_west(i,j).j=j-1;
            zone(k).point_west(i,j).k=k;
        end
        if(zones.zone(k).Neighbour.west>0)
            zone(k).point_west(i,1).i=i;
            zone(k).point_west(i,1).j=zones.zone(zones.zone(k).Neighbour.west).Nz;
            zone(k).point_west(i,1).k=zones.zone(k).Neighbour.west;
        else
            zone(k).point_west(i,1).i=0;
            zone(k).point_west(i,1).j=0;
            zone(k).point_west(i,1).k=0;
        end
    end
    for j=1:Nz
        for i=1:Nx-1
            zone(k).point_nord(i,j).i=i+1;
            zone(k).point_nord(i,j).j=j;
            zone(k).point_nord(i,j).k=k;
        end
        if(zones.zone(k).Neighbour.north>0)
            zone(k).point_nord(Nx,j).i=1;
            zone(k).point_nord(Nx,j).j=j;
            zone(k).point_nord(Nx,j).k=zones.zone(k).Neighbour.north;
        else
            zone(k).point_nord(Nx,j).i=0;
            zone(k).point_nord(Nx,j).j=0;
            zone(k).point_nord(Nx,j).k=0;
        end
    end
    for j=1:Nz
        for i=2:Nx
            zone(k).point_sud(i,j).i=i-1;
            zone(k).point_sud(i,j).j=j;
            zone(k).point_sud(i,j).k=k;
        end
        if(zones.zone(k).Neighbour.south>0)
            zone(k).point_sud(1,j).i=zones.zone(zones.zone(k).Neighbour.south).Nx;
            zone(k).point_sud(1,j).j=j;
            zone(k).point_sud(1,j).k=zones.zone(k).Neighbour.south;
        else
            zone(k).point_sud(1,j).i=0;
            zone(k).point_sud(1,j).j=0;
            zone(k).point_sud(1,j).k=0;
        end
    end
end