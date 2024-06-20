global zones;
global X_points;
global eirene;

nzones=zones.num;

nknots=0;

for k=1:nzones
    m=zones.zone(k).Nx;
    p=zones.zone(k).Nz;
    zones.zone(k).nknot=zeros(m+1,p+1);
    zones.zone(k).Rknot=zeros(m+1,p+1);
    zones.zone(k).Zknot=zeros(m+1,p+1);
    % tables to define if X point in zone corners
    zones.zone(k).V=zeros(1,4);
    % tables to define if normal point in zone corners
    zones.zone(k).V2=zeros(1,4);
end

tol=1.e-5;
for nx=1:X_points.num
    for k=1:zones.num
        d=sqrt((zones.zone(k).gridR(1,1)-X_points.R(nx))^2+(zones.zone(k).gridZ(1,1)-X_points.Z(nx))^2);
        if(d<tol)
            zones.zone(k).V(1)=nx;
        end
        d=sqrt((zones.zone(k).gridR(end,1)-X_points.R(nx))^2+(zones.zone(k).gridZ(end,1)-X_points.Z(nx))^2);
        if(d<tol)
            zones.zone(k).V(2)=nx;
        end
        d=sqrt((zones.zone(k).gridR(end,end)-X_points.R(nx))^2+(zones.zone(k).gridZ(end,end)-X_points.Z(nx))^2);
        if(d<tol)
            zones.zone(k).V(3)=nx;
        end
        d=sqrt((zones.zone(k).gridR(1,end)-X_points.R(nx))^2+(zones.zone(k).gridZ(1,end)-X_points.Z(nx))^2);
        if(d<tol)
            zones.zone(k).V(4)=nx;
        end
    end
end

%count normal points
RN=[];
ZN=[];
for k=1:zones.num
    if((zones.zone(k).Neighbour.north>0)&&(zones.zone(k).Neighbour.east>0))
        RN=[RN,zones.zone(k).gridR(end,end)];
        ZN=[ZN,zones.zone(k).gridZ(end,end)];
    end
end
%remove doublon
RN_=[RN(1)];
ZN_=[ZN(1)];
for k=2:length(RN)
   d=sqrt((RN(k)-RN_).^2+(ZN(k)-ZN_).^2);
   d=min(d);
   if(d>tol)
      RN_=[RN_,RN(k)];
      ZN_=[ZN_,ZN(k)];
   end
end
%remove Xpoints
RN=[];
ZN=[];
for k=1:length(RN_)
    d=1e10;
    for nx=1:X_points.num
        d=min(d,sqrt((RN_(k)-X_points.R(nx))^2+(ZN_(k)-X_points.Z(nx))^2));
    end
    if(d>tol)
        RN=[RN,RN_(k)];
        ZN=[ZN,ZN_(k)];
    end
end
Nnormal=length(RN);

for nn=1:Nnormal
    for k=1:zones.num
        d=sqrt((zones.zone(k).gridR(1,1)-RN(nn))^2+(zones.zone(k).gridZ(1,1)-ZN(nn))^2);
        if(d<tol)
            zones.zone(k).V2(1)=nn;
        end
        d=sqrt((zones.zone(k).gridR(end,1)-RN(nn))^2+(zones.zone(k).gridZ(end,1)-ZN(nn))^2);
        if(d<tol)
            zones.zone(k).V2(2)=nn;
        end
        d=sqrt((zones.zone(k).gridR(end,end)-RN(nn))^2+(zones.zone(k).gridZ(end,end)-ZN(nn))^2);
        if(d<tol)
            zones.zone(k).V2(3)=nn;
        end
        d=sqrt((zones.zone(k).gridR(1,end)-RN(nn))^2+(zones.zone(k).gridZ(1,end)-ZN(nn))^2);
        if(d<tol)
            zones.zone(k).V2(4)=nn;
        end
    end
end

X=[1:X_points.num];
Xnknot=[1:X_points.num];
RX=zeros(1,X_points.num);
ZX=zeros(1,X_points.num);
for k=1:X_points.num
    RX(k)=X_points.R(k);
    ZX(k)=X_points.Z(k);
end

N=[1:Nnormal];
Nnknot=N+X_points.num;

nknots=X_points.num+Nnormal; % 1 xpoints and 8 normal intersection

nX=X_points.num;
nN=Nnormal;

for k=1:nzones
    for i=2:zones.zone(k).Nx
        for j=2:zones.zone(k).Nz
            zones.zone(k).Rknot(i,j)=zones.zone(k).gridR(i,j);
            zones.zone(k).Zknot(i,j)=zones.zone(k).gridZ(i,j);
            nknots=nknots+1;
            zones.zone(k).nknot(i,j)=nknots;
        end
    end
    
    
    %treating north
    if(zones.zone(k).MagNeighbour.north==0)
        north=zones.zone(k).Neighbour.north;
        if(north<k)
            for j=2:zones.zone(k).Nz
                zones.zone(k).Rknot(zones.zone(k).Nx+1,j)=zones.zone(north).Rknot(1,j);
                zones.zone(k).Zknot(zones.zone(k).Nx+1,j)=zones.zone(north).Zknot(1,j);
                zones.zone(k).nknot(zones.zone(k).Nx+1,j)=zones.zone(north).nknot(1,j);
            end
        else
            for j=2:zones.zone(k).Nz
                zones.zone(k).Rknot(zones.zone(k).Nx+1,j)=zones.zone(k).gridR(zones.zone(k).Nx+1,j);
                zones.zone(k).Zknot(zones.zone(k).Nx+1,j)=zones.zone(k).gridZ(zones.zone(k).Nx+1,j);
                nknots=nknots+1;
                zones.zone(k).nknot(zones.zone(k).Nx+1,j)=nknots;
            end
        end
    else
        for j=2:zones.zone(k).Nz
            zones.zone(k).Rknot(zones.zone(k).Nx+1,j)=zones.zone(k).gridR(zones.zone(k).Nx+1,j);
            zones.zone(k).Zknot(zones.zone(k).Nx+1,j)=zones.zone(k).gridZ(zones.zone(k).Nx+1,j);
            nknots=nknots+1;
            zones.zone(k).nknot(zones.zone(k).Nx+1,j)=nknots;
        end
    end
    
    %treating south
    if(zones.zone(k).MagNeighbour.south==0)
        south=zones.zone(k).Neighbour.south;
        if(south<k)
            nz=south;
            for j=2:zones.zone(k).Nz
                zones.zone(k).Rknot(1,j)=zones.zone(nz).Rknot(zones.zone(nz).Nx+1,j);
                zones.zone(k).Zknot(1,j)=zones.zone(nz).Zknot(zones.zone(nz).Nx+1,j);
                zones.zone(k).nknot(1,j)=zones.zone(nz).nknot(zones.zone(nz).Nx+1,j);
            end
        else
            for j=2:zones.zone(k).Nz
                zones.zone(k).Rknot(1,j)=zones.zone(k).gridR(1,j);
                zones.zone(k).Zknot(1,j)=zones.zone(k).gridZ(1,j);
                nknots=nknots+1;
                zones.zone(k).nknot(1,j)=nknots;
            end
        end
    else
        for j=2:zones.zone(k).Nz
            zones.zone(k).Rknot(1,j)=zones.zone(k).gridR(1,j);
            zones.zone(k).Zknot(1,j)=zones.zone(k).gridZ(1,j);
            nknots=nknots+1;
            zones.zone(k).nknot(1,j)=nknots;
        end
    end
    
    %treating east
    if(zones.zone(k).MagNeighbour.east==0)
        east=zones.zone(k).Neighbour.east;
        if(east<k)
            for i=2:zones.zone(k).Nx
                zones.zone(k).Rknot(i,zones.zone(k).Nz+1)=zones.zone(east).Rknot(i,1);
                zones.zone(k).Zknot(i,zones.zone(k).Nz+1)=zones.zone(east).Zknot(i,1);
                zones.zone(k).nknot(i,zones.zone(k).Nz+1)=zones.zone(east).nknot(i,1);
            end
        else
            for i=2:zones.zone(k).Nx
                zones.zone(k).Rknot(i,zones.zone(k).Nz+1)=zones.zone(k).gridR(i,zones.zone(k).Nz+1);
                zones.zone(k).Zknot(i,zones.zone(k).Nz+1)=zones.zone(k).gridZ(i,zones.zone(k).Nz+1);
                nknots=nknots+1;
                zones.zone(k).nknot(i,zones.zone(k).Nz+1)=nknots;
            end
        end
    else
        for i=2:zones.zone(k).Nx
            zones.zone(k).Rknot(i,zones.zone(k).Nz+1)=zones.zone(k).gridR(i,zones.zone(k).Nz+1);
            zones.zone(k).Zknot(i,zones.zone(k).Nz+1)=zones.zone(k).gridZ(i,zones.zone(k).Nz+1);
            nknots=nknots+1;
            zones.zone(k).nknot(i,zones.zone(k).Nz+1)=nknots;
        end
    end
    
    %treating west
    if(zones.zone(k).MagNeighbour.west==0)
        west=zones.zone(k).Neighbour.west;
        if(west<=k)
            nz=west;
            for i=2:zones.zone(k).Nx
                zones.zone(k).Rknot(i,1)=zones.zone(nz).Rknot(i,zones.zone(nz).Nz+1);
                zones.zone(k).Zknot(i,1)=zones.zone(nz).Zknot(i,zones.zone(nz).Nz+1);
                zones.zone(k).nknot(i,1)=zones.zone(nz).nknot(i,zones.zone(nz).Nz+1);
            end
        else
            for i=2:zones.zone(k).Nx
                zones.zone(k).Rknot(i,1)=zones.zone(k).gridR(i,1);
                zones.zone(k).Zknot(i,1)=zones.zone(k).gridZ(i,1);
                nknots=nknots+1;
                zones.zone(k).nknot(i,1)=nknots;
            end
        end
    else
        for i=2:zones.zone(k).Nx
            zones.zone(k).Rknot(i,1)=zones.zone(k).gridR(i,1);
            zones.zone(k).Zknot(i,1)=zones.zone(k).gridZ(i,1);
            nknots=nknots+1;
            zones.zone(k).nknot(i,1)=nknots;
        end
    end
    
    
    %treating corner SW
    if((zones.zone(k).V2(1)==0)&&(zones.zone(k).V(1)==0))
        nex=0;
        if(zones.zone(k).MagNeighbour.west==0) %something west
            west=zones.zone(k).Neighbour.west;
            if(west<k)
                nex=1;
                nz=west;
                zones.zone(k).Rknot(1,1)=zones.zone(nz).Rknot(1,zones.zone(nz).Nz+1);
                zones.zone(k).Zknot(1,1)=zones.zone(nz).Zknot(1,zones.zone(nz).Nz+1);
                zones.zone(k).nknot(1,1)=zones.zone(nz).nknot(1,zones.zone(nz).Nz+1);
            end
        end
        if(zones.zone(k).MagNeighbour.south==0) %something south
            south=zones.zone(k).Neighbour.south;
            if(south<k)
                nex=1;
                nz=south;
                zones.zone(k).Rknot(1,1)=zones.zone(nz).Rknot(zones.zone(nz).Nx+1,1);
                zones.zone(k).Zknot(1,1)=zones.zone(nz).Zknot(zones.zone(nz).Nx+1,1);
                zones.zone(k).nknot(1,1)=zones.zone(nz).nknot(zones.zone(nz).Nx+1,1);
            end
        end
        if(nex==0)
            zones.zone(k).Rknot(1,1)=zones.zone(k).gridR(1,1);
            zones.zone(k).Zknot(1,1)=zones.zone(k).gridZ(1,1);
            nknots=nknots+1;
            zones.zone(k).nknot(1,1)=nknots;
        end
    else
        if(zones.zone(k).V2(1)~=0)
            zones.zone(k).Rknot(1,1)=RN(zones.zone(k).V2(1));
            zones.zone(k).Zknot(1,1)=ZN(zones.zone(k).V2(1));
            zones.zone(k).nknot(1,1)=Nnknot(zones.zone(k).V2(1));
        end
        if(zones.zone(k).V(1)~=0)
            zones.zone(k).Rknot(1,1)=RX(zones.zone(k).V(1));
            zones.zone(k).Zknot(1,1)=ZX(zones.zone(k).V(1));
            zones.zone(k).nknot(1,1)=Xnknot(zones.zone(k).V(1));
        end
    end
    
    %treating corner SE
    if((zones.zone(k).V2(4)==0)&&(zones.zone(k).V(4)==0))
        nex=0;
        if(zones.zone(k).MagNeighbour.east==0) %something east
            east=zones.zone(k).Neighbour.east;
            if(east<k)
                nex=1;
                nz=east;
                zones.zone(k).Rknot(1,zones.zone(k).Nz+1)=zones.zone(nz).Rknot(1,1);
                zones.zone(k).Zknot(1,zones.zone(k).Nz+1)=zones.zone(nz).Zknot(1,1);
                zones.zone(k).nknot(1,zones.zone(k).Nz+1)=zones.zone(nz).nknot(1,1);
            end
        end
        if(zones.zone(k).MagNeighbour.south==0) %something south
            south=zones.zone(k).Neighbour.south;
            if(south<k)
                nex=1;
                nz=south;
                zones.zone(k).Rknot(1,zones.zone(k).Nz+1)=zones.zone(nz).Rknot(zones.zone(nz).Nx+1,zones.zone(k).Nz+1);
                zones.zone(k).Zknot(1,zones.zone(k).Nz+1)=zones.zone(nz).Zknot(zones.zone(nz).Nx+1,zones.zone(k).Nz+1);
                zones.zone(k).nknot(1,zones.zone(k).Nz+1)=zones.zone(nz).nknot(zones.zone(nz).Nx+1,zones.zone(k).Nz+1);
            end
        end
        if(nex==0)
            zones.zone(k).Rknot(1,zones.zone(k).Nz+1)=zones.zone(k).gridR(1,zones.zone(k).Nz+1);
            zones.zone(k).Zknot(1,zones.zone(k).Nz+1)=zones.zone(k).gridZ(1,zones.zone(k).Nz+1);
            nknots=nknots+1;
            zones.zone(k).nknot(1,zones.zone(k).Nz+1)=nknots;
        end
    else
        if(zones.zone(k).V2(4)~=0)
            zones.zone(k).Rknot(1,zones.zone(k).Nz+1)=RN(zones.zone(k).V2(4));
            zones.zone(k).Zknot(1,zones.zone(k).Nz+1)=ZN(zones.zone(k).V2(4));
            zones.zone(k).nknot(1,zones.zone(k).Nz+1)=Nnknot(zones.zone(k).V2(4));
        end
        if(zones.zone(k).V(4)~=0)
            zones.zone(k).Rknot(1,zones.zone(k).Nz+1)=RX(zones.zone(k).V(4));
            zones.zone(k).Zknot(1,zones.zone(k).Nz+1)=ZX(zones.zone(k).V(4));
            zones.zone(k).nknot(1,zones.zone(k).Nz+1)=Xnknot(zones.zone(k).V(4));
        end
    end
    
    %treating corner NE
    if((zones.zone(k).V2(3)==0)&&(zones.zone(k).V(3)==0))
        nex=0;
        if(zones.zone(k).MagNeighbour.east==0) %something east
            east=zones.zone(k).Neighbour.east;
            if(east<=k)
                nex=1;
                nz=east;
                zones.zone(k).Rknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=zones.zone(nz).Rknot(zones.zone(k).Nx+1,1);
                zones.zone(k).Zknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=zones.zone(nz).Zknot(zones.zone(k).Nx+1,1);
                zones.zone(k).nknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=zones.zone(nz).nknot(zones.zone(k).Nx+1,1);
            end
        end
        if(zones.zone(k).MagNeighbour.north==0) %something north
            north=zones.zone(k).Neighbour.north;
            if(north<k)
                nex=1;
                nz=north;
                zones.zone(k).Rknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=zones.zone(nz).Rknot(1,zones.zone(k).Nz+1);
                zones.zone(k).Zknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=zones.zone(nz).Zknot(1,zones.zone(k).Nz+1);
                zones.zone(k).nknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=zones.zone(nz).nknot(1,zones.zone(k).Nz+1);
            end
        end
        if(nex==0)
            zones.zone(k).Rknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=zones.zone(k).gridR(zones.zone(k).Nx+1,zones.zone(k).Nz+1);
            zones.zone(k).Zknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=zones.zone(k).gridZ(zones.zone(k).Nx+1,zones.zone(k).Nz+1);
            nknots=nknots+1;
            zones.zone(k).nknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=nknots;
        end
    else
        if(zones.zone(k).V2(3)~=0)
            zones.zone(k).Rknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=RN(zones.zone(k).V2(3));
            zones.zone(k).Zknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=ZN(zones.zone(k).V2(3));
            zones.zone(k).nknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=Nnknot(zones.zone(k).V2(3));
        end
        if(zones.zone(k).V(3)~=0)
            zones.zone(k).Rknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=RX(zones.zone(k).V(3));
            zones.zone(k).Zknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=ZX(zones.zone(k).V(3));
            zones.zone(k).nknot(zones.zone(k).Nx+1,zones.zone(k).Nz+1)=Xnknot(zones.zone(k).V(3));
        end
    end
    
    %treating corner NW
    if((zones.zone(k).V2(2)==0)&&(zones.zone(k).V(2)==0))
        nex=0;
        if(zones.zone(k).MagNeighbour.west==0) %something west
            west=zones.zone(k).Neighbour.west;
            if(west<k)
                nex=1;
                nz=west;
                zones.zone(k).Rknot(zones.zone(k).Nx+1,1)=zones.zone(nz).Rknot(zones.zone(k).Nx+1,zones.zone(nz).Nz+1);
                zones.zone(k).Zknot(zones.zone(k).Nx+1,1)=zones.zone(nz).Zknot(zones.zone(k).Nx+1,zones.zone(nz).Nz+1);
                zones.zone(k).nknot(zones.zone(k).Nx+1,1)=zones.zone(nz).nknot(zones.zone(k).Nx+1,zones.zone(nz).Nz+1);
            end
        end
        if(zones.zone(k).MagNeighbour.north==0) %something north
            north=zones.zone(k).Neighbour.north;
            if(north<k)
                nex=1;
                nz=north;
                zones.zone(k).Rknot(zones.zone(k).Nx+1,1)=zones.zone(nz).Rknot(1,1);
                zones.zone(k).Zknot(zones.zone(k).Nx+1,1)=zones.zone(nz).Zknot(1,1);
                zones.zone(k).nknot(zones.zone(k).Nx+1,1)=zones.zone(nz).nknot(1,1);
            end
        end
        if(nex==0)
            zones.zone(k).Rknot(zones.zone(k).Nx+1,1)=zones.zone(k).gridR(zones.zone(k).Nx+1,1);
            zones.zone(k).Zknot(zones.zone(k).Nx+1,1)=zones.zone(k).gridZ(zones.zone(k).Nx+1,1);
            nknots=nknots+1;
            zones.zone(k).nknot(zones.zone(k).Nx+1,1)=nknots;
        end
    else
        if(zones.zone(k).V2(2)~=0)
            zones.zone(k).Rknot(zones.zone(k).Nx+1,1)=RN(zones.zone(k).V2(2));
            zones.zone(k).Zknot(zones.zone(k).Nx+1,1)=ZN(zones.zone(k).V2(2));
            zones.zone(k).nknot(zones.zone(k).Nx+1,1)=Nnknot(zones.zone(k).V2(2));
        end
        if(zones.zone(k).V(2)~=0)
            zones.zone(k).Rknot(zones.zone(k).Nx+1,1)=RX(zones.zone(k).V(2));
            zones.zone(k).Zknot(zones.zone(k).Nx+1,1)=ZX(zones.zone(k).V(2));
            zones.zone(k).nknot(zones.zone(k).Nx+1,1)=Xnknot(zones.zone(k).V(2));
        end
    end
    
end


% noeud qui restent (intersection normale de 4 zones)


knots=[];
R=[];
Z=[];
for i=1:nzones
    Nx=zones.zone(i).Nx;
    Nz=zones.zone(i).Nz;
    for j=1:Nx+1
        for k=1:Nz+1
            if(zones.zone(i).nknot(j,k)>0)
                if(length(find(knots==zones.zone(i).nknot(j,k)))==0)
                    knots=[knots,zones.zone(i).nknot(j,k)];
                    R=[R,zones.zone(i).Rknot(j,k)];
                    Z=[Z,zones.zone(i).Zknot(j,k)];
                end
            end
        end
    end
end




% on sauvegarde le numero des noeud eirene dans les coins des mailles
% soledge
for k=1:nzones
    zones.zone(k).knotA=zones.zone(k).nknot(1:zones.zone(k).Nx,1:zones.zone(k).Nz);
    zones.zone(k).knotB=zones.zone(k).nknot(2:zones.zone(k).Nx+1,1:zones.zone(k).Nz);
    zones.zone(k).knotC=zones.zone(k).nknot(2:zones.zone(k).Nx+1,2:zones.zone(k).Nz+1);
    zones.zone(k).knotD=zones.zone(k).nknot(1:zones.zone(k).Nx,2:zones.zone(k).Nz+1);
end

[A,B]=sort(knots);
knots=A;
eirene.R=R(B);
eirene.Z=Z(B);
eirene.nknots=nknots;
