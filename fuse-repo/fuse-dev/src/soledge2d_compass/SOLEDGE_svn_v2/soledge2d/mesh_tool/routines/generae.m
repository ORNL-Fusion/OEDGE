nzones=18;

nknots=0;
AA=load('Neighbors.txt');
BA=load('MagNeighbors.txt');

fid=fopen('Neighbors.txt');
fid2=fopen('MagNeighbors.txt');
for k=1:nzones
    cd Magnetic_input
    filename=strcat('R_',num2str(k,'%.3d'));
    filename=strcat(filename,'.txt');
    R=load(filename);
    [m,p]=size(R);
    filename=strcat('Z_',num2str(k,'%.3d'));
    filename=strcat(filename,'.txt');
    Z=load(filename);
    filename=strcat('Rr_',num2str(k,'%.3d'));
    filename=strcat(filename,'.txt');
    Rr=load(filename);
    filename=strcat('Zr_',num2str(k,'%.3d'));
    filename=strcat(filename,'.txt');
    Zr=load(filename);
    cd ../Mesh
    filename=strcat('chi_',num2str(k,'%.3d'));
    filename=strcat(filename,'.txt');
    chi=load(filename);
    cd ..
    zone(k).R=R;
    zone(k).Z=Z;
    zone(k).Rr=Rr;
    zone(k).Zr=Zr;
    zone(k).chi=chi;
    zone(k).Nx=m;
    zone(k).Nz=p;
    A=fscanf(fid,'%d',4);
    B=fscanf(fid2,'%d',4);
    zone(k).vois(1)=A(1);
    zone(k).vois(2)=A(2);
    zone(k).vois(3)=A(3);
    zone(k).vois(4)=A(4);
    zone(k).isvois(1)=B(1);
    zone(k).isvois(2)=B(2);
    zone(k).isvois(3)=B(3);
    zone(k).isvois(4)=B(4);
    zone(k).nknot=zeros(m+1,p+1);
    zone(k).Rknot=zeros(m+1,p+1);
    zone(k).Zknot=zeros(m+1,p+1);
    % tables to define if X point in zone corners
    zone(k).V=zeros(1,4);
    % tables to define if normal point in zone corners
    zone(k).V2=zeros(1,4);
end
fclose(fid);
fclose(fid2);

% 1 X points in zone corners
zone(2).V(4)=1;
zone(3).V(2)=1;
zone(4).V(1)=1;
zone(6).V(1)=1;
zone(7).V(3)=1;
zone(8).V(4)=1;
zone(17).V(3)=1;
zone(18).V(2)=1;


% 8 normal zone intersections in zone corners
zone(1).V2(3)=8;
zone(2).V2(3)=1;
zone(3).V2(3)=3;
zone(4).V2(2)=1;
zone(4).V2(3)=4;
zone(4).V2(4)=3;
zone(5).V2(2)=8;
zone(6).V2(2)=6;
zone(7).V2(2)=3;
zone(8).V2(2)=4;
zone(8).V2(3)=6;
zone(8).V2(1)=3;
zone(9).V2(4)=1;
zone(9).V2(3)=2;
zone(10).V2(4)=2;
zone(11).V2(1)=1;
zone(11).V2(2)=2;
zone(11).V2(4)=4;
zone(11).V2(3)=5;
zone(12).V2(1)=2;
zone(12).V2(4)=5;
zone(13).V2(2)=5;
zone(13).V2(1)=4;
zone(13).V2(4)=6;
zone(13).V2(3)=7;
zone(14).V2(1)=5;
zone(14).V2(4)=7;
zone(15).V2(1)=6;
zone(15).V2(2)=7;
zone(16).V2(1)=7;
zone(17).V2(4)=8;
zone(18).V2(1)=8;

X=[1];
Xnknot=[1];
RX=zeros(1,1);
ZX=zeros(1,1);

N=[1 2 3 4 5 6 7 8];
Nnknot=[2 3 4 5 6 7 8 9];
RN=zeros(1,8);
ZN=zeros(1,8);

nknots=9; % 1 xpoints and 8 normal intersection

nX=1;
nN=8;

figure(1)
hold on

for n=1:nX
    for k=1:nzones
        if(zone(k).V(1)==n)
            break
        end
    end
    RX(n)=zone(k).Rr(1,1);
    ZX(n)=zone(k).Zr(1,1);
    plot(RX(n),ZX(n),'ko','MarkerSize',8,'LineWidth',4)
end

for n=1:nN
    for k=1:nzones
        if(zone(k).V2(1)==n)
            break
        end
    end
    RN(n)=zone(k).Rr(1,1);
    ZN(n)=zone(k).Zr(1,1);
    plot(RN(n),ZN(n),'co','MarkerSize',8,'LineWidth',4)
    
end



for k=1:nzones
    for i=2:zone(k).Nx
        for j=2:zone(k).Nz
            zone(k).Rknot(i,j)=zone(k).Rr(i,j);
            zone(k).Zknot(i,j)=zone(k).Zr(i,j);
            nknots=nknots+1;
            zone(k).nknot(i,j)=nknots;
        end
    end
    
    
    %treating north
    if(BA(k,1)==0)
        if(AA(k,1)<k)
            for j=2:zone(k).Nz
                zone(k).Rknot(zone(k).Nx+1,j)=zone(AA(k,1)).Rknot(1,j);
                zone(k).Zknot(zone(k).Nx+1,j)=zone(AA(k,1)).Zknot(1,j);
                zone(k).nknot(zone(k).Nx+1,j)=zone(AA(k,1)).nknot(1,j);
            end
        else
            for j=2:zone(k).Nz
                zone(k).Rknot(zone(k).Nx+1,j)=zone(k).Rr(zone(k).Nx+1,j);
                zone(k).Zknot(zone(k).Nx+1,j)=zone(k).Zr(zone(k).Nx+1,j);
                nknots=nknots+1;
                zone(k).nknot(zone(k).Nx+1,j)=nknots;
            end
        end
    else
        for j=2:zone(k).Nz
            zone(k).Rknot(zone(k).Nx+1,j)=zone(k).Rr(zone(k).Nx+1,j);
            zone(k).Zknot(zone(k).Nx+1,j)=zone(k).Zr(zone(k).Nx+1,j);
            nknots=nknots+1;
            zone(k).nknot(zone(k).Nx+1,j)=nknots;
        end
    end
    
    %treating south
    if(BA(k,2)==0)
        if(AA(k,2)<k)
            nz=AA(k,2);
            for j=2:zone(k).Nz
                zone(k).Rknot(1,j)=zone(nz).Rknot(zone(nz).Nx+1,j);
                zone(k).Zknot(1,j)=zone(nz).Zknot(zone(nz).Nx+1,j);
                zone(k).nknot(1,j)=zone(nz).nknot(zone(nz).Nx+1,j);
            end
        else
            for j=2:zone(k).Nz
                zone(k).Rknot(1,j)=zone(k).Rr(1,j);
                zone(k).Zknot(1,j)=zone(k).Zr(1,j);
                nknots=nknots+1;
                zone(k).nknot(1,j)=nknots;
            end
        end
    else
        for j=2:zone(k).Nz
            zone(k).Rknot(1,j)=zone(k).Rr(1,j);
            zone(k).Zknot(1,j)=zone(k).Zr(1,j);
            nknots=nknots+1;
            zone(k).nknot(1,j)=nknots;
        end
    end
    
    %treating east
    if(BA(k,3)==0)
        if(AA(k,3)<k)
            for i=2:zone(k).Nx
                zone(k).Rknot(i,zone(k).Nz+1)=zone(AA(k,3)).Rknot(i,1);
                zone(k).Zknot(i,zone(k).Nz+1)=zone(AA(k,3)).Zknot(i,1);
                zone(k).nknot(i,zone(k).Nz+1)=zone(AA(k,3)).nknot(i,1);
            end
        else
            for i=2:zone(k).Nx
                zone(k).Rknot(i,zone(k).Nz+1)=zone(k).Rr(i,zone(k).Nz+1);
                zone(k).Zknot(i,zone(k).Nz+1)=zone(k).Zr(i,zone(k).Nz+1);
                nknots=nknots+1;
                zone(k).nknot(i,zone(k).Nz+1)=nknots;
            end
        end
    else
        for i=2:zone(k).Nx
            zone(k).Rknot(i,zone(k).Nz+1)=zone(k).Rr(i,zone(k).Nz+1);
            zone(k).Zknot(i,zone(k).Nz+1)=zone(k).Zr(i,zone(k).Nz+1);
            nknots=nknots+1;
            zone(k).nknot(i,zone(k).Nz+1)=nknots;
        end
    end
    
    %treating west
    if(BA(k,4)==0)
        if(AA(k,4)<k)
            nz=AA(k,4);
            for i=2:zone(k).Nx
                zone(k).Rknot(i,1)=zone(nz).Rknot(i,zone(nz).Nz+1);
                zone(k).Zknot(i,1)=zone(nz).Zknot(i,zone(nz).Nz+1);
                zone(k).nknot(i,1)=zone(nz).nknot(i,zone(nz).Nz+1);
            end
        else
            for i=2:zone(k).Nx
                zone(k).Rknot(i,1)=zone(k).Rr(i,1);
                zone(k).Zknot(i,1)=zone(k).Zr(i,1);
                nknots=nknots+1;
                zone(k).nknot(i,1)=nknots;
            end
        end
    else
        for i=2:zone(k).Nx
            zone(k).Rknot(i,1)=zone(k).Rr(i,1);
            zone(k).Zknot(i,1)=zone(k).Zr(i,1);
            nknots=nknots+1;
            zone(k).nknot(i,1)=nknots;
        end
    end
    
    
    %treating corner SW
    if((zone(k).V2(1)==0)&&(zone(k).V(1)==0))
        nex=0;
        if(BA(k,4)==0) %something west
            if(AA(k,4)<k)
                nex=1;
                nz=AA(k,4);
                zone(k).Rknot(1,1)=zone(nz).Rknot(1,zone(nz).Nz+1);
                zone(k).Zknot(1,1)=zone(nz).Zknot(1,zone(nz).Nz+1);
                zone(k).nknot(1,1)=zone(nz).nknot(1,zone(nz).Nz+1);
            end
        end
        if(BA(k,2)==0) %something south
            if(AA(k,2)<k)
                nex=1;
                nz=AA(k,2);
                zone(k).Rknot(1,1)=zone(nz).Rknot(zone(nz).Nx+1,1);
                zone(k).Zknot(1,1)=zone(nz).Zknot(zone(nz).Nx+1,1);
                zone(k).nknot(1,1)=zone(nz).nknot(zone(nz).Nx+1,1);
            end
        end
        if(nex==0)
            zone(k).Rknot(1,1)=zone(k).Rr(1,1);
            zone(k).Zknot(1,1)=zone(k).Zr(1,1);
            nknots=nknots+1;
            zone(k).nknot(1,1)=nknots;
        end
    else
        if(zone(k).V2(1)~=0)
            zone(k).Rknot(1,1)=RN(zone(k).V2(1));
            zone(k).Zknot(1,1)=ZN(zone(k).V2(1));
            zone(k).nknot(1,1)=Nnknot(zone(k).V2(1));
        end
        if(zone(k).V(1)~=0)
            zone(k).Rknot(1,1)=RX(zone(k).V(1));
            zone(k).Zknot(1,1)=ZX(zone(k).V(1));
            zone(k).nknot(1,1)=Xnknot(zone(k).V(1));
        end
    end
    
    %treating corner SE
    if((zone(k).V2(4)==0)&&(zone(k).V(4)==0))
        nex=0;
        if(BA(k,3)==0) %something east
            if(AA(k,3)<k)
                nex=1;
                nz=AA(k,3);
                zone(k).Rknot(1,zone(k).Nz+1)=zone(nz).Rknot(1,1);
                zone(k).Zknot(1,zone(k).Nz+1)=zone(nz).Zknot(1,1);
                zone(k).nknot(1,zone(k).Nz+1)=zone(nz).nknot(1,1);
            end
        end
        if(BA(k,2)==0) %something south
            if(AA(k,2)<k)
                nex=1;
                nz=AA(k,2);
                zone(k).Rknot(1,zone(k).Nz+1)=zone(nz).Rknot(zone(nz).Nx+1,zone(k).Nz+1);
                zone(k).Zknot(1,zone(k).Nz+1)=zone(nz).Zknot(zone(nz).Nx+1,zone(k).Nz+1);
                zone(k).nknot(1,zone(k).Nz+1)=zone(nz).nknot(zone(nz).Nx+1,zone(k).Nz+1);
            end
        end
        if(nex==0)
            zone(k).Rknot(1,zone(k).Nz+1)=zone(k).Rr(1,zone(k).Nz+1);
            zone(k).Zknot(1,zone(k).Nz+1)=zone(k).Zr(1,zone(k).Nz+1);
            nknots=nknots+1;
            zone(k).nknot(1,zone(k).Nz+1)=nknots;
        end
    else
        if(zone(k).V2(4)~=0)
            zone(k).Rknot(1,zone(k).Nz+1)=RN(zone(k).V2(4));
            zone(k).Zknot(1,zone(k).Nz+1)=ZN(zone(k).V2(4));
            zone(k).nknot(1,zone(k).Nz+1)=Nnknot(zone(k).V2(4));
        end
        if(zone(k).V(4)~=0)
            zone(k).Rknot(1,zone(k).Nz+1)=RX(zone(k).V(4));
            zone(k).Zknot(1,zone(k).Nz+1)=ZX(zone(k).V(4));
            zone(k).nknot(1,zone(k).Nz+1)=Xnknot(zone(k).V(4));
        end
    end
    
    %treating corner NE
    if((zone(k).V2(3)==0)&&(zone(k).V(3)==0))
        nex=0;
        if(BA(k,3)==0) %something east
            if(AA(k,3)<k)
                nex=1;
                nz=AA(k,3);
                zone(k).Rknot(zone(k).Nx+1,zone(k).Nz+1)=zone(nz).Rknot(zone(k).Nx+1,1);
                zone(k).Zknot(zone(k).Nx+1,zone(k).Nz+1)=zone(nz).Zknot(zone(k).Nx+1,1);
                zone(k).nknot(zone(k).Nx+1,zone(k).Nz+1)=zone(nz).nknot(zone(k).Nx+1,1);
            end
        end
        if(BA(k,1)==0) %something north
            if(AA(k,1)<k)
                nex=1;
                nz=AA(k,1);
                zone(k).Rknot(zone(k).Nx+1,zone(k).Nz+1)=zone(nz).Rknot(1,zone(k).Nz+1);
                zone(k).Zknot(zone(k).Nx+1,zone(k).Nz+1)=zone(nz).Zknot(1,zone(k).Nz+1);
                zone(k).nknot(zone(k).Nx+1,zone(k).Nz+1)=zone(nz).nknot(1,zone(k).Nz+1);
            end
        end
        if(nex==0)
            zone(k).Rknot(zone(k).Nx+1,zone(k).Nz+1)=zone(k).Rr(zone(k).Nx+1,zone(k).Nz+1);
            zone(k).Zknot(zone(k).Nx+1,zone(k).Nz+1)=zone(k).Zr(zone(k).Nx+1,zone(k).Nz+1);
            nknots=nknots+1;
            zone(k).nknot(zone(k).Nx+1,zone(k).Nz+1)=nknots;
        end
    else
        if(zone(k).V2(3)~=0)
            zone(k).Rknot(zone(k).Nx+1,zone(k).Nz+1)=RN(zone(k).V2(3));
            zone(k).Zknot(zone(k).Nx+1,zone(k).Nz+1)=ZN(zone(k).V2(3));
            zone(k).nknot(zone(k).Nx+1,zone(k).Nz+1)=Nnknot(zone(k).V2(3));
        end
        if(zone(k).V(3)~=0)
            zone(k).Rknot(zone(k).Nx+1,zone(k).Nz+1)=RX(zone(k).V(3));
            zone(k).Zknot(zone(k).Nx+1,zone(k).Nz+1)=ZX(zone(k).V(3));
            zone(k).nknot(zone(k).Nx+1,zone(k).Nz+1)=Xnknot(zone(k).V(3));
        end
    end
    
    %treating corner NW
    if((zone(k).V2(2)==0)&&(zone(k).V(2)==0))
        nex=0;
        if(BA(k,4)==0) %something west
            if(AA(k,4)<k)
                nex=1;
                nz=AA(k,4);
                zone(k).Rknot(zone(k).Nx+1,1)=zone(nz).Rknot(zone(k).Nx+1,zone(nz).Nz+1);
                zone(k).Zknot(zone(k).Nx+1,1)=zone(nz).Zknot(zone(k).Nx+1,zone(nz).Nz+1);
                zone(k).nknot(zone(k).Nx+1,1)=zone(nz).nknot(zone(k).Nx+1,zone(nz).Nz+1);
            end
        end
        if(BA(k,1)==0) %something north
            if(AA(k,1)<k)
                nex=1;
                nz=AA(k,1);
                zone(k).Rknot(zone(k).Nx+1,1)=zone(nz).Rknot(1,1);
                zone(k).Zknot(zone(k).Nx+1,1)=zone(nz).Zknot(1,1);
                zone(k).nknot(zone(k).Nx+1,1)=zone(nz).nknot(1,1);
            end
        end
        if(nex==0)
            zone(k).Rknot(zone(k).Nx+1,1)=zone(k).Rr(zone(k).Nx+1,1);
            zone(k).Zknot(zone(k).Nx+1,1)=zone(k).Zr(zone(k).Nx+1,1);
            nknots=nknots+1;
            zone(k).nknot(zone(k).Nx+1,1)=nknots;
        end
    else
        if(zone(k).V2(2)~=0)
            zone(k).Rknot(zone(k).Nx+1,1)=RN(zone(k).V2(2));
            zone(k).Zknot(zone(k).Nx+1,1)=ZN(zone(k).V2(2));
            zone(k).nknot(zone(k).Nx+1,1)=Nnknot(zone(k).V2(2));
        end
        if(zone(k).V(2)~=0)
            zone(k).Rknot(zone(k).Nx+1,1)=RX(zone(k).V(2));
            zone(k).Zknot(zone(k).Nx+1,1)=ZX(zone(k).V(2));
            zone(k).nknot(zone(k).Nx+1,1)=Xnknot(zone(k).V(2));
        end
    end
    
    figure(1)
    hold on
    plot(zone(k).R,zone(k).Z,'b.')
    
end


% noeud qui restent (intersection normale de 4 zones)


knots=[];
R=[];
Z=[];
for i=1:nzones
    Nx=zone(i).Nx;
    Nz=zone(i).Nz;
    for j=1:Nx+1
        for k=1:Nz+1
            if(zone(i).nknot(j,k)>0)
                if(length(find(knots==zone(i).nknot(j,k)))==0)
                    knots=[knots,zone(i).nknot(j,k)];
                    R=[R,zone(i).Rknot(j,k)];
                    Z=[Z,zone(i).Zknot(j,k)];
                end
            end
        end
    end
end


figure(1)
hold on
plot(R,Z,'rx')
axis equal


% on sauvegarde le numero des noeud eirene dans les coins des mailles
% soledge
for k=1:nzones
    zone(k).knotA=zone(k).nknot(1:zone(k).Nx,1:zone(k).Nz);
    zone(k).knotB=zone(k).nknot(2:zone(k).Nx+1,1:zone(k).Nz);
    zone(k).knotC=zone(k).nknot(2:zone(k).Nx+1,2:zone(k).Nz+1);
    zone(k).knotD=zone(k).nknot(1:zone(k).Nx,2:zone(k).Nz+1);
end

[A,B]=sort(knots);
knots=A;
R=R(B);
Z=Z(B);
