clear

load ../WEST_paroi.txt

nzones=16;

cd ..
fid=fopen('Neighbors.txt')
fid2=fopen('MagNeighbors.txt')
cd Magnetic_input
for i=1:nzones
    filename=strcat('R_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    R=load(filename);
    [m,p]=size(R);
    filename=strcat('Z_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    Z=load(filename);
    zone(i).R=R;
    zone(i).Z=Z;
    zone(i).Nx=m;
    zone(i).Nz=p;
    A=fscanf(fid,'%d',4);
    B=fscanf(fid2,'%d',4);
    zone(i).vois(1)=A(1);
    zone(i).vois(2)=A(2);
    zone(i).vois(3)=A(3);
    zone(i).vois(4)=A(4);
    zone(i).isvois(1)=B(1);
    zone(i).isvois(2)=B(2);
    zone(i).isvois(3)=B(3);
    zone(i).isvois(4)=B(4);
    zone(i).ntrinum=zeros(m+1,p+1);
    zone(i).Rtri=zeros(m+1,p+1);
    zone(i).Ztri=zeros(m+1,p+1);
    zone(i).V=zeros(1,4);
    zone(i).V2=zeros(1,4);
end
fclose(fid)
fclose(fid2)

zone(1).V(3)=1;
zone(2).V(4)=1;
zone(4).V(2)=1;
zone(5).V(1)=1;
zone(5).V(3)=2;
zone(6).V(4)=2;
zone(7).V(2)=2;
zone(8).V(1)=2;
zone(9).V(3)=2;
zone(10).V(4)=2;
zone(11).V(3)=1;
zone(12).V(2)=2;
zone(12).V(4)=1;
zone(13).V(1)=2;
zone(14).V(2)=1;
zone(15).V(1)=1;

zone(2).V2(3)=1;
zone(3).V2(4)=1;
zone(4).V2(3)=2;
zone(5).V2(2)=1;
zone(5).V2(4)=2;
zone(6).V2(1)=1;
zone(11).V2(2)=2;
zone(12).V2(1)=2;
zone(12).V2(3)=3;
zone(13).V2(4)=3;
zone(15).V2(2)=3;
zone(16).V2(1)=3;

% zone(12
figure(1)
hold on

X=[1 2];
Xnumtri=[1 2];
RX=zeros(1,2);
ZX=zeros(1,2);

N=[1 2 3];
Nnumtri=[3 4 5];
RN=zeros(1,3);
ZN=zeros(1,3);

ntri=5;

nPX=2;
nPN=3;

for n=1:nPX
    for i=1:nzones
        if(zone(i).V(1)==n)
            RX(n)=RX(n)+zone(i).R(1,1);
            ZX(n)=ZX(n)+zone(i).Z(1,1);
        end
        if(zone(i).V(2)==n)
            RX(n)=RX(n)+zone(i).R(zone(i).Nx,1);
            ZX(n)=ZX(n)+zone(i).Z(zone(i).Nx,1);
        end
        if(zone(i).V(3)==n)
            RX(n)=RX(n)+zone(i).R(zone(i).Nx,zone(i).Nz);
            ZX(n)=ZX(n)+zone(i).Z(zone(i).Nx,zone(i).Nz);
        end
        if(zone(i).V(4)==n)
            RX(n)=RX(n)+zone(i).R(1,zone(i).Nz);
            ZX(n)=ZX(n)+zone(i).Z(1,zone(i).Nz);
        end
    end
end
RX=RX/8;
ZX=ZX/8;


for n=1:nPN
    for i=1:nzones
        if(zone(i).V2(1)==n)
            RN(n)=RN(n)+zone(i).R(1,1);
            ZN(n)=ZN(n)+zone(i).Z(1,1);
        end
        if(zone(i).V2(2)==n)
            RN(n)=RN(n)+zone(i).R(zone(i).Nx,1);
            ZN(n)=ZN(n)+zone(i).Z(zone(i).Nx,1);
        end
        if(zone(i).V2(3)==n)
            RN(n)=RN(n)+zone(i).R(zone(i).Nx,zone(i).Nz);
            ZN(n)=ZN(n)+zone(i).Z(zone(i).Nx,zone(i).Nz);
        end
        if(zone(i).V2(4)==n)
            RN(n)=RN(n)+zone(i).R(1,zone(i).Nz);
            ZN(n)=ZN(n)+zone(i).Z(1,zone(i).Nz);
        end
    end
end
RN=RN/4;
ZN=ZN/4;


for i=1:nzones
    %detection des points X
    if(zone(i).V(1)~=0)
        zone(i).Rtri(1,1)=RX(zone(i).V(1));
        zone(i).Ztri(1,1)=ZX(zone(i).V(1));
        zone(i).ntrinum(1,1)=Xnumtri(zone(i).V(1));
    end
    if(zone(i).V(2)~=0)
        zone(i).Rtri(zone(i).Nx+1,1)=RX(zone(i).V(2));
        zone(i).Ztri(zone(i).Nx+1,1)=ZX(zone(i).V(2));
        zone(i).ntrinum(zone(i).Nx+1,1)=Xnumtri(zone(i).V(2));
    end
    if(zone(i).V(3)~=0)
        zone(i).Rtri(zone(i).Nx+1,zone(i).Nz+1)=RX(zone(i).V(3));
        zone(i).Ztri(zone(i).Nx+1,zone(i).Nz+1)=ZX(zone(i).V(3));
        zone(i).ntrinum(zone(i).Nx+1,zone(i).Nz+1)=Xnumtri(zone(i).V(3));
    end
    if(zone(i).V(4)~=0)
        zone(i).Rtri(1,zone(i).Nz+1)=RX(zone(i).V(4));
        zone(i).Ztri(1,zone(i).Nz+1)=ZX(zone(i).V(4));
        zone(i).ntrinum(1,zone(i).Nz+1)=Xnumtri(zone(i).V(4));
    end
    
    %detection des points normaux
    if(zone(i).V2(1)~=0)
        zone(i).Rtri(1,1)=RN(zone(i).V2(1));
        zone(i).Ztri(1,1)=ZN(zone(i).V2(1));
        zone(i).ntrinum(1,1)=Nnumtri(zone(i).V2(1));
    end
    if(zone(i).V2(2)~=0)
        zone(i).Rtri(zone(i).Nx+1,1)=RN(zone(i).V2(2));
        zone(i).Ztri(zone(i).Nx+1,1)=ZN(zone(i).V2(2));
        zone(i).ntrinum(zone(i).Nx+1,1)=Nnumtri(zone(i).V2(2));
    end
    if(zone(i).V2(3)~=0)
        zone(i).Rtri(zone(i).Nx+1,zone(i).Nz+1)=RN(zone(i).V2(3));
        zone(i).Ztri(zone(i).Nx+1,zone(i).Nz+1)=ZN(zone(i).V2(3));
        zone(i).ntrinum(zone(i).Nx+1,zone(i).Nz+1)=Nnumtri(zone(i).V2(3));
    end
    if(zone(i).V2(4)~=0)
        zone(i).Rtri(1,zone(i).Nz+1)=RN(zone(i).V2(4));
        zone(i).Ztri(1,zone(i).Nz+1)=ZN(zone(i).V2(4));
        zone(i).ntrinum(1,zone(i).Nz+1)=Nnumtri(zone(i).V2(4));
    end
    
    
    % au sud
    if(zone(i).isvois(2)==1)
        %les triangles ne sont pas encore créés au sud (bord du domaine)
        cd Neighs
        filename=strcat('South_',num2str(i,'%.3d'));
        filename=strcat(filename,'.txt');
        Neig=load(filename);
        RNeig=Neig(:,3);
        ZNeig=Neig(:,4);
        if(zone(i).isvois(3)==1)
            % treat corner
            R1=2*RNeig(end)-RNeig(end-1);
            Z1=2*ZNeig(end)-ZNeig(end-1);
            filename=strcat('East_',num2str(i,'%.3d'));
            filename=strcat(filename,'.txt');
            Neig2=load(filename);
            RNeig2=Neig2(:,3);
            ZNeig2=Neig2(:,4);
            R2=2*RNeig2(1)-RNeig2(2);
            Z2=2*ZNeig2(1)-ZNeig2(2);
            zone(i).Rtri(1,zone(i).Nz+1)=(R1+R2)/2;
            zone(i).Ztri(1,zone(i).Nz+1)=(Z1+Z2)/2;
            zone(i).ntrinum(1,zone(i).Nz+1)=ntri+1;
            ntri=ntri+1;
            plot((R1+R2)/2,(Z1+Z2)/2,'co')
        else
            if(zone(i).vois(3)<i)
                %copy knots
                zn=zone(i).vois(3);
                zone(i).Rtri(1,zone(i).Nz+1)=zone(zn).Rtri(1,1);
                zone(i).Ztri(1,zone(i).Nz+1)=zone(zn).Ztri(1,1);
                zone(i).ntrinum(1,zone(i).Nz+1)=zone(zn).ntrinum(1,1);
            else
                if(zone(i).vois(3)>i)
                    %create knots
                    zn=zone(i).vois(3);
                    filename=strcat('South_',num2str(zn,'%.3d'));
                    filename=strcat(filename,'.txt');
                    Neig2=load(filename);
                    RNeig2=Neig2(:,3);
                    ZNeig2=Neig2(:,4);
                    zone(i).Rtri(1,zone(i).Nz+1)=(RNeig(end)+RNeig2(1))/2;
                    zone(i).Ztri(1,zone(i).Nz+1)=(ZNeig(end)+ZNeig2(1))/2;
                    zone(i).ntrinum(1,zone(i).Nz+1)=ntri+1;
                    ntri=ntri+1;
                end
            end
        end
        if(zone(i).isvois(4)==1)
            % treat corner
            R1=2*RNeig(1)-RNeig(2);
            Z1=2*ZNeig(1)-ZNeig(2);
            filename=strcat('West_',num2str(i,'%.3d'));
            filename=strcat(filename,'.txt');
            Neig2=load(filename);
            RNeig2=Neig2(:,3);
            ZNeig2=Neig2(:,4);
            R2=2*RNeig2(1)-RNeig2(2);
            Z2=2*ZNeig2(1)-ZNeig2(2);
            zone(i).Rtri(1,1)=(R2+R1)/2;
            zone(i).Ztri(1,1)=(Z2+Z1)/2;
            zone(i).ntrinum(1,1)=ntri+1;
            ntri=ntri+1;
            plot((R1+R2)/2,(Z1+Z2)/2,'co')
        else
            if(zone(i).vois(4)<i)
                %copy knots
                zn=zone(i).vois(4);
                zone(i).Rtri(1,1)=zone(zn).Rtri(1,zone(zn).Nz+1);
                zone(i).Ztri(1,1)=zone(zn).Ztri(1,zone(zn).Nz+1);
                zone(i).ntrinum(1,1)=zone(zn).ntrinum(1,zone(zn).Nz+1);
            else
                if(zone(i).vois(4)==i)
                    %create knots
                    zn=zone(i).vois(4);
                    filename=strcat('South_',num2str(zn,'%.3d'));
                    filename=strcat(filename,'.txt');
                    Neig2=load(filename);
                    RNeig2=Neig2(:,3);
                    ZNeig2=Neig2(:,4);
                    zone(i).Rtri(1,1)=(RNeig(1)+RNeig2(end))/2;
                    zone(i).Ztri(1,1)=(ZNeig(1)+ZNeig2(end))/2;
                    zone(i).ntrinum(1,1)=ntri+1;
                    ntri=ntri+1;
                else
                    %create knots
                    zn=zone(i).vois(4);
                    filename=strcat('South_',num2str(zn,'%.3d'));
                    filename=strcat(filename,'.txt');
                    Neig2=load(filename);
                    RNeig2=Neig2(:,3);
                    ZNeig2=Neig2(:,4);
                    zone(i).Rtri(1,1)=(RNeig(1)+RNeig2(end))/2;
                    zone(i).Ztri(1,1)=(ZNeig(1)+ZNeig2(end))/2;
                    zone(i).ntrinum(1,1)=ntri+1;
                    ntri=ntri+1;
                end
            end
        end
        
        for k=1:zone(i).Nz-1
            zone(i).Rtri(1,k+1)=(RNeig(k)+RNeig(k+1))*0.5;
            zone(i).Ztri(1,k+1)=(ZNeig(k)+ZNeig(k+1))*0.5;
            zone(i).ntrinum(1,k+1)=ntri+1;
            ntri=ntri+1;
        end
        cd ..
        figure(1)
        hold on
        plot(RNeig,ZNeig,'g.')
        
    else
        if(zone(i).vois(2)<i)
            %les triangles sont deja faits
            for k=0:zone(i).Nz
                zone(i).ntrinum(1,k+1)=zone(zone(i).vois(2)).ntrinum(zone(zone(i).vois(2)).Nx+1,k+1);
                zone(i).Rtri(1,k+1)=zone(zone(i).vois(2)).Rtri(zone(zone(i).vois(2)).Nx+1,k+1);
                zone(i).Ztri(1,k+1)=zone(zone(i).vois(2)).Ztri(zone(zone(i).vois(2)).Nx+1,k+1);
            end
        else
            %treating corners
            if(zone(i).isvois(3)==1)
                % treat corner
                zn=zone(i).vois(2);
                filename=strcat('East_',num2str(i,'%.3d'));
                filename=strcat(filename,'.txt');
                Neig=load(filename);
                RNeig=Neig(:,3);
                ZNeig=Neig(:,4);
                filename=strcat('East_',num2str(zn,'%.3d'));
                filename=strcat(filename,'.txt');
                Neig2=load(filename);
                RNeig2=Neig2(:,3);
                ZNeig2=Neig2(:,4);
                zone(i).Rtri(1,zone(i).Nz+1)=(RNeig(1)+RNeig2(end))/2;
                zone(i).Ztri(1,zone(i).Nz+1)=(ZNeig(1)+ZNeig2(end))/2;
                zone(i).ntrinum(1,zone(i).Nz+1)=ntri+1;
                ntri=ntri+1;
            end
            if(zone(i).isvois(4)==1)
                % treat corner
                zn=zone(i).vois(2);
                filename=strcat('West_',num2str(i,'%.3d'));
                filename=strcat(filename,'.txt');
                Neig=load(filename);
                RNeig=Neig(:,3);
                ZNeig=Neig(:,4);
                filename=strcat('West_',num2str(zn,'%.3d'));
                filename=strcat(filename,'.txt');
                Neig2=load(filename);
                RNeig2=Neig2(:,3);
                ZNeig2=Neig2(:,4);
                zone(i).Rtri(1,1)=(RNeig(1)+RNeig2(end))/2;
                zone(i).Ztri(1,1)=(ZNeig(1)+ZNeig2(end))/2;
                zone(i).ntrinum(1,1)=ntri+1;
                ntri=ntri+1;
            end
            
            %les triangles sont a faire
            Nxn=zone(zone(i).vois(2)).Nx;
            zn=zone(i).vois(2);
            for k=1:zone(i).Nz-1
                zone(i).Rtri(1,k+1)=(zone(i).R(1,k)+zone(i).R(1,k+1)+zone(zn).R(Nxn,k)+zone(zn).R(Nxn,k+1))/4;
                zone(i).Ztri(1,k+1)=(zone(i).Z(1,k)+zone(i).Z(1,k+1)+zone(zn).Z(Nxn,k)+zone(zn).Z(Nxn,k+1))/4;
                zone(i).ntrinum(1,k+1)=ntri+1;
                ntri=ntri+1;
            end
        end
    end
    
    % ############   au nord des zones  ##########
    if(zone(i).isvois(1)==1)
        %les triangles ne sont pas encore créés au nord (bord du domaine)
        cd Neighs
        filename=strcat('North_',num2str(i,'%.3d'));
        filename=strcat(filename,'.txt');
        Neig=load(filename);
        RNeig=Neig(:,3);
        ZNeig=Neig(:,4);
        Nx=zone(i).Nx;
        if(zone(i).isvois(3)==1)
            % treat corner
            R1=2*RNeig(end)-RNeig(end-1);
            Z1=2*ZNeig(end)-ZNeig(end-1);
            filename=strcat('East_',num2str(i,'%.3d'));
            filename=strcat(filename,'.txt');
            Neig2=load(filename);
            RNeig2=Neig2(:,3);
            ZNeig2=Neig2(:,4);
            R2=2*RNeig2(end)-RNeig2(end-1);
            Z2=2*ZNeig2(end)-ZNeig2(end-1);
            zone(i).Rtri(zone(i).Nx+1,zone(i).Nz+1)=(R1+R2)/2;
            zone(i).Ztri(zone(i).Nx+1,zone(i).Nz+1)=(Z1+Z2)/2;
            zone(i).ntrinum(zone(i).Nx+1,zone(i).Nz+1)=ntri+1;
            ntri=ntri+1;
            plot((R1+R2)/2,(Z1+Z2)/2,'co')
        else
            if(zone(i).vois(3)<i)
                %copy knots
                zn=zone(i).vois(3);
                zone(i).Rtri(zone(i).Nx+1,zone(i).Nz+1)=zone(zn).Rtri(zone(i).Nx+1,1);
                zone(i).Ztri(zone(i).Nx+1,zone(i).Nz+1)=zone(zn).Ztri(zone(i).Nx+1,1);
                zone(i).ntrinum(zone(i).Nx+1,zone(i).Nz+1)=zone(zn).ntrinum(zone(i).Nx+1,1);
            else
                if(zone(i).vois(3)>i)
                    %create knots
                    zn=zone(i).vois(3);
                    filename=strcat('North_',num2str(zn,'%.3d'));
                    filename=strcat(filename,'.txt');
                    Neig2=load(filename);
                    RNeig2=Neig2(:,3);
                    ZNeig2=Neig2(:,4);
                    zone(i).Rtri(zone(i).Nx+1,zone(i).Nz+1)=(RNeig(end)+RNeig2(1))/2;
                    zone(i).Ztri(zone(i).Nx+1,zone(i).Nz+1)=(ZNeig(end)+ZNeig2(1))/2;
                    zone(i).ntrinum(zone(i).Nx+1,zone(i).Nz+1)=ntri+1;
                    ntri=ntri+1;
                end
            end
        end
        if(zone(i).isvois(4)==1)
            % treat corner
            R1=2*RNeig(1)-RNeig(2);
            Z1=2*ZNeig(1)-ZNeig(2);
            filename=strcat('West_',num2str(i,'%.3d'));
            filename=strcat(filename,'.txt');
            Neig2=load(filename);
            RNeig2=Neig2(:,3);
            ZNeig2=Neig2(:,4);
            R2=2*RNeig2(end)-RNeig2(end-1);
            Z2=2*ZNeig2(end)-ZNeig2(end-1);
            zone(i).Rtri(zone(i).Nx+1,1)=(R2+R1)/2;
            zone(i).Ztri(zone(i).Nx+1,1)=(Z2+Z1)/2;
            zone(i).ntrinum(zone(i).Nx+1,1)=ntri+1;
            plot((R1+R2)/2,(Z1+Z2)/2,'co')
            ntri=ntri+1;
        else
            if(zone(i).vois(4)<i)
                %copy knots
                zn=zone(i).vois(4);
                zone(i).Rtri(zone(i).Nx+1,1)=zone(zn).Rtri(zone(i).Nx+1,zone(zn).Nz+1);
                zone(i).Ztri(zone(i).Nx+1,1)=zone(zn).Ztri(zone(i).Nx+1,zone(zn).Nz+1);
                zone(i).ntrinum(zone(i).Nx+1,1)=zone(zn).ntrinum(zone(i).Nx+1,zone(zn).Nz+1);
            else
                %create knots
                zn=zone(i).vois(4);
                filename=strcat('North_',num2str(zn,'%.3d'));
                filename=strcat(filename,'.txt');
                Neig2=load(filename);
                RNeig2=Neig2(:,3);
                ZNeig2=Neig2(:,4);
                zone(i).Rtri(zone(i).Nx+1,1)=(RNeig(1)+RNeig2(end))/2;
                zone(i).Ztri(zone(i).Nx+1,1)=(ZNeig(1)+ZNeig2(end))/2;
                zone(i).ntrinum(zone(i).Nx+1,1)=ntri+1;
                ntri=ntri+1;
            end
        end
        for k=1:zone(i).Nz-1
            zone(i).Rtri(Nx+1,k+1)=(RNeig(k)+RNeig(k+1))*0.5;
            zone(i).Ztri(Nx+1,k+1)=(ZNeig(k)+ZNeig(k+1))*0.5;
            zone(i).ntrinum(Nx+1,k+1)=ntri+1;
            ntri=ntri+1;
        end
        cd ..
        figure(1)
        hold on
        plot(RNeig,ZNeig,'g.')
        
    else
        if(zone(i).vois(1)<i)
            %les triangles sont deja faits
            Nx=zone(i).Nx;
            for k=0:zone(i).Nz
                zone(i).ntrinum(Nx+1,k+1)=zone(zone(i).vois(1)).ntrinum(1,k+1);
                zone(i).Rtri(Nx+1,k+1)=zone(zone(i).vois(1)).Rtri(1,k+1);
                zone(i).Ztri(Nx+1,k+1)=zone(zone(i).vois(1)).Ztri(1,k+1);
            end
        else
            %les triangles sont a faire
            %treating corners
            cd Neighs
            if(zone(i).isvois(3)==1)
                % treat corner
                zn=zone(i).vois(1);
                filename=strcat('East_',num2str(i,'%.3d'));
                filename=strcat(filename,'.txt');
                Neig=load(filename);
                RNeig=Neig(:,3);
                ZNeig=Neig(:,4);
                filename=strcat('East_',num2str(zn,'%.3d'));
                filename=strcat(filename,'.txt');
                Neig2=load(filename);
                RNeig2=Neig2(:,3);
                ZNeig2=Neig2(:,4);
                zone(i).Rtri(zone(i).Nx+1,zone(i).Nz+1)=(RNeig(end)+RNeig2(1))/2;
                zone(i).Ztri(zone(i).Nx+1,zone(i).Nz+1)=(ZNeig(end)+ZNeig2(1))/2;
                zone(i).ntrinum(zone(i).Nx+1,zone(i).Nz+1)=ntri+1;
                ntri=ntri+1;
            end
            if(zone(i).isvois(4)==1)
                % treat corner
                zn=zone(i).vois(1);
                filename=strcat('West_',num2str(i,'%.3d'));
                filename=strcat(filename,'.txt');
                Neig=load(filename);
                RNeig=Neig(:,3);
                ZNeig=Neig(:,4);
                filename=strcat('West_',num2str(zn,'%.3d'));
                filename=strcat(filename,'.txt');
                Neig2=load(filename);
                RNeig2=Neig2(:,3);
                ZNeig2=Neig2(:,4);
                zone(i).Rtri(zone(i).Nx+1,1)=(RNeig(end)+RNeig2(1))/2;
                zone(i).Ztri(zone(i).Nx+1,1)=(ZNeig(end)+ZNeig2(1))/2;
                zone(i).ntrinum(zone(i).Nx+1,1)=ntri+1;
                ntri=ntri+1;
            end
            cd ..
            Nx=zone(i).Nx;
            zn=zone(i).vois(1);
            for k=1:zone(i).Nz-1
                zone(i).Rtri(Nx+1,k+1)=(zone(i).R(Nx,k)+zone(i).R(Nx,k+1)+zone(zn).R(1,k)+zone(zn).R(1,k+1))/4;
                zone(i).Ztri(Nx+1,k+1)=(zone(i).Z(Nx,k)+zone(i).Z(Nx,k+1)+zone(zn).Z(1,k)+zone(zn).Z(1,k+1))/4;
                zone(i).ntrinum(Nx+1,k+1)=ntri+1;
                ntri=ntri+1;
            end
        end
    end
    
    %a l'est
    % au nord des zones
    if(zone(i).isvois(3)==1)
        %les triangles ne sont pas encore créés au nord (bord du domaine)
        cd Neighs
        filename=strcat('East_',num2str(i,'%.3d'));
        filename=strcat(filename,'.txt');
        Neig=load(filename);
        RNeig=Neig(:,3);
        ZNeig=Neig(:,4);
        Nx=zone(i).Nx;
        Nz=zone(i).Nz;
        for k=1:zone(i).Nx-1
            zone(i).Rtri(k+1,Nz+1)=(RNeig(k)+RNeig(k+1))*0.5;
            zone(i).Ztri(k+1,Nz+1)=(ZNeig(k)+ZNeig(k+1))*0.5;
            zone(i).ntrinum(k+1,Nz+1)=ntri+1;
            ntri=ntri+1;
        end
        cd ..
        figure(1)
        hold on
        plot(RNeig,ZNeig,'g.')
        
    else
        if(zone(i).vois(3)<i)
            %les triangles sont deja faits
            Nx=zone(i).Nx;
            Nz=zone(i).Nz;
            for k=1:zone(i).Nx-1
                zone(i).ntrinum(k+1,Nz+1)=zone(zone(i).vois(3)).ntrinum(k+1,1);
                zone(i).Rtri(k+1,Nz+1)=zone(zone(i).vois(3)).Rtri(k+1,1);
                zone(i).Ztri(k+1,Nz+1)=zone(zone(i).vois(3)).Ztri(k+1,1);
            end
        else
            if(zone(i).vois(3)>i)
                %les triangles sont a faire
                Nx=zone(i).Nx;
                Nz=zone(i).Nz
                zn=zone(i).vois(3);
                for k=1:zone(i).Nx-1
                    zone(i).Rtri(k+1,Nz+1)=(zone(i).R(k,Nz)+zone(i).R(k+1,Nz)+zone(zn).R(k,1)+zone(zn).R(k+1,1))/4;
                    zone(i).Ztri(k+1,Nz+1)=(zone(i).Z(k,Nz)+zone(i).Z(k+1,Nz)+zone(zn).Z(k,1)+zone(zn).Z(k+1,1))/4;
                    zone(i).ntrinum(k+1,Nz+1)=ntri+1;
                    ntri=ntri+1;
                end
            end
        end
    end
    
    %a l'ouest
    if(zone(i).isvois(4)==1)
        %les triangles ne sont pas encore créés au nord (bord du domaine)
        cd Neighs
        filename=strcat('West_',num2str(i,'%.3d'));
        filename=strcat(filename,'.txt');
        Neig=load(filename);
        RNeig=Neig(:,3);
        ZNeig=Neig(:,4);
        Nx=zone(i).Nx;
        Nz=zone(i).Nz;
        for k=1:zone(i).Nx-1
            zone(i).Rtri(k+1,1)=(RNeig(k)+RNeig(k+1))*0.5;
            zone(i).Ztri(k+1,1)=(ZNeig(k)+ZNeig(k+1))*0.5;
            zone(i).ntrinum(k+1,1)=ntri+1;
            ntri=ntri+1;
        end
        cd ..
        figure(1)
        hold on
        plot(RNeig,ZNeig,'g.')
        
    else
        if(zone(i).vois(4)<i)
            %les triangles sont deja faits
            Nx=zone(i).Nx;
            Nz=zone(i).Nz;
            Nzn=zone(zone(i).vois(4)).Nz
            for k=1:zone(i).Nx-1
                zone(i).ntrinum(k+1,1)=zone(zone(i).vois(4)).ntrinum(k+1,Nzn+1);
                zone(i).Rtri(k+1,1)=zone(zone(i).vois(4)).Rtri(k+1,Nzn+1);
                zone(i).Ztri(k+1,1)=zone(zone(i).vois(4)).Ztri(k+1,Nzn+1);
            end
        else
            %les triangles sont a faire
            Nx=zone(i).Nx;
            Nz=zone(i).Nz
            zn=zone(i).vois(4);
            Nzn=zone(zn).Nz
            for k=1:zone(i).Nx-1
                zone(i).Rtri(k+1,1)=(zone(i).R(k,1)+zone(i).R(k+1,1)+zone(zn).R(k,Nzn)+zone(zn).R(k+1,Nzn))/4;
                zone(i).Ztri(k+1,1)=(zone(i).Z(k,1)+zone(i).Z(k+1,1)+zone(zn).Z(k,Nzn)+zone(zn).Z(k+1,Nzn))/4;
                zone(i).ntrinum(k+1,1)=ntri+1;
                ntri=ntri+1;
            end
        end
    end
    
    %triangles at east when autoperiodic zone
    if(zone(i).vois(3)==i)
        Nx=zone(i).Nx;
        Nz=zone(i).Nz;
        for k=1:zone(i).Nx+1
            zone(i).ntrinum(k,Nz+1)=zone(zone(i).vois(3)).ntrinum(k,1);
            zone(i).Rtri(k,Nz+1)=zone(zone(i).vois(3)).Rtri(k,1);
            zone(i).Ztri(k,Nz+1)=zone(zone(i).vois(3)).Ztri(k,1);
        end
    end
    
    %knots for the middle
    Nx=zone(i).Nx;
    Nz=zone(i).Nz;
    for j=1:Nx-1
        for k=1:Nz-1
            zone(i).Rtri(j+1,k+1)=(zone(i).R(j,k)+zone(i).R(j,k+1)+zone(i).R(j+1,k)+zone(i).R(j+1,k+1))/4;
            zone(i).Ztri(j+1,k+1)=(zone(i).Z(j,k)+zone(i).Z(j,k+1)+zone(i).Z(j+1,k)+zone(i).Z(j+1,k+1))/4;
            zone(i).ntrinum(j+1,k+1)=ntri+1;
            ntri=ntri+1;
        end
    end
    figure(1)
    hold on
    plot(zone(i).R,zone(i).Z,'b.')
    
end



% noeud qui restent (intersection normale de 4 zones)


tri=[];
R=[];
Z=[];
for i=1:nzones
    Nx=zone(i).Nx;
    Nz=zone(i).Nz;
    for j=1:Nx+1
        for k=1:Nz+1
            if(zone(i).ntrinum(j,k)>0)
                if(length(find(tri==zone(i).ntrinum(j,k)))==0)
                    tri=[tri,zone(i).ntrinum(j,k)];
                    R=[R,zone(i).Rtri(j,k)];
                    Z=[Z,zone(i).Ztri(j,k)];
                end
            end
        end
    end
end


figure(1)
hold on
plot(R,Z,'rx')
axis equal

for k=1:nzones
    zone(k).knotA=zone(k).ntrinum(1:zone(k).Nx,1:zone(k).Nz);
    zone(k).knotB=zone(k).ntrinum(2:zone(k).Nx+1,1:zone(k).Nz);
    zone(k).knotC=zone(k).ntrinum(2:zone(k).Nx+1,2:zone(k).Nz+1);
    zone(k).knotD=zone(k).ntrinum(1:zone(k).Nx,2:zone(k).Nz+1);
end

[A,B]=sort(tri);
tri=A;
R=R(B);
Z=Z(B);

tri_w=[];
R_w=[];
Z_w=[];
tri_wn=[];
ntri_w=0
for k=1:nzones
    zone(k).iscrossed=zeros(zone(k).Nx,zone(k).Nz);
    zone(k).knotE=zeros(zone(k).Nx,zone(k).Nz);
    zone(k).knotF=zeros(zone(k).Nx,zone(k).Nz);
    zone(k).knotG=zeros(zone(k).Nx,zone(k).Nz);
    zone(k).knotH=zeros(zone(k).Nx,zone(k).Nz);
    k
    for i=1:zone(k).Nx
        for j=1:zone(k).Nz
            %segment1
            P1_R=R(zone(k).knotA(i,j));
            P1_Z=Z(zone(k).knotA(i,j));
            P2_R=R(zone(k).knotB(i,j));
            P2_Z=Z(zone(k).knotB(i,j));
            if(inpolygon(P1_R,P1_Z,WEST_paroi(:,1),WEST_paroi(:,2))-inpolygon(P2_R,P2_Z,WEST_paroi(:,1),WEST_paroi(:,2))~=0)
                zone(k).iscrossed(i,j)=1;
                ntri_w=ntri_w+1;
                tri_w=[tri_w,ntri_w];
                zone(k).knotE(i,j)=ntri_w;
                [r_w,z_w]=intersect_pol(P1_R,P1_Z,P2_R,P2_Z,WEST_paroi(:,1),WEST_paroi(:,2));
                R_w=[R_w,r_w];
                Z_w=[Z_w,z_w];
            end
            %segment2
            P1_R=R(zone(k).knotB(i,j));
            P1_Z=Z(zone(k).knotB(i,j));
            P2_R=R(zone(k).knotC(i,j));
            P2_Z=Z(zone(k).knotC(i,j));
            if(inpolygon(P1_R,P1_Z,WEST_paroi(:,1),WEST_paroi(:,2))-inpolygon(P2_R,P2_Z,WEST_paroi(:,1),WEST_paroi(:,2))~=0)
                zone(k).iscrossed(i,j)=1;
                ntri_w=ntri_w+1;
                tri_w=[tri_w,ntri_w];
                zone(k).knotF(i,j)=ntri_w;
                [r_w,z_w]=intersect_pol(P1_R,P1_Z,P2_R,P2_Z,WEST_paroi(:,1),WEST_paroi(:,2));
                R_w=[R_w,r_w];
                Z_w=[Z_w,z_w];
            end
            %segment3
            P1_R=R(zone(k).knotC(i,j));
            P1_Z=Z(zone(k).knotC(i,j));
            P2_R=R(zone(k).knotD(i,j));
            P2_Z=Z(zone(k).knotD(i,j));
            if(inpolygon(P1_R,P1_Z,WEST_paroi(:,1),WEST_paroi(:,2))-inpolygon(P2_R,P2_Z,WEST_paroi(:,1),WEST_paroi(:,2))~=0)
                zone(k).iscrossed(i,j)=1;
                ntri_w=ntri_w+1;
                tri_w=[tri_w,ntri_w];
                zone(k).knotG(i,j)=ntri_w;
                [r_w,z_w]=intersect_pol(P1_R,P1_Z,P2_R,P2_Z,WEST_paroi(:,1),WEST_paroi(:,2));
                R_w=[R_w,r_w];
                Z_w=[Z_w,z_w];
            end
            %segment4
            P1_R=R(zone(k).knotD(i,j));
            P1_Z=Z(zone(k).knotD(i,j));
            P2_R=R(zone(k).knotA(i,j));
            P2_Z=Z(zone(k).knotA(i,j));
            if(inpolygon(P1_R,P1_Z,WEST_paroi(:,1),WEST_paroi(:,2))-inpolygon(P2_R,P2_Z,WEST_paroi(:,1),WEST_paroi(:,2))~=0)
                zone(k).iscrossed(i,j)=1;
                ntri_w=ntri_w+1;
                tri_w=[tri_w,ntri_w];
                zone(k).knotH(i,j)=ntri_w;
                [r_w,z_w]=intersect_pol(P1_R,P1_Z,P2_R,P2_Z,WEST_paroi(:,1),WEST_paroi(:,2));
                R_w=[R_w,r_w];
                Z_w=[Z_w,z_w];
            end
        end
    end
end

%renumérotation des points
npt=1;
for n=1:ntri_w
    tri_wn(n)=npt;
    npt=npt+1;
    for n2=1:n-1
        if(sqrt((R_w(n)-R_w(n2))^2+(Z_w(n)-Z_w(n2))^2)<1e-6)
            tri_wn(n)=tri_wn(n2);
            npt=npt-1
            break
        end
    end
end
ntri_w2=max(tri_wn);
ntri_w=ntri_w+ntri;

for n=1:ntri_w2
    A=find(tri_wn,n);
    tri=[tri,n+ntri];
    R=[R,R_w(A(1))];
    Z=[Z,Z_w(A(1))];
end

[A,B]=sort(tri);
tri=A;
R=R(B);
Z=Z(B);

for k=1:nzones
    for i=1:zone(k).Nx
        for j=1:zone(k).Nz
            if(zone(k).knotE(i,j)~=0)
                zone(k).knotE(i,j)=tri_wn(zone(k).knotE(i,j));
            end
            if(zone(k).knotF(i,j)~=0)
                zone(k).knotF(i,j)=tri_wn(zone(k).knotF(i,j));
            end
            if(zone(k).knotG(i,j)~=0)
                zone(k).knotG(i,j)=tri_wn(zone(k).knotG(i,j));
            end
            if(zone(k).knotH(i,j)~=0)
                zone(k).knotH(i,j)=tri_wn(zone(k).knotH(i,j));
            end
        end
    end
end

plot(R_w,Z_w,'mv')