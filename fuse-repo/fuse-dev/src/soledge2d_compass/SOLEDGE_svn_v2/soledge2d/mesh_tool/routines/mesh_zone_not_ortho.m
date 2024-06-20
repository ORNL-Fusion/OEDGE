function mesh_zone_not_ortho(nz,refR,refZ)

global zones;
global megazone;
global r2D;
global z2D;
global flux2D;

long=length(refR);
dist=zeros(1,long);
for k=2:long
    dist(k)=dist(k-1)+sqrt((refR(k)-refR(k-1))^2+(refZ(k)-refZ(k-1))^2);
end
dist=dist/dist(end);

if(megazone.mz(zones.zone(nz).mz).ismeshed)
    psilist=megazone.mz(zones.zone(nz).mz).refpoints.psi;
    
    %first one
    Rm=zones.zone(nz).south.R;
    Zm=zones.zone(nz).south.Z;
    dm=zeros(size(Rm));
    for j=2:length(Rm)
        dm(j)=dm(j-1)+sqrt((Rm(j)-Rm(j-1))^2+(Zm(j)-Zm(j-1))^2);
    end
    dm=dm/dm(end);
    for j=1:length(dist)
        zones.zone(nz).gridR(1,j)=interp1(dm,Rm,dist(j));
        zones.zone(nz).gridZ(1,j)=interp1(dm,Zm,dist(j));
    end
    zones.zone(nz).gridR(1,1)=zones.zone(nz).south.R(1);
    zones.zone(nz).gridZ(1,1)=zones.zone(nz).south.Z(1);
    zones.zone(nz).gridR(1,length(dist))=zones.zone(nz).south.R(end);
    zones.zone(nz).gridZ(1,length(dist))=zones.zone(nz).south.Z(end);
    for i=2:length(psilist)-1
        c1=contour_better(r2D,z2D,flux2D,psilist(i));
        cE.num=1;
        cE.arc(1).x=zones.zone(nz).east.R;
        cE.arc(1).y=zones.zone(nz).east.Z;
        cW.num=1;
        cW.arc(1).x=zones.zone(nz).west.R;
        cW.arc(1).y=zones.zone(nz).west.Z;
        XE=intersect_contour(c1,cE);
        XW=intersect_contour(c1,cW);
        cnum=XE.arc1(1);
        cin.x=c1.arc(cnum).x;
        cin.y=c1.arc(cnum).y;
        p1.x=XE.x(1);
        p1.y=XE.y(1);
        p2.x=XW.x(1);
        p2.y=XW.y(1);
        cout=part_contour2(cin,p1,p2);
        d1=sqrt((cout.x(1)-zones.zone(nz).gridR(i-1,1))^2+...
            (cout.y(1)-zones.zone(nz).gridZ(i-1,1))^2);
        d2=sqrt((cout.x(end)-zones.zone(nz).gridR(i-1,1))^2+...
            (cout.y(end)-zones.zone(nz).gridZ(i-1,1))^2);
        if(d1>d2) %return
            cout.x=cout.x(end:-1:1);
            cout.y=cout.y(end:-1:1);
        end
        Rm=cout.x;
        Zm=cout.y;
        dm=zeros(size(Rm));
        for j=2:length(Rm)
            dm(j)=dm(j-1)+sqrt((Rm(j)-Rm(j-1))^2+(Zm(j)-Zm(j-1))^2);
        end
        dm=dm/dm(end);
        for j=1:length(dist)
            zones.zone(nz).gridR(i,j)=interp1(dm,Rm,dist(j));
            zones.zone(nz).gridZ(i,j)=interp1(dm,Zm,dist(j));
        end
    end
    %last one
    Rm=zones.zone(nz).north.R;
    Zm=zones.zone(nz).north.Z;
    dm=zeros(size(Rm));
    for j=2:length(Rm)
        dm(j)=dm(j-1)+sqrt((Rm(j)-Rm(j-1))^2+(Zm(j)-Zm(j-1))^2);
    end
    dm=dm/dm(end);
    for j=1:length(dist)
        zones.zone(nz).gridR(length(psilist),j)=interp1(dm,Rm,dist(j));
        zones.zone(nz).gridZ(length(psilist),j)=interp1(dm,Zm,dist(j));
    end
    zones.zone(nz).gridR(length(psilist),1)=zones.zone(nz).north.R(1);
    zones.zone(nz).gridZ(length(psilist),1)=zones.zone(nz).north.Z(1);
    zones.zone(nz).gridR(length(psilist),length(dist))=zones.zone(nz).north.R(end);
    zones.zone(nz).gridZ(length(psilist),length(dist))=zones.zone(nz).north.Z(end);
    
end

end