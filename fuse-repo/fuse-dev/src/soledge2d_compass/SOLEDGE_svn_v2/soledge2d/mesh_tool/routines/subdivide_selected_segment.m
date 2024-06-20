global zone_sel;
global side_sel;
global Rwall;
global Zwall;
global zones;

if(side_sel==1)
    R=zones.zone(zone_sel).north.R;
    Z=zones.zone(zone_sel).north.Z;
    c1.num=1;
    c1.arc(1).x=R;
    c1.arc(1).y=Z;
    c2.num=1;
    c2.arc(1).x=Rwall;
    c2.arc(1).y=Zwall;
    X=intersect_contour(c1,c2);
    if(X.num>=1)
        d=zeros(X.num,1);
        if(inpolygon(c1.arc(1).x(1),c1.arc(1).y(1),Rwall,Zwall))
            for k=1:X.num
               d(k)=sqrt((X.x(k)- c1.arc(1).x(1))^2+(X.y(k)-c1.arc(1).y(1))^2);
            end
            a=find(d==min(d));
            X.num=1;
            X.x(1)=X.x(a);
            X.y(1)=X.y(a);
        else
            for k=1:X.num
               d(k)=sqrt((X.x(k)- c1.arc(1).x(end))^2+(X.y(k)-c1.arc(1).y(end))^2);
            end
            a=find(d==min(d));
            X.num=1;
            X.x(1)=X.x(a);
            X.y(1)=X.y(a);
        end
    end
    if(X.num==1)
        cin.x=R;
        cin.y=Z;
        p1.x=X.x(1);
        p1.y=X.y(1);
        cout=part_contour(cin,p1);
        subR1=cout.arc(1).x;
        subZ1=cout.arc(1).y;
        subR2=cout.arc(2).x;
        subZ2=cout.arc(2).y;
        %         plot(subR1,subZ1,'k-','LineWidth',2)
        %         plot(subR2,subZ2,'r-','LineWidth',2)
        zones.zone(zone_sel).northaligned=true;
        zones.zone(zone_sel).subNorth(1).R=subR1;
        zones.zone(zone_sel).subNorth(1).Z=subZ1;
        zones.zone(zone_sel).subNorth(2).R=subR2;
        zones.zone(zone_sel).subNorth(2).Z=subZ2;
        zones.zone(zone_sel).subNorth(1).ismeshed=false;
        zones.zone(zone_sel).subNorth(2).ismeshed=false;
        if(zones.zone(zone_sel).Neighbour.north>0)
            nz=zones.zone(zone_sel).Neighbour.north;
            zones.zone(nz).southaligned=true;
            zones.zone(nz).subSouth(1).R=subR1;
            zones.zone(nz).subSouth(1).Z=subZ1;
            zones.zone(nz).subSouth(2).R=subR2;
            zones.zone(nz).subSouth(2).Z=subZ2;
            zones.zone(nz).subSouth(1).ismeshed=false;
            zones.zone(nz).subSouth(2).ismeshed=false;
        end
    end
end

if(side_sel==2)
    R=zones.zone(zone_sel).south.R;
    Z=zones.zone(zone_sel).south.Z;
    c1.num=1;
    c1.arc(1).x=R;
    c1.arc(1).y=Z;
    c2.num=1;
    c2.arc(1).x=Rwall;
    c2.arc(1).y=Zwall;
    X=intersect_contour(c1,c2);
    if(X.num==1)
        cin.x=R;
        cin.y=Z;
        p1.x=X.x(1);
        p1.y=X.y(1);
        cout=part_contour(cin,p1);
        subR1=cout.arc(1).x;
        subZ1=cout.arc(1).y;
        subR2=cout.arc(2).x;
        subZ2=cout.arc(2).y;
        %         plot(subR1,subZ1,'k-','LineWidth',2)
        %         plot(subR2,subZ2,'r-','LineWidth',2)
        zones.zone(zone_sel).southaligned=true;
        zones.zone(zone_sel).subSouth(1).R=subR1;
        zones.zone(zone_sel).subSouth(1).Z=subZ1;
        zones.zone(zone_sel).subSouth(2).R=subR2;
        zones.zone(zone_sel).subSouth(2).Z=subZ2;
        zones.zone(zone_sel).subSouth(1).ismeshed=false;
        zones.zone(zone_sel).subSouth(2).ismeshed=false;
        if(zones.zone(zone_sel).Neighbour.south>0)
            nz=zones.zone(zone_sel).Neighbour.south;
            zones.zone(nz).northaligned=true;
            zones.zone(nz).subNorth(1).R=subR1;
            zones.zone(nz).subNorth(1).Z=subZ1;
            zones.zone(nz).subNorth(2).R=subR2;
            zones.zone(nz).subNorth(2).Z=subZ2;
            zones.zone(nz).subNorth(1).ismeshed=false;
            zones.zone(nz).subNorth(2).ismeshed=false;
        end
    end
    
end

