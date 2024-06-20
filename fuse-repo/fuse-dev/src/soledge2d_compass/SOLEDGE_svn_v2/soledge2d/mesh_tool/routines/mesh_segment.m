function mesh_segment(pt_num,refineType,refineSide,adjustMode,paramL,paramR)

global zones;
global megazone;
global Pmegazone;
global zone_sel;
global side_sel;
global r2D;
global z2D;
global flux2D;
global sub_sel;
global d_MI;

switch(side_sel)
    case 1
        if(zones.zone(zone_sel).northaligned)
            R=zones.zone(zone_sel).subNorth(sub_sel).R;
            Z=zones.zone(zone_sel).subNorth(sub_sel).Z;
        else
            R=zones.zone(zone_sel).north.R;
            Z=zones.zone(zone_sel).north.Z;
        end
    case 2
        if(zones.zone(zone_sel).southaligned)
            R=zones.zone(zone_sel).subSouth(sub_sel).R;
            Z=zones.zone(zone_sel).subSouth(sub_sel).Z;
        else
            R=zones.zone(zone_sel).south.R;
            Z=zones.zone(zone_sel).south.Z;
        end
    case 3
        R=zones.zone(zone_sel).east.R;
        Z=zones.zone(zone_sel).east.Z;
    case 4
        R=zones.zone(zone_sel).west.R;
        Z=zones.zone(zone_sel).west.Z;
end

if(side_sel>=3) %poloidal
    mz=zones.zone(zone_sel).mz;
    for k=1:length(megazone.mz(mz).list)
        zones.zone(megazone.mz(mz).list(k)).east.ismeshed=true;
        zones.zone(megazone.mz(mz).list(k)).west.ismeshed=true;
    end
    dist=zeros(size(R));
    for k=1:length(R)-1
        dist(k+1)=dist(k)+sqrt((R(k+1)-R(k))^2+(Z(k+1)-Z(k))^2) ;
    end
    %normalize
    longueur=dist(end);
    dist=dist/dist(end);
    if(refineType==0) %none
        pt=linspace(0,1,pt_num);
    else
        if(refineType==1) %linear
            if(refineSide==0) %left
                if(adjustMode==1) %absolute
                    M=zeros(2,2);
                    RHS=zeros(2,1);
                    M(1,1)=(pt_num-1)*pt_num/2;
                    M(1,2)=pt_num-1;
                    RHS(1)=1;
                    M(2,1)=1;
                    M(2,2)=1;
                    RHS(2)=paramL/100/longueur;
                    A=M^-1*RHS;
                    dpt=max(A(1)*(1:pt_num-1)+A(2),1e-3);
                    pt=zeros(1,pt_num);
                    for k=2:pt_num
                        pt(k)=pt(k-1)+dpt(k-1);
                    end
                    pt=(pt-pt(1))/(pt(end)-pt(1));
                else
                    dpt=1+paramL*(1:pt_num-1);
                    pt=zeros(1,pt_num);
                    for k=2:pt_num
                        pt(k)=pt(k-1)+dpt(k-1);
                    end
                    pt=(pt-pt(1))/(pt(end)-pt(1));
                end
            else
                if(refineSide==1) %right
                    if(adjustMode==1) %absolute
                        M=zeros(2,2);
                        RHS=zeros(2,1);
                        M(1,1)=(pt_num-1)*pt_num/2;
                        M(1,2)=pt_num-1;
                        RHS(1)=1;
                        M(2,1)=pt_num-1;
                        M(2,2)=1;
                        RHS(2)=paramR/100/longueur;
                        A=M^-1*RHS;
                        dpt=max(A(1)*(1:pt_num-1)+A(2),1e-3);
                        pt=zeros(1,pt_num);
                        for k=2:pt_num
                            pt(k)=pt(k-1)+dpt(k-1);
                        end
                        pt=(pt-pt(1))/(pt(end)-pt(1));
                    else
                        dpt=1+paramR*(1:pt_num-1);
                        pt=zeros(1,pt_num);
                        for k=2:pt_num
                            pt(k)=pt(k-1)+dpt(k-1);
                        end
                        pt=(pt-pt(1))/(pt(end)-pt(1));
                        pt=1-pt(end:-1:1);
                    end
                else %both
                    if(adjustMode==1) %absolute
                        %part1
                        M=zeros(2,2);
                        RHS=zeros(2,1);
                        M(1,1)=(pt_num/2-1)*pt_num/2/2;
                        M(1,2)=pt_num/2-1;
                        RHS(1)=0.5;
                        M(2,1)=1;
                        M(2,2)=1;
                        RHS(2)=paramL/100/longueur;
                        A=M^-1*RHS;
                        dpt=max(A(1)*(1:floor(pt_num/2)-1)+A(2),1e-3);
                        pt1=zeros(1,floor(pt_num/2));
                        for k=2:floor(pt_num/2)
                            pt1(k)=pt1(k-1)+dpt(k-1);
                        end
                        %part2
                        M=zeros(2,2);
                        RHS=zeros(2,1);
                        M(1,1)=(pt_num/2-1)*pt_num/2/2;
                        M(1,2)=pt_num/2-1;
                        RHS(1)=0.5;
                        M(2,1)=pt_num/2-1;
                        M(2,2)=1;
                        RHS(2)=paramR/100/longueur;
                        A=M^-1*RHS;
                        dpt=max(A(1)*(1:pt_num-floor(pt_num/2))+A(2),1e-3);
                        pt2=zeros(1,pt_num-floor(pt_num/2)+1);
                        for k=2:pt_num-floor(pt_num/2)+1
                            pt2(k)=pt2(k-1)+dpt(k-1);
                        end
                        pt=[pt1,0.5+pt2(2:end)];
                        pt=(pt-pt(1))/(pt(end)-pt(1));
                    else
                        dpt=1+paramL*(1:floor(pt_num/2)-1);
                        pt1=zeros(1,floor(pt_num/2));
                        for k=2:floor(pt_num/2)
                            pt1(k)=pt1(k-1)+dpt(k-1);
                        end
                        pt1=(pt1-pt1(1))/(pt1(end)-pt1(1))*0.5;
                        dpt=1+paramR*(1:pt_num-floor(pt_num/2));
                        pt2=zeros(1,pt_num-floor(pt_num/2)+1);
                        for k=2:pt_num-floor(pt_num/2)+1
                            pt2(k)=pt2(k-1)+dpt(k-1);
                        end
                        pt2=(pt2-pt2(1))/(pt2(end)-pt2(1));
                        pt2=(1-pt2(end:-1:1))*0.5;
                        pt=[pt1,0.5+pt2(2:end)];
                    end
                end
            end
        else %exponential
            waitfor(custom_mesh_interface(pt_num,R,Z));
            pt=d_MI;
        end
    end
    %check
    if(length(pt)~=pt_num)
        disp('problem')
    end
    if(min(abs(pt(2:end)-pt(1:end-1)))==0)
        disp('problem2')
    end
    megazone.mz(mz).refpoints.R=interp1(dist,R,pt);
    megazone.mz(mz).refpoints.Z=interp1(dist,Z,pt);
    megazone.mz(mz).refpoints.psi=interp2(r2D,z2D,flux2D,...
        megazone.mz(mz).refpoints.R,megazone.mz(mz).refpoints.Z);
    megazone.mz(mz).refpoints.psi(1)=zones.zone(zone_sel).pA.coord(1);
    megazone.mz(mz).refpoints.psi(end)=zones.zone(zone_sel).pB.coord(1);
    megazone.mz(mz).ismeshed=true;
    for nz=1:length(megazone.mz(mz).list)
        zones.zone(megazone.mz(mz).list(nz)).orthomeshchanged=true;
    end
end


if(side_sel<3) %radial
    mz=zones.zone(zone_sel).pmz;
    if(((side_sel==1)&&(zones.zone(zone_sel).northaligned))||...
            ((side_sel==2)&&(zones.zone(zone_sel).southaligned)))
        if(side_sel==1)
            zones.zone(zone_sel).subNorth(sub_sel).ismeshed=true;
            doit=true;
            for n=1:2
                doit=doit&&zones.zone(zone_sel).subNorth(n).ismeshed;
            end
        end
        if(side_sel==2)
            zones.zone(zone_sel).subSouth(sub_sel).ismeshed=true;
            doit=true;
            for n=1:2
                doit=doit&&zones.zone(zone_sel).subSouth(n).ismeshed;
            end
        end
        if(doit)
            for k=1:length(Pmegazone.mz(mz).list)
                zones.zone(Pmegazone.mz(mz).list(k)).north.ismeshed=true;
                zones.zone(Pmegazone.mz(mz).list(k)).south.ismeshed=true;
            end
            Pmegazone.mz(mz).ismeshed=true;
        end
    else
        for k=1:length(Pmegazone.mz(mz).list)
            zones.zone(Pmegazone.mz(mz).list(k)).north.ismeshed=true;
            zones.zone(Pmegazone.mz(mz).list(k)).south.ismeshed=true;
        end
        Pmegazone.mz(mz).ismeshed=true;
    end
    dist=zeros(size(R));
    for k=1:length(R)-1
        dist(k+1)=dist(k)+sqrt((R(k+1)-R(k))^2+(Z(k+1)-Z(k))^2) ;
    end
    dist_=dist(1);
    R_=R(1);
    Z_=Z(1);
    for np1=2:length(dist)
        if(abs(dist(np1)-dist(np1-1))>0)
            dist_=[dist_,dist(np1)];
            R_=[R_,R(np1)];
            Z_=[Z_,Z(np1)];
        end
    end
    dist=dist_;
    R=R_;
    Z=Z_;
    %normalize
    longueur=dist(end);
    dist=dist/dist(end);
    if(refineType==0) %none
        pt=linspace(0,1,pt_num);
    else
        if(refineType==1) %linear
            if(refineSide==0) %left
                if(adjustMode==1) %absolute
                    M=zeros(2,2);
                    RHS=zeros(2,1);
                    M(1,1)=(pt_num-1)*pt_num/2;
                    M(1,2)=pt_num-1;
                    RHS(1)=1;
                    M(2,1)=1;
                    M(2,2)=1;
                    RHS(2)=paramL/100/longueur;
                    A=M^-1*RHS;
                    dpt=max(A(1)*(1:pt_num-1)+A(2),1e-3);
                    pt=zeros(1,pt_num);
                    for k=2:pt_num
                        pt(k)=pt(k-1)+dpt(k-1);
                    end
                    pt=(pt-pt(1))/(pt(end)-pt(1));
                else
                    dpt=1+paramL*(1:pt_num-1);
                    pt=zeros(1,pt_num);
                    for k=2:pt_num
                        pt(k)=pt(k-1)+dpt(k-1);
                    end
                    pt=(pt-pt(1))/(pt(end)-pt(1));
                end
            else
                if(refineSide==1) %right
                    if(adjustMode==1) %absolute
                        M=zeros(2,2);
                        RHS=zeros(2,1);
                        M(1,1)=(pt_num-1)*pt_num/2;
                        M(1,2)=pt_num-1;
                        RHS(1)=1;
                        M(2,1)=pt_num-1;
                        M(2,2)=1;
                        RHS(2)=paramR/100/longueur;
                        A=M^-1*RHS;
                        dpt=max(A(1)*(1:pt_num-1)+A(2),1e-3);
                        pt=zeros(1,pt_num);
                        for k=2:pt_num
                            pt(k)=pt(k-1)+dpt(k-1);
                        end
                        pt=(pt-pt(1))/(pt(end)-pt(1));
                    else
                        dpt=1+paramR*(1:pt_num-1);
                        pt=zeros(1,pt_num);
                        for k=2:pt_num
                            pt(k)=pt(k-1)+dpt(k-1);
                        end
                        pt=(pt-pt(1))/(pt(end)-pt(1));
                        pt=1-pt(end:-1:1);
                    end
                else %both
                    if(adjustMode==1) %absolute
                        %part1
                        M=zeros(2,2);
                        RHS=zeros(2,1);
                        M(1,1)=(pt_num/2-1)*pt_num/2/2;
                        M(1,2)=pt_num/2-1;
                        RHS(1)=0.5;
                        M(2,1)=1;
                        M(2,2)=1;
                        RHS(2)=paramL/100/longueur;
                        A=M^-1*RHS;
                        dpt=max(A(1)*(1:floor(pt_num/2)-1)+A(2),1e-3);
                        pt1=zeros(1,floor(pt_num/2));
                        for k=2:floor(pt_num/2)
                            pt1(k)=pt1(k-1)+dpt(k-1);
                        end
                        %part2
                        M=zeros(2,2);
                        RHS=zeros(2,1);
                        M(1,1)=(pt_num/2-1)*pt_num/2/2;
                        M(1,2)=pt_num/2-1;
                        RHS(1)=0.5;
                        M(2,1)=pt_num/2-1;
                        M(2,2)=1;
                        RHS(2)=paramR/100/longueur;
                        A=M^-1*RHS;
                        dpt=max(A(1)*(1:pt_num-floor(pt_num/2))+A(2),1e-3);
                        pt2=zeros(1,pt_num-floor(pt_num/2)+1);
                        for k=2:pt_num-floor(pt_num/2)+1
                            pt2(k)=pt2(k-1)+dpt(k-1);
                        end
                        pt=[pt1,0.5+pt2(2:end)];
                        pt=(pt-pt(1))/(pt(end)-pt(1));
                    else
                        dpt=1+paramL*(1:floor(pt_num/2)-1);
                        pt1=zeros(1,floor(pt_num/2));
                        for k=2:floor(pt_num/2)
                            pt1(k)=pt1(k-1)+dpt(k-1);
                        end
                        pt1=(pt1-pt1(1))/(pt1(end)-pt1(1))*0.5;
                        dpt=1+paramR*(1:pt_num-floor(pt_num/2));
                        pt2=zeros(1,pt_num-floor(pt_num/2)+1);
                        for k=2:pt_num-floor(pt_num/2)+1
                            pt2(k)=pt2(k-1)+dpt(k-1);
                        end
                        pt2=(pt2-pt2(1))/(pt2(end)-pt2(1));
                        pt2=(1-pt2(end:-1:1))*0.5;
                        pt=[pt1,0.5+pt2(2:end)];
                    end
                end
            end
        else %exponential
            waitfor(custom_mesh_interface(pt_num,R,Z));
            pt=d_MI;
        end
    end
    %check
    if(length(pt)~=pt_num)
        disp('problem')
    end
    if(min(abs(pt(2:end)-pt(1:end-1)))==0)
        disp('problem2')
    end
    if(((side_sel==1)&&(zones.zone(zone_sel).northaligned))||...
            ((side_sel==2)&&(zones.zone(zone_sel).southaligned)))
        Pmegazone.mz(mz).subrefpoints(sub_sel).R=interp1(dist,R,pt);
        Pmegazone.mz(mz).subrefpoints(sub_sel).Z=interp1(dist,Z,pt);    
        Pmegazone.mz(mz).refpoints.R=[Pmegazone.mz(mz).subrefpoints(1).R(1:end-1),...
            Pmegazone.mz(mz).subrefpoints(2).R]
        Pmegazone.mz(mz).refpoints.Z=[Pmegazone.mz(mz).subrefpoints(1).Z(1:end-1),...
            Pmegazone.mz(mz).subrefpoints(2).Z]
    else
        Pmegazone.mz(mz).refpoints.R=interp1(dist,R,pt);
        Pmegazone.mz(mz).refpoints.Z=interp1(dist,Z,pt);
    end
    Pmegazone.mz(mz).refpoints.nz=zone_sel;
    Pmegazone.mz(mz).refpoints.nzB=side_sel;
    Pmegazone.mz(mz).meshchanged=true;
    for nz=1:length(Pmegazone.mz(mz).list)
        zones.zone(Pmegazone.mz(mz).list(nz)).orthomeshchanged=true;
    end
end


end