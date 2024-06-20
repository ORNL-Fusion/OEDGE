global megazone;
global Pmegazone;
global zones;
global r2D;
global z2D;
global flux2D;
global segments;

numwait=0;
for k=1:Pmegazone.num
    if(Pmegazone.mz(k).ismeshed)
        numwait=numwait+length(Pmegazone.mz(k).refpoints.R);
    end
end

for k=1:zones.num
%     zones.zone(k).gridR=[];
%     zones.zone(k).gridZ=[];
    t1=zones.zone(k).Xtype_east;
    t2=zones.zone(k).Xtype_west;
    if((t1==1)&&(t2==1))
        zones.zone(k).meshortho=true;
    else
        zones.zone(k).meshortho=false;
    end
end

for k=1:Pmegazone.num
    Pmegazone.mz(k).meshortho=false;
    for k2=1:length(Pmegazone.mz(k).list)
        if(zones.zone(Pmegazone.mz(k).list(k2)).meshortho)
            Pmegazone.mz(k).meshortho=true;
        end
    end
end

hwait=waitbar(0,'meshing... (coffee break!)');

step=0;
for k=1:Pmegazone.num
    
    if((Pmegazone.mz(k).ismeshed)&&(Pmegazone.mz(k).meshortho))
        
        R=Pmegazone.mz(k).refpoints.R;
        Z=Pmegazone.mz(k).refpoints.Z;
        
        psimin=zones.zone(Pmegazone.mz(k).list(1)).pA.coord(1);
        psimax=zones.zone(Pmegazone.mz(k).list(end)).pB.coord(1);
        
        step=step+1; %first one
        Pmegazone.mz(k).mesh(1).R=[];
        Pmegazone.mz(k).mesh(1).Z=[];
        for k2=1:length(Pmegazone.mz(k).list)
            Pmegazone.mz(k).mesh(1).R=[Pmegazone.mz(k).mesh(1).R,...
                zones.zone(Pmegazone.mz(k).list(k2)).west.R];
            Pmegazone.mz(k).mesh(1).Z=[Pmegazone.mz(k).mesh(1).Z,...
                zones.zone(Pmegazone.mz(k).list(k2)).west.Z];
        end
        c1.num=1;
        c1.arc(1).x=Pmegazone.mz(k).mesh(1).R;
        c1.arc(1).y=Pmegazone.mz(k).mesh(1).Z;
        x1=c1.arc(1).x;
        y1=c1.arc(1).y;
        d=sqrt((x1(2:end)-x1(1:end-1)).^2+(y1(2:end)-y1(1:end-1)).^2);
        x1=x1(find(d~=0)+1);
        y1=y1(find(d~=0)+1);
        c1.arc(1).x=x1;
        c1.arc(1).y=y1;
        for k3=1:length(Pmegazone.mz(k).list)
            nz=Pmegazone.mz(k).list(k3);
            if(megazone.mz(zones.zone(nz).mz).ismeshed)
                psilist=megazone.mz(zones.zone(nz).mz).refpoints.psi;
                for k4=2:length(psilist)-1
                    c2=contour_better(r2D,z2D,flux2D,psilist(k4));
                    X=intersect_contour(c1,c2);
                    if(X.num>1)
                        disp(['Warning double point in zone ',num2str(nz)])
                    end
                    zones.zone(nz).gridR(k4,1)=X.x(1);
                    zones.zone(nz).gridZ(k4,1)=X.y(1);
                end
                zones.zone(nz).gridR(1,1)=zones.zone(nz).west.R(1);
                zones.zone(nz).gridR(length(psilist),1)=zones.zone(nz).west.R(end);
                zones.zone(nz).gridZ(1,1)=zones.zone(nz).west.Z(1);
                zones.zone(nz).gridZ(length(psilist),1)=zones.zone(nz).west.Z(end);
            end
        end
        waitbar(step/numwait,hwait,'meshing... (coffee break!)');
        
        dro=(max(max(r2D))-min(min(r2D)));
        dpsi=max(max(flux2D))-min(min(flux2D));
        psim=min(min(flux2D));
        psiM=max(max(flux2D));
        for k2=2:length(R)-1
            step=step+1;
            if(Pmegazone.mz(k).meshchanged)
                [R1,Z1]=follow_grad(R(k2),Z(k2),max(psimin-abs(dpsi)*0.05,psim),r2D,z2D,flux2D,5e-3);
                [R2,Z2]=follow_grad(R(k2),Z(k2),min(psimax+abs(dpsi)*0.05,psiM),r2D,z2D,flux2D,5e-3);
                Rm=[R1(end:-1:2),R2];
                Zm=[Z1(end:-1:2),Z2];
                Pmegazone.mz(k).mesh(k2).R=Rm;
                Pmegazone.mz(k).mesh(k2).Z=Zm;
            else
                Rm=Pmegazone.mz(k).mesh(k2).R;
                Zm=Pmegazone.mz(k).mesh(k2).Z;
            end
            c1.num=1;
            c1.arc(1).x=Rm;
            c1.arc(1).y=Zm;
            x1=c1.arc(1).x;
            y1=c1.arc(1).y;
            d=sqrt((x1(2:end)-x1(1:end-1)).^2+(y1(2:end)-y1(1:end-1)).^2);
            x1=x1(find(d~=0)+1);
            y1=y1(find(d~=0)+1);
            c1.arc(1).x=x1;
            c1.arc(1).y=y1;
            for k3=1:length(Pmegazone.mz(k).list)
                nz=Pmegazone.mz(k).list(k3);
                if(zones.zone(nz).meshortho)
                    if(megazone.mz(zones.zone(nz).mz).ismeshed)
                        psilist=megazone.mz(zones.zone(nz).mz).refpoints.psi;
                        for k4=1:length(psilist)
                            c2=contour_better(r2D,z2D,flux2D,psilist(k4));
                            X=intersect_contour(c1,c2);
                            if(X.num==0)
                                disp('f')
                                figure
                                hold on
                                contour(r2D,z2D,flux2D,psilist(k4),'b-')
                                plot(x1,y1,'m-')
                            end
                            zones.zone(nz).gridR(k4,k2)=X.x(1);
                            zones.zone(nz).gridZ(k4,k2)=X.y(1);
                        end
                    end
                end
            end
            waitbar(step/numwait,hwait,'meshing... (coffee break!)');
        end
        
        step=step+1; %last one
        Pmegazone.mz(k).mesh(length(R)).R=[];
        Pmegazone.mz(k).mesh(length(R)).Z=[];
        for k2=1:length(Pmegazone.mz(k).list)
            Pmegazone.mz(k).mesh(length(R)).R=[Pmegazone.mz(k).mesh(length(R)).R,...
                zones.zone(Pmegazone.mz(k).list(k2)).east.R];
            Pmegazone.mz(k).mesh(length(R)).Z=[Pmegazone.mz(k).mesh(length(R)).Z,...
                zones.zone(Pmegazone.mz(k).list(k2)).east.Z];
        end
        c1.num=1;
        c1.arc(1).x=Pmegazone.mz(k).mesh(length(R)).R;
        c1.arc(1).y=Pmegazone.mz(k).mesh(length(R)).Z;
        x1=c1.arc(1).x;
        y1=c1.arc(1).y;
        d=sqrt((x1(2:end)-x1(1:end-1)).^2+(y1(2:end)-y1(1:end-1)).^2);
        x1=x1(find(d~=0)+1);
        y1=y1(find(d~=0)+1);
        c1.arc(1).x=x1;
        c1.arc(1).y=y1;
        for k3=1:length(Pmegazone.mz(k).list)
            nz=Pmegazone.mz(k).list(k3);
            if(megazone.mz(zones.zone(nz).mz).ismeshed)
                psilist=megazone.mz(zones.zone(nz).mz).refpoints.psi;
                for k4=2:length(psilist)-1
                    c2=contour_better(r2D,z2D,flux2D,psilist(k4));
                    X=intersect_contour(c1,c2);
                    zones.zone(nz).gridR(k4,length(R))=X.x(1);
                    zones.zone(nz).gridZ(k4,length(R))=X.y(1);
                end
                zones.zone(nz).gridR(1,end)=zones.zone(nz).east.R(1);
                zones.zone(nz).gridR(length(psilist),end)=zones.zone(nz).east.R(end);
                zones.zone(nz).gridZ(1,end)=zones.zone(nz).east.Z(1);
                zones.zone(nz).gridZ(length(psilist),end)=zones.zone(nz).east.Z(end);
            end
        end
        waitbar(step/numwait,hwait,'meshing... (coffee break!)');
        
    end
    
    Pmegazone.mz(k).meshchanged=false;
end



% mesh non ortho
for k=1:Pmegazone.num
    if(Pmegazone.mz(k).ismeshed)
        nl=0;
        for k1=1:length(Pmegazone.mz(k).list)
            nz=Pmegazone.mz(k).list(k1);
            if(zones.zone(nz).meshortho)
                nl=k1;
            end
        end
        if(nl>0)
            containsortho=true;
        else
            containsortho=false;
        end
        
        if(~containsortho)
            meshstart=Pmegazone.mz(k).refpoints.nz;
            meshstartB=Pmegazone.mz(k).refpoints.nzB; %1north 2south
            
            num=find(Pmegazone.mz(k).list==meshstart);
            if(meshstartB==1) %north
                refR=Pmegazone.mz(k).refpoints.R;
                refZ=Pmegazone.mz(k).refpoints.Z;
                for k1=num+1:length(Pmegazone.mz(k).list)
                    nz=Pmegazone.mz(k).list(k1);
                    if(megazone.mz(zones.zone(nz).mz).ismeshed)
                        if(zones.zone(nz).orthomeshchanged)
                            optimized_mesh_zone_not_ortho(nz,refR,refZ,1);
                            refR=zones.zone(nz).gridR(end,:);
                            refZ=zones.zone(nz).gridZ(end,:);
                            zones.zone(nz).orthomeshchanged=false;
                        end
                    else
                        break;
                    end
                end
                refR=Pmegazone.mz(k).refpoints.R;
                refZ=Pmegazone.mz(k).refpoints.Z;
                for k1=num:-1:1
                    nz=Pmegazone.mz(k).list(k1);
                    if(megazone.mz(zones.zone(nz).mz).ismeshed)
                        if(zones.zone(nz).orthomeshchanged)
                            optimized_mesh_zone_not_ortho(nz,refR,refZ,-1);
                            refR=zones.zone(nz).gridR(1,:);
                            refZ=zones.zone(nz).gridZ(1,:);
                            zones.zone(nz).orthomeshchanged=false;
                        end
                    else
                        break;
                    end
                end
            else
                refR=Pmegazone.mz(k).refpoints.R;
                refZ=Pmegazone.mz(k).refpoints.Z;
                for k1=num:length(Pmegazone.mz(k).list)
                    nz=Pmegazone.mz(k).list(k1);
                    if(megazone.mz(zones.zone(nz).mz).ismeshed)
                        if(zones.zone(nz).orthomeshchanged)
                            optimized_mesh_zone_not_ortho(nz,refR,refZ,1);
                            refR=zones.zone(nz).gridR(end,:);
                            refZ=zones.zone(nz).gridZ(end,:);
                            zones.zone(nz).orthomeshchanged=false;
                        end
                    else
                        break;
                    end
                end
                refR=Pmegazone.mz(k).refpoints.R;
                refZ=Pmegazone.mz(k).refpoints.Z;
                for k1=num-1:-1:1
                    nz=Pmegazone.mz(k).list(k1);
                    if(megazone.mz(zones.zone(nz).mz).ismeshed)
                        if(zones.zone(nz).orthomeshchanged)
                            optimized_mesh_zone_not_ortho(nz,refR,refZ,-1);
                            refR=zones.zone(nz).gridR(1,:);
                            refZ=zones.zone(nz).gridZ(1,:);
                            zones.zone(nz).orthomeshchanged=false;
                        end
                    else
                        break;
                    end
                end
            end
            
        else
            num=nl+1;
            nz=Pmegazone.mz(k).list(nl);
            refR=zones.zone(nz).gridR(end,:);
            refZ=zones.zone(nz).gridZ(end,:);
            for k1=num:length(Pmegazone.mz(k).list)
                nz=Pmegazone.mz(k).list(k1);
                if(megazone.mz(zones.zone(nz).mz).ismeshed)
                    if(zones.zone(nz).orthomeshchanged)
                        optimized_mesh_zone_not_ortho(nz,refR,refZ,1);
                        refR=zones.zone(nz).gridR(end,:);
                        refZ=zones.zone(nz).gridZ(end,:);
                        zones.zone(nz).orthomeshchanged=false;
                    end
                else
                    break;
                end
            end
        end
        step=step+length(Pmegazone.mz(k).refpoints.R);
        waitbar(step/numwait,hwait,'meshing... (coffee break!)');
    end
    
    Pmegazone.mz(k).meshchanged=false;
end


% save_align_table();
delete(hwait)