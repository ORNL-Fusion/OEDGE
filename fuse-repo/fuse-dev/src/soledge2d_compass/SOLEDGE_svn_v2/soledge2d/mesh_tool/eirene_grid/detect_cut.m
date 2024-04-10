global eirene;
global zones;
% global Rwall;
% global Zwall;

for k=1:zones.num
    zones.zone(k).iscrossed=zeros(zones.zone(k).Nx,zones.zone(k).Nz,4);
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            ncut=0;
            %segment1
            P1_R=eirene.R(zones.zone(k).knotA(i,j));
            P1_Z=eirene.Z(zones.zone(k).knotA(i,j));
            P2_R=eirene.R(zones.zone(k).knotB(i,j));
            P2_Z=eirene.Z(zones.zone(k).knotB(i,j));
            if(inpolygon(P1_R,P1_Z,eirene.Rwall,eirene.Zwall)-inpolygon(P2_R,P2_Z,eirene.Rwall,eirene.Zwall)~=0)
                zones.zone(k).iscrossed(i,j,1)=1;
            end
            %segment2
            P1_R=eirene.R(zones.zone(k).knotB(i,j));
            P1_Z=eirene.Z(zones.zone(k).knotB(i,j));
            P2_R=eirene.R(zones.zone(k).knotC(i,j));
            P2_Z=eirene.Z(zones.zone(k).knotC(i,j));
            if(inpolygon(P1_R,P1_Z,eirene.Rwall,eirene.Zwall)-inpolygon(P2_R,P2_Z,eirene.Rwall,eirene.Zwall)~=0)
                zones.zone(k).iscrossed(i,j,2)=1;
            end
            %segment3
            P1_R=eirene.R(zones.zone(k).knotC(i,j));
            P1_Z=eirene.Z(zones.zone(k).knotC(i,j));
            P2_R=eirene.R(zones.zone(k).knotD(i,j));
            P2_Z=eirene.Z(zones.zone(k).knotD(i,j));
            if(inpolygon(P1_R,P1_Z,eirene.Rwall,eirene.Zwall)-inpolygon(P2_R,P2_Z,eirene.Rwall,eirene.Zwall)~=0)
                zones.zone(k).iscrossed(i,j,3)=1;
            end
            %             segment4
            P1_R=eirene.R(zones.zone(k).knotD(i,j));
            P1_Z=eirene.Z(zones.zone(k).knotD(i,j));
            P2_R=eirene.R(zones.zone(k).knotA(i,j));
            P2_Z=eirene.Z(zones.zone(k).knotA(i,j));
            if(inpolygon(P1_R,P1_Z,eirene.Rwall,eirene.Zwall)-inpolygon(P2_R,P2_Z,eirene.Rwall,eirene.Zwall)~=0)
                zones.zone(k).iscrossed(i,j,4)=1;
            end
        end
    end
end

for k=1:zones.num
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(sum(zones.zone(k).iscrossed(i,j,:))>=3)
                %needs refinement
                nres=10; % raffinement
                ncut1=0;
                %segment1
                P1_R=eirene.R(zones.zone(k).knotA(i,j));
                P1_Z=eirene.Z(zones.zone(k).knotA(i,j));
                P2_R=eirene.R(zones.zone(k).knotB(i,j));
                P2_Z=eirene.Z(zones.zone(k).knotB(i,j));
                %             [Rx,Zx]=polyxpoly([P1_R P2_R],[P1_Z,P2_Z],eirene.Rwall,eirene.Zwall);
                %             ncut=ncut+length(Rx);
                segR=linspace(P1_R,P2_R,nres);
                segZ=linspace(P1_Z,P2_Z,nres);
                for nst=1:nres-1
                    if(inpolygon(segR(nst),segZ(nst),eirene.Rwall,eirene.Zwall)-inpolygon(segR(nst+1),segZ(nst+1),eirene.Rwall,eirene.Zwall)~=0)
                        ncut1=ncut1+1;
                    end
                end
                %segment2
                ncut2=0;
                P1_R=eirene.R(zones.zone(k).knotB(i,j));
                P1_Z=eirene.Z(zones.zone(k).knotB(i,j));
                P2_R=eirene.R(zones.zone(k).knotC(i,j));
                P2_Z=eirene.Z(zones.zone(k).knotC(i,j));
                %             [Rx,Zx]=polyxpoly([P1_R P2_R],[P1_Z,P2_Z],eirene.Rwall,eirene.Zwall);
                %             ncut=ncut+length(Rx);
                segR=linspace(P1_R,P2_R,nres);
                segZ=linspace(P1_Z,P2_Z,nres);
                for nst=1:nres-1
                    if(inpolygon(segR(nst),segZ(nst),eirene.Rwall,eirene.Zwall)-inpolygon(segR(nst+1),segZ(nst+1),eirene.Rwall,eirene.Zwall)~=0)
                        ncut2=ncut2+1;
                    end
                end
                %segment3
                ncut3=0;
                P1_R=eirene.R(zones.zone(k).knotC(i,j));
                P1_Z=eirene.Z(zones.zone(k).knotC(i,j));
                P2_R=eirene.R(zones.zone(k).knotD(i,j));
                P2_Z=eirene.Z(zones.zone(k).knotD(i,j));
                %             [Rx,Zx]=polyxpoly([P1_R P2_R],[P1_Z,P2_Z],eirene.Rwall,eirene.Zwall);
                %             ncut=ncut+length(Rx);
                segR=linspace(P1_R,P2_R,nres);
                segZ=linspace(P1_Z,P2_Z,nres);
                for nst=1:nres-1
                    if(inpolygon(segR(nst),segZ(nst),eirene.Rwall,eirene.Zwall)-inpolygon(segR(nst+1),segZ(nst+1),eirene.Rwall,eirene.Zwall)~=0)
                        ncut3=ncut3+1;
                    end
                end
                %             segment4
                ncut4=0;
                P1_R=eirene.R(zones.zone(k).knotD(i,j));
                P1_Z=eirene.Z(zones.zone(k).knotD(i,j));
                P2_R=eirene.R(zones.zone(k).knotA(i,j));
                P2_Z=eirene.Z(zones.zone(k).knotA(i,j));
                %             [Rx,Zx]=polyxpoly([P1_R P2_R],[P1_Z,P2_Z],eirene.Rwall,eirene.Zwall);
                %             ncut=ncut+length(Rx);
                segR=linspace(P1_R,P2_R,nres);
                segZ=linspace(P1_Z,P2_Z,nres);
                for nst=1:nres-1
                    if(inpolygon(segR(nst),segZ(nst),eirene.Rwall,eirene.Zwall)-inpolygon(segR(nst+1),segZ(nst+1),eirene.Rwall,eirene.Zwall)~=0)
                        ncut4=ncut4+1;
                    end
                end
                %take decision
                if(ncut1==1)
                    zones.zone(k).iscrossed(i,j,1)=1;
                else
                    zones.zone(k).iscrossed(i,j,1)=0;
                end
                if(ncut2==1)
                    zones.zone(k).iscrossed(i,j,2)=1;
                else
                    zones.zone(k).iscrossed(i,j,2)=0;
                end
                if(ncut3==1)
                    zones.zone(k).iscrossed(i,j,3)=1;
                else
                    zones.zone(k).iscrossed(i,j,3)=0;
                end
                if(ncut4==1)
                    zones.zone(k).iscrossed(i,j,4)=1;
                else
                    zones.zone(k).iscrossed(i,j,4)=0;
                end
                if(sum(zones.zone(k).iscrossed(i,j,:))~=2)
                    %manual
                    hf=figure()
                    hold on
                    plot(eirene.Rwall,eirene.Zwall,'k-')
                    nA=zones.zone(k).knotA(i,j);
                    nB=zones.zone(k).knotB(i,j);
                    nC=zones.zone(k).knotC(i,j);
                    nD=zones.zone(k).knotD(i,j);
                    plot([eirene.R(nA) eirene.R(nB) eirene.R(nC) eirene.R(nD) eirene.R(nA)],...
                        [eirene.Z(nA) eirene.Z(nB) eirene.Z(nC) eirene.Z(nD) eirene.Z(nA)],'b.-');
                    text(eirene.R(nA),eirene.Z(nA),'A');
                    text(eirene.R(nB),eirene.Z(nB),'B');
                    text(eirene.R(nC),eirene.Z(nC),'C');
                    text(eirene.R(nD),eirene.Z(nD),'D');
                    hold off
                    choice1 = questdlg('Select segment 1','Select segment 1',...
                        'AB','BC','CD','AD','AB');
                    choice2 = questdlg('Select segment 2','Select segment 2 (must be different of segment 1)',...
                        'AB','BC','CD','AD','AB');
                    close(hf);
                    zones.zone(k).iscrossed(i,j,1)=0;
                    zones.zone(k).iscrossed(i,j,2)=0;
                    zones.zone(k).iscrossed(i,j,3)=0;
                    zones.zone(k).iscrossed(i,j,4)=0;
                    switch choice1
                        case 'AB'
                            zones.zone(k).iscrossed(i,j,1)=1;
                        case 'BC' 
                            zones.zone(k).iscrossed(i,j,2)=1;
                        case 'CD'
                            zones.zone(k).iscrossed(i,j,3)=1;
                        case 'AD'
                            zones.zone(k).iscrossed(i,j,4)=1;
                    end
                    switch choice2
                        case 'AB'
                            zones.zone(k).iscrossed(i,j,1)=1;
                        case 'BC' 
                            zones.zone(k).iscrossed(i,j,2)=1;
                        case 'CD'
                            zones.zone(k).iscrossed(i,j,3)=1;
                        case 'AD'
                            zones.zone(k).iscrossed(i,j,4)=1;
                    end
                end
            end
        end
    end
end