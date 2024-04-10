global eirene
global zones

% on s'occupe des noeuds a l'intersection des mailles soledge et de la
% paroi
knot_w=[];
R_w=[];
Z_w=[];
knot_wn=[];
nknot_w=0;



for k=1:zones.num
    zones.zone(k).knotE=zeros(zones.zone(k).Nx,zones.zone(k).Nz);
    zones.zone(k).knotF=zeros(zones.zone(k).Nx,zones.zone(k).Nz);
    zones.zone(k).knotG=zeros(zones.zone(k).Nx,zones.zone(k).Nz);
    zones.zone(k).knotH=zeros(zones.zone(k).Nx,zones.zone(k).Nz);
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(~(zones.zone(k).isaligned(i,j)==1))
                if((sum(zones.zone(k).iscrossed(i,j,:))==2)&&(zones.zone(k).inplasma(i,j)==1))
                    %segment1
                    P1_R=eirene.R(zones.zone(k).knotA(i,j));
                    P1_Z=eirene.Z(zones.zone(k).knotA(i,j));
                    P2_R=eirene.R(zones.zone(k).knotB(i,j));
                    P2_Z=eirene.Z(zones.zone(k).knotB(i,j));
                    if(zones.zone(k).iscrossed(i,j,1)==1)
                        [r_w,z_w]=intersect_pol(P1_R,P1_Z,P2_R,P2_Z,eirene.Rwall,eirene.Zwall);
                        dist1=sqrt((P1_R-r_w)^2+(P1_Z-z_w)^2);
                        dist2=sqrt((P1_R-r_w)^2+(P1_Z-z_w)^2);
%                         if((dist1>1e-5)&&(dist2>1e-5)) % knot created if distinct
                            nknot_w=nknot_w+1;
                            knot_w=[knot_w,nknot_w];
                            zones.zone(k).knotE(i,j)=nknot_w;
                            R_w=[R_w,r_w];
                            Z_w=[Z_w,z_w];
%                         else
%                             disp(['I have done something ',num2str(dist1),' ',num2str(dist2)])
%                         end
                    end
                    %segment2
                    P1_R=eirene.R(zones.zone(k).knotB(i,j));
                    P1_Z=eirene.Z(zones.zone(k).knotB(i,j));
                    P2_R=eirene.R(zones.zone(k).knotC(i,j));
                    P2_Z=eirene.Z(zones.zone(k).knotC(i,j));
                    if(zones.zone(k).iscrossed(i,j,2)==1)
                        [r_w,z_w]=intersect_pol(P1_R,P1_Z,P2_R,P2_Z,eirene.Rwall,eirene.Zwall);
                        dist1=sqrt((P1_R-r_w)^2+(P1_Z-z_w)^2);
                        dist2=sqrt((P1_R-r_w)^2+(P1_Z-z_w)^2);
%                         if((dist1>1e-5)&&(dist2>1e-5)) % knot created if distinct
                            nknot_w=nknot_w+1;
                            knot_w=[knot_w,nknot_w];
                            zones.zone(k).knotF(i,j)=nknot_w;
                            R_w=[R_w,r_w];
                            Z_w=[Z_w,z_w];
%                         else
%                             disp(['I have done something ',num2str(dist1),' ',num2str(dist2)])
%                         end
                    end
                    %segment3
                    P1_R=eirene.R(zones.zone(k).knotC(i,j));
                    P1_Z=eirene.Z(zones.zone(k).knotC(i,j));
                    P2_R=eirene.R(zones.zone(k).knotD(i,j));
                    P2_Z=eirene.Z(zones.zone(k).knotD(i,j));
                    if(zones.zone(k).iscrossed(i,j,3)==1)
                        [r_w,z_w]=intersect_pol(P1_R,P1_Z,P2_R,P2_Z,eirene.Rwall,eirene.Zwall);
                        dist1=sqrt((P1_R-r_w)^2+(P1_Z-z_w)^2);
                        dist2=sqrt((P1_R-r_w)^2+(P1_Z-z_w)^2);
%                         if((dist1>1e-5)&&(dist2>1e-5)) % knot created if distinct
                            nknot_w=nknot_w+1;
                            knot_w=[knot_w,nknot_w];
                            zones.zone(k).knotG(i,j)=nknot_w;
                            R_w=[R_w,r_w];
                            Z_w=[Z_w,z_w];
%                         else
%                             disp(['I have done something ',num2str(dist1),' ',num2str(dist2)])
%                         end
                    end
                    %segment4
                    P1_R=eirene.R(zones.zone(k).knotD(i,j));
                    P1_Z=eirene.Z(zones.zone(k).knotD(i,j));
                    P2_R=eirene.R(zones.zone(k).knotA(i,j));
                    P2_Z=eirene.Z(zones.zone(k).knotA(i,j));
                    if(zones.zone(k).iscrossed(i,j,4)==1)
                        [r_w,z_w]=intersect_pol(P1_R,P1_Z,P2_R,P2_Z,eirene.Rwall,eirene.Zwall);
                        dist1=sqrt((P1_R-r_w)^2+(P1_Z-z_w)^2);
                        dist2=sqrt((P1_R-r_w)^2+(P1_Z-z_w)^2);
%                         if((dist1>1e-5)&&(dist2>1e-5)) % knot created if distinct
                            nknot_w=nknot_w+1;
                            knot_w=[knot_w,nknot_w];
                            zones.zone(k).knotH(i,j)=nknot_w;
                            R_w=[R_w,r_w];
                            Z_w=[Z_w,z_w];
%                         else
%                             disp(['I have done something ',num2str(dist1),' ',num2str(dist2)])
%                         end
                    end
                end
            end
        end
    end
end

%renumérotation des points
npt=1;
for n=1:nknot_w
    knot_wn(n)=npt;
    npt=npt+1;
    for n2=1:n-1
        if(sqrt((R_w(n)-R_w(n2))^2+(Z_w(n)-Z_w(n2))^2)<1e-6)
            knot_wn(n)=knot_wn(n2);
            npt=npt-1;
            break
        end
    end
end
nknot_w2=max(knot_wn);

Rnewknots=[];
Znewknots=[];
for n=1:nknot_w2
    A=find(knot_wn==n);
    knots=[knots,n+nknots];
    Rnewknots=[Rnewknots,R_w(A(1))];
    Znewknots=[Znewknots,Z_w(A(1))];
end


% on sauvegarde le numero des noeuds eirene a l'intersection des maille
% soledge et de la paroi.
for k=1:zones.num
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(zones.zone(k).knotE(i,j)~=0)
                zones.zone(k).knotE(i,j)=knot_wn(zones.zone(k).knotE(i,j))+eirene.nknots;
            end
            if(zones.zone(k).knotF(i,j)~=0)
                zones.zone(k).knotF(i,j)=knot_wn(zones.zone(k).knotF(i,j))+eirene.nknots;
            end
            if(zones.zone(k).knotG(i,j)~=0)
                zones.zone(k).knotG(i,j)=knot_wn(zones.zone(k).knotG(i,j))+eirene.nknots;
            end
            if(zones.zone(k).knotH(i,j)~=0)
                zones.zone(k).knotH(i,j)=knot_wn(zones.zone(k).knotH(i,j))+eirene.nknots;
            end
        end
    end
end

eirene.nknots=nknot_w2+eirene.nknots;
