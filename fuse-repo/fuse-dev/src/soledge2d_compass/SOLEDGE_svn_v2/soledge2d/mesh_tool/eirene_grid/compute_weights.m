global zones
global eirene


for k=1:zones.num
    Nx=zones.zone(k).Nx;
    Nz=zones.zone(k).Nz;
    for i=1:Nx
        for j=1:Nz
            zones.zone(k).list_tri_e(i,j).nums=[];
            zones.zone(k).list_tri_w(i,j).nums=[];
            zones.zone(k).list_tri_e(i,j).weights=[];
            zones.zone(k).list_tri_w(i,j).weights=[];
            if((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i+1,j+2)==1)) % cas 01 East
                list=zone(k).list_tri(i,j).nums;
                m=length(list);
                bR=zone(k).R2(i+1,j+2)-zone(k).R2(i+1,j+1);
                bZ=zone(k).Z2(i+1,j+2)-zone(k).Z2(i+1,j+1);
                norm=sqrt(bR^2+bZ^2);
                bR_=bR/norm;
                bZ_=bZ/norm;
                weights=zeros(1,m);
                for n=1:m
                    ntri=list(n);
                    ntri_vrai=triangle(ntri).ntri;
                    if(triangle(ntri).side==1)
                        Rmed=(R(ntri_vrai,2)+R(ntri_vrai,1))*0.5;
                        dSR=-(Z(ntri_vrai,2)-Z(ntri_vrai,1))*2*pi*Rmed;
                        dSZ=(R(ntri_vrai,2)-R(ntri_vrai,1))*2*pi*Rmed;
                        triangle(ntri).Surf=sqrt(dSR^2+dSZ^2);
                        norm=sqrt(dSR^2+dSZ^2);
                        dSR_=dSR/norm;
                        dSZ_=dSZ/norm;
                        cR=(R(ntri_vrai,2)+R(ntri_vrai,1))*0.5;
                        cZ=(Z(ntri_vrai,2)+Z(ntri_vrai,1))*0.5;
                    end
                    if(triangle(ntri).side==2)
                        Rmed=(R(ntri_vrai,3)+R(ntri_vrai,2))*0.5;
                        dSR=-(Z(ntri_vrai,3)-Z(ntri_vrai,2))*2*pi*Rmed;
                        dSZ=(R(ntri_vrai,3)-R(ntri_vrai,2))*2*pi*Rmed;
                        triangle(ntri).Surf=sqrt(dSR^2+dSZ^2);
                        norm=sqrt(dSR^2+dSZ^2);
                        dSR_=dSR/norm;
                        dSZ_=dSZ/norm;
                        cR=(R(ntri_vrai,3)+R(ntri_vrai,2))*0.5;
                        cZ=(Z(ntri_vrai,3)+Z(ntri_vrai,2))*0.5;
                    end
                    if(triangle(ntri).side==3)
                        Rmed=(R(ntri_vrai,3)+R(ntri_vrai,1))*0.5;
                        dSR=-(Z(ntri_vrai,1)-Z(ntri_vrai,3))*2*pi*Rmed;
                        dSZ=(R(ntri_vrai,1)-R(ntri_vrai,3))*2*pi*Rmed;
                        triangle(ntri).Surf=sqrt(dSR^2+dSZ^2);
                        norm=sqrt(dSR^2+dSZ^2);
                        dSR_=dSR/norm;
                        dSZ_=dSZ/norm;
                        cR=(R(ntri_vrai,3)+R(ntri_vrai,1))*0.5;
                        cZ=(Z(ntri_vrai,3)+Z(ntri_vrai,1))*0.5;
                    end
                    weights(n)=max(-(bR*dSR+bZ*dSZ),0);
                end
                if(sum(weights)~=0)
                    weights=weights/sum(weights);
                    zones.zone(k).list_tri_e(i,j).nums=zone(k).list_tri(i,j).nums;
                    zones.zone(k).list_tri_e(i,j).weights=weights;
                else
                    zones.zone(k).list_tri_e(i,j).nums=zone(k).list_tri(i,j).nums;
                    zones.zone(k).list_tri_e(i,j).weights=zeros(size(weights));
                end
            end
            
            if((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i+1,j)==1)) % cas 10 west
                list=zone(k).list_tri(i,j).nums;
                m=length(list);
                bR=zone(k).R2(i+1,j)-zone(k).R2(i+1,j+1);
                bZ=zone(k).Z2(i+1,j)-zone(k).Z2(i+1,j+1);
                norm=sqrt(bR^2+bZ^2);
                bR_=bR/norm;
                bZ_=bZ/norm;
                weights=zeros(1,m);
                for n=1:m
                    ntri=list(n);
                    ntri_vrai=triangle(ntri).ntri;
                    if(triangle(ntri).side==1)
                        Rmed=(R(ntri_vrai,2)+R(ntri_vrai,1))*0.5;
                        dSR=-(Z(ntri_vrai,2)-Z(ntri_vrai,1))*2*pi*Rmed;
                        dSZ=(R(ntri_vrai,2)-R(ntri_vrai,1))*2*pi*Rmed;
                        triangle(ntri).Surf=sqrt(dSR^2+dSZ^2);
                        norm=sqrt(dSR^2+dSZ^2);
                        dSR_=dSR/norm;
                        dSZ_=dSZ/norm;
                        cR=(R(ntri_vrai,2)+R(ntri_vrai,1))*0.5;
                        cZ=(Z(ntri_vrai,2)+Z(ntri_vrai,1))*0.5;
                    end
                    if(triangle(ntri).side==2)
                        Rmed=(R(ntri_vrai,3)+R(ntri_vrai,2))*0.5;
                        dSR=-(Z(ntri_vrai,3)-Z(ntri_vrai,2))*2*pi*Rmed;
                        dSZ=(R(ntri_vrai,3)-R(ntri_vrai,2))*2*pi*Rmed;
                        triangle(ntri).Surf=sqrt(dSR^2+dSZ^2);
                        norm=sqrt(dSR^2+dSZ^2);
                        dSR_=dSR/norm;
                        dSZ_=dSZ/norm;
                        cR=(R(ntri_vrai,3)+R(ntri_vrai,2))*0.5;
                        cZ=(Z(ntri_vrai,3)+Z(ntri_vrai,2))*0.5;
                    end
                    if(triangle(ntri).side==3)
                        Rmed=(R(ntri_vrai,3)+R(ntri_vrai,1))*0.5;
                        dSR=-(Z(ntri_vrai,1)-Z(ntri_vrai,3))*2*pi*Rmed;
                        dSZ=(R(ntri_vrai,1)-R(ntri_vrai,3))*2*pi*Rmed;
                        triangle(ntri).Surf=sqrt(dSR^2+dSZ^2);
                        norm=sqrt(dSR^2+dSZ^2);
                        dSR_=dSR/norm;
                        dSZ_=dSZ/norm;
                        cR=(R(ntri_vrai,3)+R(ntri_vrai,1))*0.5;
                        cZ=(Z(ntri_vrai,3)+Z(ntri_vrai,1))*0.5;
                    end
                    weights(n)=max(-(bR*dSR+bZ*dSZ),0);
                end
                if(sum(weights)~=0)
                    weights=weights/sum(weights);
                    zones.zone(k).list_tri_w(i,j).nums=zone(k).list_tri(i,j).nums;
                    zones.zone(k).list_tri_w(i,j).weights=weights;
                else
                    zones.zone(k).list_tri_w(i,j).nums=zone(k).list_tri(i,j).nums;
                    zones.zone(k).list_tri_w(i,j).weights=zeros(size(weights));
                end
            end
        end
    end
end
