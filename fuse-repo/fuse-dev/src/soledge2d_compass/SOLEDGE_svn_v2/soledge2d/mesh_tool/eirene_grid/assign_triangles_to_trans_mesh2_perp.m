global zones
global eirene
global Rwall
global Zwall


for k=1:zones.num
    Nx=zones.zone(k).Nx;
    Nz=zones.zone(k).Nz;
    for i=1:Nx
        for j=1:Nz
            if(((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i+2,j+1)==1))||...
                    ((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i,j+1)==1))) %cas 01 et 10 perp
                pointOK=0;
                ncase=0;
                while(pointOK==0)
                    ncase=ncase+3;
                    %mesh a considerer: %essayer avec n_case maille autour
                    list=[i;j;k];
                    i_=i;
                    j_=j;
                    k_=k;
                    for n=1:ncase
                        i_2=i_;
                        j_2=j_;
                        k_2=k_;
                        i_=zone(k_2).point_nord(i_2,j_2).i;
                        j_=zone(k_2).point_nord(i_2,j_2).j;
                        k_=zone(k_2).point_nord(i_2,j_2).k;
                        %test proximite
                        if(k_==0)
                            break
                        end
                        d=sqrt((zones.zone(k_2).gridRc(i_2,j_2)-zones.zone(k_).gridRc(i_,j_))^2+...
                            (zones.zone(k_2).gridZc(i_2,j_2)-zones.zone(k_).gridZc(i_,j_))^2);
                        list=[[i_;j_;k_],list];
                        if(d>0.1)
                            break
                        end
                    end
                    i_=i;
                    j_=j;
                    k_=k;
                    for n=1:ncase
                        i_2=i_;
                        j_2=j_;
                        k_2=k_;
                        i_=zone(k_2).point_sud(i_2,j_2).i;
                        j_=zone(k_2).point_sud(i_2,j_2).j;
                        k_=zone(k_2).point_sud(i_2,j_2).k;
                        %test proximite
                        if(k_==0)
                            break
                        end
                        d=sqrt((zones.zone(k_2).gridRc(i_2,j_2)-zones.zone(k_).gridRc(i_,j_))^2+...
                            (zones.zone(k_2).gridZc(i_2,j_2)-zones.zone(k_).gridZc(i_,j_))^2);
                        list=[list,[i_;j_;k_]];
                        if(d>0.1)
                            break
                        end
                    end
                                        
                    list_tri=[];
                    for n=1:ntri_wall
                        [m2,p2]=size(list);
                        for n2=1:p2
                            test=[triangle(n).i;triangle(n).j;triangle(n).k]==list(:,n2);
                            if(sum(test)==3)
                                list_tri=[list_tri,n];
                            end
                        end
                    end
                    
                    
                    m=length(list_tri);
                    if(m==0)
                        %planB
                        disp(['error no triangle found ',num2str(i),' ',num2str(j),' ',num2str(k),'  ncase=',num2str(ncase)])
                        if(ncase>50)
                            return %stop here
                        end
                    else
                        pointOK=1;    
                        zone(k).list_tri_perp(i,j).nums=list_tri;
                    end
                end
            end
        end
    end
end