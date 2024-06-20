global eirene
global zones

for n=1:ntri_wall
    triangle(n).list_w=[];
    triangle(n).weight_w=[];
    triangle(n).list_e=[];
    triangle(n).weight_e=[];
    
    for k=1:zones.num
        Nx=zones.zone(k).Nx;
        Nz=zones.zone(k).Nz;
        for i=1:Nx
            for j=1:Nz
                if(((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i+1,j+2)==1))||...
                    ((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i+1,j)==1))) %cas 01 et 10 parallel
                    [a,b]=find(zones.zone(k).list_tri_e(i,j).nums==n);
                    if(length(b)==1)
                        triangle(n).list_e=[triangle(n).list_e,[i;j;k]];
                    end
                    triangle(n).weight_e=[triangle(n).weight_e,zones.zone(k).list_tri_e(i,j).weights(b)];
                end
            end
        end
    end
    
    for k=1:zones.num
        Nx=zones.zone(k).Nx;
        Nz=zones.zone(k).Nz;
        for i=1:Nx
            for j=1:Nz
                if(((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i+1,j+2)==1))||...
                    ((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i+1,j)==1))) %cas 01 et 10 parallel
                    [a,b]=find(zones.zone(k).list_tri_w(i,j).nums==n);
                    if(length(b)==1)
                        triangle(n).list_w=[triangle(n).list_w,[i;j;k]];
                    end
                    triangle(n).weight_w=[triangle(n).weight_w,zones.zone(k).list_tri_w(i,j).weights(b)];
                end
            end
        end
    end
    
    disp([num2str(n/ntri_wall*100,'%5.1f'),' %'])
    
end