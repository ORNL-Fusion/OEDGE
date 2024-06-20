global eirene
global zones

for n=1:ntri_wall
    triangle(n).list_n=[];
    triangle(n).weight_n=[];
    triangle(n).list_s=[];
    triangle(n).weight_s=[];
    
    for k=1:zones.num
        Nx=zones.zone(k).Nx;
        Nz=zones.zone(k).Nz;
        for i=1:Nx
            for j=1:Nz
                if(((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i+2,j+1)==1))||...
                    ((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i,j+1)==1))) %cas 01 et 10 parallel
                    [a,b]=find(zones.zone(k).list_tri_s(i,j).nums==n);
                    if(length(b)==1)
                        triangle(n).list_s=[triangle(n).list_s,[i;j;k]];
                    end
                    triangle(n).weight_s=[triangle(n).weight_s,zones.zone(k).list_tri_s(i,j).weights(b)];
                end
            end
        end
    end
    
    for k=1:zones.num
        Nx=zones.zone(k).Nx;
        Nz=zones.zone(k).Nz;
        for i=1:Nx
            for j=1:Nz
                if(((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i+2,j+1)==1))||...
                    ((zone(k).chi2(i+1,j+1)==0)&&(zone(k).chi2(i,j+1)==1))) %cas 01 et 10 parallel
                    [a,b]=find(zones.zone(k).list_tri_n(i,j).nums==n);
                    if(length(b)==1)
                        triangle(n).list_n=[triangle(n).list_n,[i;j;k]];
                    end
                    triangle(n).weight_n=[triangle(n).weight_n,zones.zone(k).list_tri_n(i,j).weights(b)];
                end
            end
        end
    end
    
    disp([num2str(n/ntri_wall*100,'%5.1f'),' %'])
    
end