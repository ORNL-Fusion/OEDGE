global zones
global X_points
global Rwall
global Zwall


for k=1:zones.num
    
    if((zones.zone(k).pA.coord(3)==-1)&&(zones.zone(k).pA.coord(4)==1))
       if(~inpolygon(X_points.R(zones.zone(k).pA.coord(2)),...
               X_points.Z(zones.zone(k).pA.coord(2)),Rwall,Zwall))
           zones.zone(k).chi(1,1)=1;
       end
   end
   if((zones.zone(k).pB.coord(3)==-1)&&(zones.zone(k).pB.coord(4)==1))
       if(~inpolygon(X_points.R(zones.zone(k).pB.coord(2)),...
               X_points.Z(zones.zone(k).pB.coord(2)),Rwall,Zwall))
           zones.zone(k).chi(end,1)=1;
       end
   end
   if((zones.zone(k).pC.coord(3)==-1)&&(zones.zone(k).pC.coord(4)==1))
       if(~inpolygon(X_points.R(zones.zone(k).pC.coord(2)),...
               X_points.Z(zones.zone(k).pC.coord(2)),Rwall,Zwall))
           zones.zone(k).chi(end,end)=1;
       end
   end
   if((zones.zone(k).pD.coord(3)==-1)&&(zones.zone(k).pD.coord(4)==1))
       if(~inpolygon(X_points.R(zones.zone(k).pD.coord(2)),...
               X_points.Z(zones.zone(k).pD.coord(2)),Rwall,Zwall))
           zones.zone(k).chi(1,end)=1;
       end
   end
    
    chi=zones.zone(k).chi;
    [m,p]=size(chi);
    for i=1:m
        l=0;
        waitingforchi1=1;
        for j=1:p
            if((chi(i,j)==1)&&(waitingforchi1==1))
                if(l==1)
                    chi(i,j-1)=1;
                end
                if(l==2)
                    chi(i,j-1)=1;
                    chi(i,j-2)=1;
                end
                waitingforchi1=0;
                l=0;
            end
            if(chi(i,j)==0)
                l=l+1;
                waitingforchi1=1;
            end
        end
        if((chi(i,p)==0)&&(chi(i,p-1)==1))
            chi(i,p)=1;
        end
    end
    
    
    for i=1:m
        for j=2:p-1
            if((chi(i,j-1)==0)&&(chi(i,j)==1)&&(chi(i,j+1)==0))
               chi(i,j)=0; 
            end
        end
        for j=2:p-2
            if((chi(i,j-1)==0)&&(chi(i,j)==1)&&(chi(i,j+1)==1)&&(chi(i,j+2)==0))
               chi(i,j)=0; 
               chi(i,j+1)=0; 
            end
        end
    end
    
    zones.zone(k).chi=chi;
end

for k=1:zones.num
    chi=zones.zone(k).chi;
    chi_=zeros(size(chi)+[2 2]);
    chi_(2:end-1,2:end-1)=chi;
    North=zones.zone(k).Neighbour.north;
    South=zones.zone(k).Neighbour.south;
    East=zones.zone(k).Neighbour.east;
    West=zones.zone(k).Neighbour.west;
    if(West>0)
        chi_(2:end-1,1)=zones.zone(West).chi(:,end);
    else
        chi_(2:end-1,1)=chi_(2:end-1,2);
    end
    if(East>0)
        chi_(2:end-1,end)=zones.zone(East).chi(:,1);
    else
        chi_(2:end-1,end)=chi_(2:end-1,end-1);
    end
    if(North>0)
        chi_(end,2:end-1)=zones.zone(North).chi(1,:);
    else
        chi_(end,2:end-1)=chi_(end-1,2:end-1);
    end
    if(South>0)
        chi_(1,2:end-1)=zones.zone(South).chi(end,:);
    else
        chi_(1,2:end-1)=chi_(2,2:end-1);
    end
    
    [m,p]=size(chi);
    for i=1:m
        for j=1:p
            if(chi(i,j)==0)
                if(sum(sum(chi_(i:i+2,j:j+2)))>=7)
                    chi(i,j)=1;
                else
                    if(sum(sum(chi_(i:i+2,j:j+2)))>=6)
                       if((chi_(i,j+1)==0)&&(chi_(i+2,j+1)==0))
                           chi(i,j)=1;
                       end
                       if((chi_(i+1,j)==0)&&(chi_(i+1,j+2)==0))
                           chi(i,j)=1;
                       end
                    end
                end
            end
        end
    end
    
    zones.zone(k).chi=chi;
   
end
