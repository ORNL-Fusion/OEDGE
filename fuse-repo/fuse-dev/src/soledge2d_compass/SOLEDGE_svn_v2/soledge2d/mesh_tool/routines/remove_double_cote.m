clear cotes
ncotes=0;
list=[];

tol=0.01; %1cm

for nc=1:ncote
    
    [m,p]=size(list);
    in=false;
    for k=1:m
        if((sum(list(k,:)==[cote(nc).p1.coord([1,2,4]),cote(nc).p2.coord([1,2,4])])==6)||(sum(list(k,:)==[cote(nc).p2.coord([1,2,4]),cote(nc).p1.coord([1,2,4])])==6))
            d=zeros(1,length(cote(nc).R));
            for k2=1:length(cote(nc).R)
                d(k2)=min(sqrt((cotes(k).R-cote(nc).R(k2)).^2+(cotes(k).Z-cote(nc).Z(k2)).^2));
            end
            dm=sum(d==0)
            nc
            if(length(d)==2)
                if(dm==2)
                    in=true
                end
            else
                if(dm>2) %more than 2 points overlap
                    in=true;    
                end
            end
        end
    end
    
    d=0;
    for k=1:length(cote(nc).R)-1
        d=d+sqrt((cote(nc).R(k+1)-cote(nc).R(k))^2+(cote(nc).Z(k+1)-cote(nc).Z(k))^2);
    end
    if(d==0)
        zerol=true;
    else
        zerol=false;
    end
    
    if((~in)&&(~zerol))
        list=[list;cote(nc).p1.coord([1,2,4]),cote(nc).p2.coord([1,2,4])];
        ncotes=ncotes+1;
        cotes(ncotes).p1=cote(nc).p1;
        cotes(ncotes).p2=cote(nc).p2;
        cotes(ncotes).R=cote(nc).R;
        cotes(ncotes).Z=cote(nc).Z;
        cotes(ncotes).cote=nc;
    end
    
end


list=[];
for nc=1:ncotes
    if(((cotes(nc).p1.coord(4)==1)&&(cotes(nc).p1.coord(3)==-1))&&((cotes(nc).p2.coord(4)==1)&&(cotes(nc).p2.coord(3)==-1))) % boucle
        list=[list,nc];
    end
end

dist=zeros(1,length(list));
if(length(list)>2)
    disp('should look at remove_double_cote.m - might be a problem')
end
% if(length(list)==1)
%     for k=1:length(list)
%         dist(k)=sum(sqrt((cotes(list(k)).R(2:end)-cotes(list(k)).R(1:end-1)).^2+(cotes(list(k)).Z(2:end)-cotes(list(k)).Z(1:end-1)).^2));
%     end
%     [a,b]=min(dist);
%     %remove b(1)
%     ncr=list(b(1));
%     for k=ncr+1:ncotes
%        cotes(k-1).R=cotes(k).R; 
%        cotes(k-1).Z=cotes(k).Z; 
%        cotes(k-1).p1=cotes(k).p1; 
%        cotes(k-1).p2=cotes(k).p2; 
%     end
%     ncotes=ncotes-1;
% end
