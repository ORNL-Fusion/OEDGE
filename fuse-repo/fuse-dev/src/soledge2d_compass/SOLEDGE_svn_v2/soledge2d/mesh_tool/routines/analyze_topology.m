global megazone;
global Pmegazone;

remaining=[1:zones.num];

megazone.num=0;

while(length(remaining>0))
    
    megazone.num=megazone.num+1;
    ref=remaining(1);
    megazone.mz(megazone.num).list=[ref];
    
    %look east
    cur=ref;
    num=0;
    type=zones.zone(cur).pC.coord(4);
    while((type~=2)&&(num~=ref))
        east=[zones.zone(cur).pD.coord;...
            zones.zone(cur).pC.coord];
        
        for k=1:length(remaining)
            num=remaining(k);
            west=[zones.zone(num).pA.coord;...
                zones.zone(num).pB.coord];
            somme=0;
            for i=1:2
                for j=1:4
                    if(west(i,j)==east(i,j))
                        somme=somme+1;
                    else
                        if((j==3)||(j==7))
                            if((west(j)==-1)||(east(j)==-1))
                                somme=somme+1;
                                %joker
                            end
                        end
                    end
                end
            end
            if(somme==8)
                megazone.mz(megazone.num).list=[megazone.mz(megazone.num).list,num];
                cur=num;
                type=zones.zone(num).pC.coord(4);
                break;
            end
            
            if(num~=cur)
                west=[zones.zone(num).pD.coord;...
                    zones.zone(num).pC.coord];
                somme=0;
                for i=1:2
                    for j=1:4
                        if(west(i,j)==east(i,j))
                            somme=somme+1;
                        else
                            if((j==3)||(j==7))
                                if((west(j)==-1)||(east(j)==-1))
                                    somme=somme+1;
                                    %joker
                                end
                            end
                        end
                    end
                end
                if(somme==8)
                    megazone.mz(megazone.num).list=[megazone.mz(megazone.num).list,num];
                    cur=num;
                    type=zones.zone(num).pB.coord(4);
                    %return zone
                    temp.R=zones.zone(num).east.R;
                    temp.Z=zones.zone(num).east.Z;
                    zones.zone(num).east.R=zones.zone(num).west.R;
                    zones.zone(num).east.Z=zones.zone(num).west.Z;
                    zones.zone(num).west.R=temp.R;
                    zones.zone(num).west.Z=temp.Z;
                    zones.zone(num).north.R=zones.zone(num).north.R(end:-1:1);
                    zones.zone(num).north.Z=zones.zone(num).north.Z(end:-1:1);
                    zones.zone(num).south.R=zones.zone(num).south.R(end:-1:1);
                    zones.zone(num).south.Z=zones.zone(num).south.Z(end:-1:1);
                    pA=zones.zone(num).pA;
                    pB=zones.zone(num).pB;
                    zones.zone(num).pA=zones.zone(num).pD;
                    zones.zone(num).pB=zones.zone(num).pC;
                    zones.zone(num).pC=pB;
                    zones.zone(num).pD=pA;
                    break;
                end
            end
        end
    end
    
    %     if((megazone.mz(megazone.num).list(end)==megazone.mz(megazone.num).list(1))&&...
    %             (~zones.zone(megazone.mz(megazone.num).list(1)).pA.coord(4)==2))
    %         megazone.mz(megazone.num).isperiodic=true;
    %         megazone.mz(megazone.num).list=megazone.mz(megazone.num).list(1:end-1);
    %
    %     else
    %         megazone.mz(megazone.num).isperiodic=false;
    
    
    %look west
    cur=ref;
    num=0;
    type=zones.zone(cur).pB.coord(4);
    while((type~=2)&&(num~=ref))
        west=[zones.zone(cur).pA.coord;...
            zones.zone(cur).pB.coord];
        
        for k=1:length(remaining)
            num=remaining(k);
            east=[zones.zone(num).pD.coord;...
                zones.zone(num).pC.coord];
            somme=0;
            for i=1:2
                for j=1:4
                    if(west(i,j)==east(i,j))
                        somme=somme+1;
                    else
                        if((j==3)||(j==7))
                            if((west(j)==-1)||(east(j)==-1))
                                somme=somme+1;
                                %joker
                            end
                        end
                    end
                end
            end
            if(somme==8)
                megazone.mz(megazone.num).list=[num,megazone.mz(megazone.num).list];
                cur=num;
                type=zones.zone(num).pB.coord(4);
                break;
            end
            
            
            if(num~=cur)
                east=[zones.zone(num).pA.coord;...
                    zones.zone(num).pB.coord];
                somme=0;
                for i=1:2
                    for j=1:4
                        if(west(i,j)==east(i,j))
                            somme=somme+1;
                        else
                            if((j==3)||(j==7))
                                if((west(j)==-1)||(east(j)==-1))
                                    somme=somme+1;
                                    %joker
                                end
                            end
                        end
                    end
                end
                if(somme==8)
                    megazone.mz(megazone.num).list=[num,megazone.mz(megazone.num).list];
                    cur=num;
                    type=zones.zone(num).pC.coord(4);
                    %return zone
                    temp.R=zones.zone(num).east.R;
                    temp.Z=zones.zone(num).east.Z;
                    zones.zone(num).east.R=zones.zone(num).west.R;
                    zones.zone(num).east.Z=zones.zone(num).west.Z;
                    zones.zone(num).west.R=temp.R;
                    zones.zone(num).west.Z=temp.Z;
                    zones.zone(num).north.R=zones.zone(num).north.R(end:-1:1);
                    zones.zone(num).north.Z=zones.zone(num).north.Z(end:-1:1);
                    zones.zone(num).south.R=zones.zone(num).south.R(end:-1:1);
                    zones.zone(num).south.Z=zones.zone(num).south.Z(end:-1:1);
                    pA=zones.zone(num).pA;
                    pB=zones.zone(num).pB;
                    zones.zone(num).pA=zones.zone(num).pD;
                    zones.zone(num).pB=zones.zone(num).pC;
                    zones.zone(num).pC=pB;
                    zones.zone(num).pD=pA;
                    break;
                end
            end
            
        end
    end
    
    %     end
    
    %remove from remaining
    remaining_=remaining;
    remaining=[];
    for k=1:length(remaining_)
        if(length(find(megazone.mz(megazone.num).list==remaining_(k)))==0)
            remaining=[remaining,remaining_(k)];
        end
    end
    
end

for k=1:megazone.num
    %     if(megazone.mz(k).isperiodic)
    megazone.mz(k).isperiodic=false;
    ref=megazone.mz(k).list(1);
    list=[ref];
    for k2=2:length(megazone.mz(k).list)
        if(megazone.mz(k).list(k2)~=ref)
            list=[list,megazone.mz(k).list(k2)];
        else
            megazone.mz(k).isperiodic=true;
            break;
        end
    end
    megazone.mz(k).list=list;
    %     end
    for nz=1:length(megazone.mz(k).list)
        zones.zone(megazone.mz(k).list(nz)).mz=k;
    end
end

%megazone1 used as a reference for east west direction
done=[];
todo=[1];
while(length(done)<megazone.num)
    
    ref=todo(1);
    
    for k=1:length(megazone.mz(ref).list)
        
        zones.zone(megazone.mz(ref).list(k)).Neighbour.north=0;
        zones.zone(megazone.mz(ref).list(k)).Neighbour.south=0;
        
        if(megazone.mz(ref).list(k)==17)
           disp('coc'); 
        end
        if(megazone.mz(ref).list(k)==30)
           disp('coc'); 
        end
        
        %north scan
        north=[zones.zone(megazone.mz(ref).list(k)).pB.coord;...
            zones.zone(megazone.mz(ref).list(k)).pC.coord];
        
        for k2=1:zones.num
            south=[zones.zone(k2).pA.coord;...
                zones.zone(k2).pD.coord];
            
            somme=0;
            for i=1:2
                for j=1:4
                    if(north(i,j)==south(i,j))
                        somme=somme+1;
                    else
                        if((j==3)||(j==7))
                            if((west(j)==-1)||(east(j)==-1))
                                somme=somme+1;
                                %joker
                            end
                        end
                    end
                end
            end
            if(somme==8)
                %check sens
                sum=0;
                vx=zones.zone(megazone.mz(ref).list(k)).north.R(end)-...
                    zones.zone(megazone.mz(ref).list(k)).north.R(1);
                vy=zones.zone(megazone.mz(ref).list(k)).north.Z(end)-...
                    zones.zone(megazone.mz(ref).list(k)).north.Z(1);
                for k5=1:length(zones.zone(megazone.mz(ref).list(k)).north.R)
                    vx1=zones.zone(megazone.mz(ref).list(k)).north.R(k5)-...
                        zones.zone(megazone.mz(ref).list(k)).north.R(1);
                    vy1=zones.zone(megazone.mz(ref).list(k)).north.Z(k5)-...
                        zones.zone(megazone.mz(ref).list(k)).north.Z(1);
                    sum=sum+vx*vy1-vy*vx1;
                end
                sum1=0;
                vx=zones.zone(k2).south.R(end)-...
                    zones.zone(k2).south.R(1);
                vy=zones.zone(k2).south.Z(end)-...
                    zones.zone(k2).south.Z(1);
                for k5=1:length(zones.zone(k2).south.R)
                    vx1=zones.zone(k2).south.R(k5)-...
                        zones.zone(k2).south.R(1);
                    vy1=zones.zone(k2).south.Z(k5)-...
                        zones.zone(k2).south.Z(1);
                    sum1=sum1+vx*vy1-vy*vx1;
                end
                if(sum1*sum>=0)
                    zones.zone(megazone.mz(ref).list(k)).Neighbour.north=k2;
                    if(isempty(find(done==zones.zone(k2).mz)))
                        todo=[todo,zones.zone(k2).mz];
                    end
                    break;
                end
            end
            
            %reverse
            south=[zones.zone(k2).pD.coord;...
                zones.zone(k2).pA.coord];
            
            somme=0;
            for i=1:2
                for j=1:4
                    if(north(i,j)==south(i,j))
                        somme=somme+1;
                    else
                        if((j==3)||(j==7))
                            if((west(j)==-1)||(east(j)==-1))
                                somme=somme+1;
                                %joker
                            end
                        end
                    end
                end
            end
            if(somme==8)
                %check sens
                sum=0;
                vx=zones.zone(megazone.mz(ref).list(k)).north.R(end)-...
                    zones.zone(megazone.mz(ref).list(k)).north.R(1);
                vy=zones.zone(megazone.mz(ref).list(k)).north.Z(end)-...
                    zones.zone(megazone.mz(ref).list(k)).north.Z(1);
                for k5=1:length(zones.zone(megazone.mz(ref).list(k)).north.R)
                    vx1=zones.zone(megazone.mz(ref).list(k)).north.R(k5)-...
                        zones.zone(megazone.mz(ref).list(k)).north.R(1);
                    vy1=zones.zone(megazone.mz(ref).list(k)).north.Z(k5)-...
                        zones.zone(megazone.mz(ref).list(k)).north.Z(1);
                    sum=sum+vx*vy1-vy*vx1;
                end
                sum1=0;
                vx=zones.zone(k2).south.R(1)-...
                    zones.zone(k2).south.R(end);
                vy=zones.zone(k2).south.Z(1)-...
                    zones.zone(k2).south.Z(end);
                for k5=1:length(zones.zone(k2).south.R)
                    vx1=zones.zone(k2).south.R(k5)-...
                        zones.zone(k2).south.R(end);
                    vy1=zones.zone(k2).south.Z(k5)-...
                        zones.zone(k2).south.Z(end);
                    sum1=sum1+vx*vy1-vy*vx1;
                end
                if(sum1*sum>=0)
                    zones.zone(megazone.mz(ref).list(k)).Neighbour.north=k2;
                    reverse_megazone(zones.zone(k2).mz);
                    if(isempty(find(done==zones.zone(k2).mz)))
                        todo=[todo,zones.zone(k2).mz];
                    end
                    break;
                end
            end
        end
        
        
        %south scan
        south=[zones.zone(megazone.mz(ref).list(k)).pA.coord;...
            zones.zone(megazone.mz(ref).list(k)).pD.coord];
        
        for k2=1:zones.num
            north=[zones.zone(k2).pB.coord;...
                zones.zone(k2).pC.coord];
            
            somme=0;
            for i=1:2
                for j=1:4
                    if(north(i,j)==south(i,j))
                        somme=somme+1;
                    else
                        if((j==3)||(j==7))
                            if((west(j)==-1)||(east(j)==-1))
                                somme=somme+1;
                                %joker
                            end
                        end
                    end
                end
            end
            if(somme==8)
                %check sens
                sum=0;
                vx=zones.zone(megazone.mz(ref).list(k)).south.R(end)-...
                    zones.zone(megazone.mz(ref).list(k)).south.R(1);
                vy=zones.zone(megazone.mz(ref).list(k)).south.Z(end)-...
                    zones.zone(megazone.mz(ref).list(k)).south.Z(1);
                for k5=1:length(zones.zone(megazone.mz(ref).list(k)).south.R)
                    vx1=zones.zone(megazone.mz(ref).list(k)).south.R(k5)-...
                        zones.zone(megazone.mz(ref).list(k)).south.R(1);
                    vy1=zones.zone(megazone.mz(ref).list(k)).south.Z(k5)-...
                        zones.zone(megazone.mz(ref).list(k)).south.Z(1);
                    sum=sum+vx*vy1-vy*vx1;
                end
                sum1=0;
                vx=zones.zone(k2).north.R(end)-...
                    zones.zone(k2).north.R(1);
                vy=zones.zone(k2).north.Z(end)-...
                    zones.zone(k2).north.Z(1);
                for k5=1:length(zones.zone(k2).north.R)
                    vx1=zones.zone(k2).north.R(k5)-...
                        zones.zone(k2).north.R(1);
                    vy1=zones.zone(k2).north.Z(k5)-...
                        zones.zone(k2).north.Z(1);
                    sum1=sum1+vx*vy1-vy*vx1;
                end
                if(sum1*sum>=0)
                    zones.zone(megazone.mz(ref).list(k)).Neighbour.south=k2;
                    if(isempty(find(done==zones.zone(k2).mz)))
                        todo=[todo,zones.zone(k2).mz];
                    end
                    break;
                end
            end
            
            %reverse
            north=[zones.zone(k2).pC.coord;...
                zones.zone(k2).pB.coord];
            
            somme=0;
            for i=1:2
                for j=1:4
                    if(north(i,j)==south(i,j))
                        somme=somme+1;
                    else
                        if((j==3)||(j==7))
                            if((west(j)==-1)||(east(j)==-1))
                                somme=somme+1;
                                %joker
                            end
                        end
                    end
                end
            end
            if(somme==8)
                %check sens
                sum=0;
                vx=zones.zone(megazone.mz(ref).list(k)).south.R(end)-...
                    zones.zone(megazone.mz(ref).list(k)).south.R(1);
                vy=zones.zone(megazone.mz(ref).list(k)).south.Z(end)-...
                    zones.zone(megazone.mz(ref).list(k)).south.Z(1);
                for k5=1:length(zones.zone(megazone.mz(ref).list(k)).south.R)
                    vx1=zones.zone(megazone.mz(ref).list(k)).south.R(k5)-...
                        zones.zone(megazone.mz(ref).list(k)).south.R(1);
                    vy1=zones.zone(megazone.mz(ref).list(k)).south.Z(k5)-...
                        zones.zone(megazone.mz(ref).list(k)).south.Z(1);
                    sum=sum+vx*vy1-vy*vx1;
                end
                sum1=0;
                vx=zones.zone(k2).north.R(1)-...
                    zones.zone(k2).north.R(end);
                vy=zones.zone(k2).north.Z(1)-...
                    zones.zone(k2).north.Z(end);
                for k5=1:length(zones.zone(k2).north.R)
                    vx1=zones.zone(k2).north.R(k5)-...
                        zones.zone(k2).north.R(end);
                    vy1=zones.zone(k2).north.Z(k5)-...
                        zones.zone(k2).north.Z(end);
                    sum1=sum1+vx*vy1-vy*vx1;
                end
                if(sum1*sum>=0)
                    zones.zone(megazone.mz(ref).list(k)).Neighbour.south=k2;
                    reverse_megazone(zones.zone(k2).mz);
                    if(isempty(find(done==zones.zone(k2).mz)))
                        todo=[todo,zones.zone(k2).mz];
                    end
                    break;
                end
            end
        end
        
    end
    
    done=[done,ref];
    todo_=[];
    for k=1:length(todo)
        if(isempty(find(done==todo(k))))
            todo_=[todo_,todo(k)];
        end
    end
    todo=todo_;
end

for k=1:megazone.num
    for nz=1:length(megazone.mz(k).list)-1
        zones.zone(megazone.mz(k).list(nz)).Neighbour.east=megazone.mz(k).list(nz+1);
        zones.zone(megazone.mz(k).list(nz+1)).Neighbour.west=megazone.mz(k).list(nz);
    end
    if(megazone.mz(k).isperiodic)
        zones.zone(megazone.mz(k).list(1)).Neighbour.west=megazone.mz(k).list(end);
        zones.zone(megazone.mz(k).list(end)).Neighbour.east=megazone.mz(k).list(1);
    else
        zones.zone(megazone.mz(k).list(1)).Neighbour.west=0;
        zones.zone(megazone.mz(k).list(end)).Neighbour.east=0;
    end
end


%weak check
for k=1:zones.num
    if(zones.zone(k).Neighbour.north>0)
        N=zones.zone(k).Neighbour.north;
        if(zones.zone(N).Neighbour.south~=k)
            Errordlg(['Problem with neighbors: Zone ',num2str(k)]);
        end
    end
    if(zones.zone(k).Neighbour.south>0)
        S=zones.zone(k).Neighbour.south;
        if(zones.zone(S).Neighbour.north~=k)
            Errordlg(['Problem with neighbors: Zone ',num2str(k)]);
        end
    end
    if(zones.zone(k).Neighbour.east>0)
        E=zones.zone(k).Neighbour.east;
        if(zones.zone(E).Neighbour.west~=k)
            Errordlg(['Problem with neighbors: Zone ',num2str(k)]);
        end
    end
    if(zones.zone(k).Neighbour.west>0)
        W=zones.zone(k).Neighbour.west;
        if(zones.zone(W).Neighbour.east~=k)
            Errordlg(['Problem with neighbors: Zone ',num2str(k)]);
        end
    end
end

for k=1:megazone.num
    megazone.mz(k).ismeshed=false;
end


Pmegazone.num=0
for k=1:zones.num
    if(zones.zone(k).Neighbour.south==0)
        Pmegazone.num=Pmegazone.num+1;
        Pmegazone.mz(Pmegazone.num).list=[k];
        zones.zone(k).pmz=Pmegazone.num;
        num=k;
        while(zones.zone(num).Neighbour.north~=0)
            num=zones.zone(num).Neighbour.north;
            zones.zone(num).pmz=Pmegazone.num;
            Pmegazone.mz(Pmegazone.num).list=[Pmegazone.mz(Pmegazone.num).list,num];
        end
    end
end

for k=1:Pmegazone.num
    Pmegazone.mz(k).isaligned=false;
    Pmegazone.mz(k).subrefpoints(1).R=[];
    Pmegazone.mz(k).subrefpoints(1).Z=[];
    Pmegazone.mz(k).subrefpoints(2).R=[];
    Pmegazone.mz(k).subrefpoints(2).Z=[];
end
for k=1:zones.num
   zones.zone(k).northaligned=false;
   zones.zone(k).southaligned=false;
end