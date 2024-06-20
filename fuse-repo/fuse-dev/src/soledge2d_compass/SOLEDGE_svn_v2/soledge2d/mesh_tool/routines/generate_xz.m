global zones;
global megazone;
global Pmegazone;
global r2D;
global z2D;
global flux2D;

for k=1:zones.num
    zones.zone(k).isx_defined=false;
    zones.zone(k).isz_defined=false;
end

minx=1e10;
maxx=-1e10;
for k=1:zones.num
    if(~zones.zone(k).isx_defined)
        R=zones.zone(k).gridR(:,1);
        Z=zones.zone(k).gridZ(:,1);
        psi=interp2(r2D,z2D,flux2D,R,Z);
        minx_=min(psi);
        maxx_=max(psi);
        if(minx>minx_)
            minx=minx_;
        end
        if(maxx<maxx_)
            maxx=maxx_;
        end
        mz=zones.zone(k).mz;
        for k1=1:length(megazone.mz(mz).list)
            nz=megazone.mz(mz).list(k1);
            [m,p]=size(zones.zone(nz).gridR);
            zones.zone(nz).xb=psi*ones(1,p);
%             zones.zone(nz).xb=[1:m]*ones(1,p);
            zones.zone(nz).x=(zones.zone(nz).xb(1:end-1,1:end-1)+zones.zone(nz).xb(2:end,1:end-1))*0.5;
            zones.zone(nz).isx_defined=true;
        end
    end
end

for k=1:zones.num
    zones.zone(k).x=(zones.zone(k).x-minx)/(maxx-minx);
    zones.zone(k).xb=(zones.zone(k).xb-minx)/(maxx-minx);
end


mzin=[];

while(~(length(mzin)==megazone.num))
    
    l=zeros(1,megazone.num-length(mzin));
    l2=zeros(1,megazone.num-length(mzin));
    k1=1;
    for k=1:megazone.num
        if(sum(find(k==mzin))==0)
            l(k1)=length(megazone.mz(k).list);
            l2(k1)=k;
            k1=k1+1;
        end
    end
    d=find(l==max(l));
    mz=l2(d(1));
    mzin=[mzin,mz];
    
    containsz=false;
    for k=1:length(megazone.mz(mz).list)
        containsz=containsz||zones.zone(megazone.mz(mz).list(k)).isz_defined;
    end
    
    if(~containsz)
        dist1=0;
        dist2=0;
        for k=1:length(megazone.mz(mz).list)
            nz=megazone.mz(mz).list(k);
            [m,p]=size(zones.zone(nz).gridZ);
            for j=1:p-1
                %middle point m/2
                dist1=[dist1,dist1(end)+sqrt((zones.zone(nz).gridR(floor(m/2),j+1)-zones.zone(nz).gridR(floor(m/2),j))^2+...
                    (zones.zone(nz).gridZ(floor(m/2),j+1)-zones.zone(nz).gridZ(floor(m/2),j))^2)];
%                 dist2=[dist2,dist2(end)+sqrt((zones.zone(nz).gridR(end,j+1)-zones.zone(nz).gridR(end,j))^2+...
%                     (zones.zone(nz).gridZ(end,j+1)-zones.zone(nz).gridZ(end,j))^2)];
            end
        end
%         if(dist1(end)>dist2(end))
%             dist=dist1;
%         else
%             dist=dist2;
%         end
        dist=dist1/dist1(end);
        
        %broadcast perp
        for k=1:length(megazone.mz(mz).list)
            nz=megazone.mz(mz).list(k);
            pmz=zones.zone(nz).pmz;
            for k1=1:length(Pmegazone.mz(pmz).list)
                nz2=Pmegazone.mz(pmz).list(k1);
                [m,p]=size(zones.zone(nz2).gridZ);
                if(~zones.zone(nz2).isz_defined)
                    zones.zone(nz2).zb=ones(m,1)*dist(1:p);
                    zones.zone(nz2).isz_defined=true;
                end
            end
            dist=dist(p:end);
        end
        
    else
         
        %begin
        tomesh=[];
        for k1=1:length(megazone.mz(mz).list)
            nz2=megazone.mz(mz).list(k1);
            if(~zones.zone(nz2).isz_defined)
                tomesh=[tomesh,nz2];
            else
                beg=nz2;
                break;
            end
        end
        dist1=0;
        dist2=0;
        for k2=1:length(tomesh)
            nz=tomesh(k2);
            [m,p]=size(zones.zone(nz).gridZ);
            for j=1:p-1
                dist1=[dist1,dist1(end)+sqrt((zones.zone(nz).gridR(floor(m/2),j+1)-zones.zone(nz).gridR(floor(m/2),j))^2+...
                    (zones.zone(nz).gridZ(floor(m/2),j+1)-zones.zone(nz).gridZ(floor(m/2),j))^2)];
%                 dist2=[dist2,dist2(end)+sqrt((zones.zone(nz).gridR(end,j+1)-zones.zone(nz).gridR(end,j))^2+...
%                     (zones.zone(nz).gridZ(end,j+1)-zones.zone(nz).gridZ(end,j))^2)];
            end
        end
%         if(dist1(end)>dist2(end))
%             dist=dist1;
%         else
%             dist=dist2;
%         end
        dist=dist1/dist1(end);
%         dist=dist*zones.zone(beg).zb(1,1);
        if(length(dist)>=2) 
            %rescale
            dist=dist/(dist(end)-dist(end-1))*(zones.zone(beg).zb(1,2)-zones.zone(beg).zb(1,1));
            %shift
            dist=dist-dist(end)+zones.zone(beg).zb(1,1);
        else
            dist=dist*zones.zone(beg).zb(1,1);
        end
        
        %broadcast perp
        for k2=1:length(tomesh)
            nz=tomesh(k2);
            pmz=zones.zone(nz).pmz;
            for k1=1:length(Pmegazone.mz(pmz).list)
                nz2=Pmegazone.mz(pmz).list(k1);
                [m,p]=size(zones.zone(nz2).gridZ);
                if(~zones.zone(nz2).isz_defined)
                    zones.zone(nz2).zb=ones(m,1)*dist(1:p);
                    zones.zone(nz2).isz_defined=true;
                end
            end
            dist=dist(p:end);
        end
      
        %end
        tomesh=[];
        for k1=length(megazone.mz(mz).list):-1:1
            nz2=megazone.mz(mz).list(k1);
            if(~zones.zone(nz2).isz_defined)
                tomesh=[tomesh,nz2];
            else
                beg=nz2;
                break;
            end
        end
        dist1=0;
        dist2=0;
        tomesh=tomesh(end:-1:1);
        for k2=1:length(tomesh)
            nz=tomesh(k2);
            [m,p]=size(zones.zone(nz).gridZ);
            for j=1:p-1
                dist1=[dist1,dist1(end)+sqrt((zones.zone(nz).gridR(floor(m/2),j+1)-zones.zone(nz).gridR(floor(m/2),j))^2+...
                    (zones.zone(nz).gridZ(floor(m/2),j+1)-zones.zone(nz).gridZ(floor(m/2),j))^2)];
%                 dist2=[dist2,dist2(end)+sqrt((zones.zone(nz).gridR(end,j+1)-zones.zone(nz).gridR(end,j))^2+...
%                     (zones.zone(nz).gridZ(end,j+1)-zones.zone(nz).gridZ(end,j))^2)];
            end
        end
%         if(dist1(end)>dist2(end))
%             dist=dist1;
%         else
%             dist=dist2;
%         end
        dist=dist1/dist1(end);
%         dist=zones.zone(beg).zb(1,end)+dist*(1-zones.zone(beg).zb(1,end));
        if(length(dist)>=2)
            %rescale
            dist=dist/(dist(2)-dist(1))*(zones.zone(beg).zb(1,end)-zones.zone(beg).zb(1,end-1));
            %shift
            dist=dist+zones.zone(beg).zb(1,end);
        else
            dist=zones.zone(beg).zb(1,end)+dist*(1-zones.zone(beg).zb(1,end));
        end
        %broadcast perp
        for k2=1:length(tomesh)
            nz=tomesh(k2);
            pmz=zones.zone(nz).pmz;
            for k1=1:length(Pmegazone.mz(pmz).list)
                nz2=Pmegazone.mz(pmz).list(k1);
                [m,p]=size(zones.zone(nz2).gridZ);
                if(~zones.zone(nz2).isz_defined)
                    zones.zone(nz2).zb=ones(m,1)*dist(1:p);
                    zones.zone(nz2).isz_defined=true;
                end
            end
            dist=dist(p:end);
        end     
    end
    
end

for k=1:zones.num
    zones.zone(k).z=(zones.zone(k).zb(1:end-1,1:end-1)+zones.zone(k).zb(1:end-1,2:end))*0.5;
end
