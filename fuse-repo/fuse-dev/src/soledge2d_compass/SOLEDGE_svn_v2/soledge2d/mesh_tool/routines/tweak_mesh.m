global Pmegazone;
global megazone;
global dmin_ad;
global drmin_ad;
global zones;

%radial tweak
for k=1:Pmegazone.num
    R=[];
    Z=[];
   for k1=1:length(Pmegazone.mz(k).list)
       nz=Pmegazone.mz(k).list(k1);
       R=[R;zones.zone(nz).gridR(1:end-1,:)];
       Z=[Z;zones.zone(nz).gridZ(1:end-1,:)];
   end
   R=[R;zones.zone(nz).gridR(end,:)];
   Z=[Z;zones.zone(nz).gridZ(end,:)];
   [m,p]=size(R);
   for j=1:p
      R_=R(:,j);
      Z_=Z(:,j);
      dist=zeros(size(R_));
      for i=2:length(dist)
         dist(i)=dist(i-1)+sqrt((R_(i)-R_(i-1))^2+(Z_(i)-Z_(i-1))^2);
      end
      dmin=drmin_ad*1.e-3/dist(end);
      dist=dist/dist(end);
      distnew=dist;
      for i=2:length(distnew)
         dd=distnew(i)-distnew(i-1);
         if(dd<dmin)
             distnew(i)=distnew(i-1)+dmin;
         end
      end
      distnew=distnew/distnew(end);
      
      R_new=interp1(dist,R_,distnew);
      Z_new=interp1(dist,Z_,distnew);
      R(:,j)=R_new;
      Z(:,j)=Z_new;
   end
   for k1=1:length(Pmegazone.mz(k).list)
       nz=Pmegazone.mz(k).list(k1);
       [m,p]=size(zones.zone(nz).gridR);
       zones.zone(nz).gridR=R(1:m,1:p);
       zones.zone(nz).gridZ=Z(1:m,1:p);
       R=R(m:end,:);
       Z=Z(m:end,:);
   end
end

%poloidal tweak
for k=1:megazone.num
    R=[];
    Z=[];
   for k1=1:length(megazone.mz(k).list)
       nz=megazone.mz(k).list(k1);
       R=[R,zones.zone(nz).gridR(:,1:end-1)];
       Z=[Z,zones.zone(nz).gridZ(:,1:end-1)];
   end
   R=[R,zones.zone(nz).gridR(:,end)];
   Z=[Z,zones.zone(nz).gridZ(:,end)];
   [m,p]=size(R);
   for i=1:m
      R_=R(i,:);
      Z_=Z(i,:);
      dist=zeros(size(R_));
      for j=2:length(dist)
         dist(j)=dist(j-1)+sqrt((R_(j)-R_(j-1))^2+(Z_(j)-Z_(j-1))^2);
      end
      dmin=dmin_ad*1.e-3/dist(end);
      dist=dist/dist(end);
      distnew=dist;
      for j=2:length(distnew)
         dd=distnew(j)-distnew(j-1);
         if(dd<dmin)
             distnew(j)=distnew(j-1)+dmin;
         end
      end
      distnew=distnew/distnew(end);
      
      R_new=interp1(dist,R_,distnew);
      Z_new=interp1(dist,Z_,distnew);
      R(i,:)=R_new;
      Z(i,:)=Z_new;
   end
   for k1=1:length(megazone.mz(k).list)
       nz=megazone.mz(k).list(k1);
       [m,p]=size(zones.zone(nz).gridR);
       zones.zone(nz).gridR=R(1:m,1:p);
       zones.zone(nz).gridZ=Z(1:m,1:p);
       R=R(:,p:end);
       Z=Z(:,p:end);
   end
end

% %radial tweak
% for k=1:Pmegazone.num
%     for k1=1:length(Pmegazone.mz(k).list)
%         nz=Pmegazone.mz(k).list(k1);
%         if(~isempty(zones.zone(nz).gridR))
%             [m,p]=size(zones.zone(nz).gridR);
%             if(zones.zone(nz).Neighbour.south>0)
%                 if(~isempty(zones.zone(zones.zone(nz).Neighbour.south).gridR))
%                     for i=m:-1:1
%                         zones.zone(nz).gridR(i,:)=zones.zone(nz).gridR(i,:)-zones.zone(nz).gridR(1,:)+zones.zone(zones.zone(nz).Neighbour.south).gridR(end,:);
%                         zones.zone(nz).gridZ(i,:)=zones.zone(nz).gridZ(i,:)-zones.zone(nz).gridZ(1,:)+zones.zone(zones.zone(nz).Neighbour.south).gridZ(end,:);
%                     end
%                 end
%             end
%             
%             vx=zones.zone(nz).gridR(m-1,2)-zones.zone(nz).gridR(m-1,1);
%             vy=zones.zone(nz).gridZ(m-1,2)-zones.zone(nz).gridZ(m-1,1);
%             vx1=zones.zone(nz).gridR(m,1)-zones.zone(nz).gridR(m-1,1);
%             vy1=zones.zone(nz).gridZ(m,1)-zones.zone(nz).gridZ(m-1,1);
%             sens=sign(vx*vy1-vx1*vy);
%             for i=2:m
%                 for j=1:p-1
%                     vx=zones.zone(nz).gridR(i-1,j+1)-zones.zone(nz).gridR(i-1,j);
%                     vy=zones.zone(nz).gridZ(i-1,j+1)-zones.zone(nz).gridZ(i-1,j);
%                     vx1=zones.zone(nz).gridR(i,j)-zones.zone(nz).gridR(i-1,j);
%                     vy1=zones.zone(nz).gridZ(i,j)-zones.zone(nz).gridZ(i-1,j);
%                     vx1n=(vx*vx1+vy*vy1)/sqrt(vx^2+vy^2)*vx;
%                     vy1n=vy1-(vx*vx1+vy*vy1)/sqrt(vx^2+vy^2)*vy;
%                     vx1o=vx1-vx1n;
%                     vy1o=vy1-vy1n;
%                     dr=(vx*vy1-vx1*vy)*sens/sqrt((vx^2+vy^2));
%                     if(dr<drmin_ad*1e-3)
%                         for i2=m:-1:i
%                             zones.zone(nz).gridR(i2,j)=zones.zone(nz).gridR(i2,j)+vx1o/sqrt(vx1o^2+vy1o^2)*(drmin_ad*1.e-3-dr);
%                             zones.zone(nz).gridZ(i2,j)=zones.zone(nz).gridZ(i2,j)+vy1o/sqrt(vx1o^2+vy1o^2)*(drmin_ad*1.e-3-dr);
%                         end
%                     end
%                 end
%                 %last
%                 vx=zones.zone(nz).gridR(i-1,p)-zones.zone(nz).gridR(i-1,p-1);
%                 vy=zones.zone(nz).gridZ(i-1,p)-zones.zone(nz).gridZ(i-1,p-1);
%                 vx1=zones.zone(nz).gridR(i,p)-zones.zone(nz).gridR(i-1,p);
%                 vy1=zones.zone(nz).gridZ(i,p)-zones.zone(nz).gridZ(i-1,p);
%                 vx1n=(vx*vx1+vy*vy1)/sqrt(vx^2+vy^2)*vx;
%                 vy1n=vy1-(vx*vx1+vy*vy1)/sqrt(vx^2+vy^2)*vy;
%                 vx1o=vx1-vx1n;
%                 vy1o=vy1-vy1n;
%                 dr=(vx*vy1-vx1*vy)*sens/sqrt((vx^2+vy^2));
%                 if(dr<drmin_ad*1e-3)
%                     for i2=m:-1:i
%                         zones.zone(nz).gridR(i2,p)=zones.zone(nz).gridR(i2,p)+vx1o/sqrt(vx1o^2+vy1o^2)*(drmin_ad*1.e-3-dr);
%                             zones.zone(nz).gridZ(i2,p)=zones.zone(nz).gridZ(i2,p)+vy1o/sqrt(vx1o^2+vy1o^2)*(drmin_ad*1.e-3-dr);
% %                         zones.zone(nz).gridR(i2,p)=zones.zone(nz).gridR(i2,p)-zones.zone(nz).gridR(i,p)+zones.zone(nz).gridR(i-1,p)+vx1o*(drmin_ad*1.e-3/dr);
% %                         zones.zone(nz).gridZ(i2,p)=zones.zone(nz).gridZ(i2,p)-zones.zone(nz).gridZ(i,p)+zones.zone(nz).gridZ(i-1,p)+vy1o*(drmin_ad*1.e-3/dr);
%                     end
%                 end
%             end
%         end
%     end
% end