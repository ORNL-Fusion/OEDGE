global megazone
global zones;

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
      dist=dist/dist(end);
      distnew=dist;
      for j=2:length(distnew)-1
         dd=distnew(j)-distnew(j-1);
         dd2=distnew(j+1)-distnew(j-1);
         rat=dd/dd2;
         if(rat<0.4)
             ratn=0.4;
         else
             if(rat>0.6)
                 ratn=0.6;
             else
                 ratn=rat;
             end
         end
         distnew(j)=distnew(j-1)+ratn*dd2;
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