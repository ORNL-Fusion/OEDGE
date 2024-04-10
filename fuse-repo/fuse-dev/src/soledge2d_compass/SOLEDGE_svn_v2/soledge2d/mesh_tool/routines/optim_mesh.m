function [dout,Rout,Zout]=optim_mesh(din,Rin,Zin,cx,cy,sens)

global dmin_ad;
global range_ad;
global length_ad;
global alpha_ad;
global orthoptim;
global target_ad;

% close all
% return optimized mesh on piece of contour c
% generated from near points R1, Z1


Npts=length(Rin);
dout = din;

dist=zeros(size(cx));
long=length(cx);
for k=2:long
    dist(k)=dist(k-1)+sqrt((cx(k)-cx(k-1))^2+(cy(k)-cy(k-1))^2);
end
delta_min=dmin_ad*1e-3; 
delta_min=delta_min/dist(end);
dist=dist/dist(end);
%first iteration
Rout=interp1(dist,cx,dout);
Zout=interp1(dist,cy,dout);

% figure(1)
% clf
% hold on
% plot([Rin;Rout],[Zin;Zout],'b-')
% plot([Rout],[Zout],'ro')
% plot([Rin;Rout]',[Zin;Zout]','b-')

long=length(Rin);
cost=zeros(1,long-1);
for k=2:long-1
    vx=Rin(k+1)-Rin(k);
    vy=Zin(k+1)-Zin(k);
    vx1=Rout(k)-Rin(k);
    vy1=Zout(k)-Zin(k);
    qual=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2));%+...
        %(vx2*vy-vy2*vx)/(sqrt(vx^2+vy^2)*sqrt(vx2^2+vy2^2));
%     cost(k)=qual*sens;
    l1=dout(k+1)-dout(k);
    l2=dout(k)-dout(k-1);
    cost(k)=qual*sens*orthoptim+((l1/l2+l2/l1)-1)*(1-orthoptim);
end

% figure(2)
% plot(cost,'b.-')
% hold on
% plot(mean(cost)*ones(size(cost)),'g-')
% plot(min(cost)*ones(size(cost)),'r-')

cost_er=[min(cost)];

n=0;
while((cost_er(end)<target_ad)&&(n<length_ad))
    
    %compute grad
    grad=zeros(1,long-2);
    for k=1:long-2
        dd_=0.01*(dout(k+2)-dout(k+1));
        range=range_ad;
        a=min(range,(dout(k+1)-dout(1)));
        dd(1:k+1)=max((dd_+dd_*(dout(1:k+1)-dout(k+1))/a),0);
        a=min(range,(dout(end)-dout(k+1)));
        dd(k+2:length(dout))=max((dd_-dd_*(dout(k+2:end)-dout(k+1))/a),0);
%         figure(8)
%         plot(dd,'r.-')
        dout_t=dout+dd;
        dout_t(1)=dout(1);
        dout_t(end)=dout(end);
%         dout_t(1:k+1)=dout(1:k+1)*(dout(k+1)+0.1*(dout(k+2)-dout(k+1)))/dout(k+1);
%         dout_t(k+2:length(dout))=dout(end)-((dout(end)-dout(k+2:end))*((dout(end)-(dout(k+1)+0.1*(dout(k+2)-dout(k+1))))...
%             /(dout(end)-dout(k+1))));
        dout_tt(k,:)=dout_t;
        Rout_t=interp1(dist,cx,dout_t);
        Zout_t=interp1(dist,cy,dout_t);
        cost_t=zeros(1,long-1);
        for k2=2:long-1
            vx=Rin(k2+1)-Rin(k2);
            vy=Zin(k2+1)-Zin(k2);
            vx1=Rout_t(k2)-Rin(k2);
            vy1=Zout_t(k2)-Zin(k2);
            qual=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2));%+...
                %(vx2*vy-vy2*vx)/(sqrt(vx^2+vy^2)*sqrt(vx2^2+vy2^2));
            l1=dout_t(k2+1)-dout_t(k2);
            l2=dout_t(k2)-dout_t(k2-1);
            cost_t(k2)=qual*sens*orthoptim+((l1/l2+l2/l1)-1)*(1-orthoptim);
        end
        grad(k)=(min(cost_t(2:end))-min(cost(2:end)))/0.01;
    end

    
    %     figure(3)
    %     plot(grad,'k.-')
    
    alpha=alpha_ad;
    eps=1e-6;
%     ddout=zeros(size(dout));
    dout_new=dout;
    for k=1:long-2
        dout_new=dout_new+grad(k)./sqrt(sum(grad.^2)+eps)*(dout_tt(k,:)-dout)*alpha;
        %         s=s-grad(k)*alpha;
    end
%     for k=2:length(ddout)-1
%         ddout(k)=max(ddout(k),0.4*(dout(k+1)-dout(k)));
%         ddout(k)=min(ddout(k),0.4*(dout(k-1)-dout(k)));
%     end
%     dout_new=dout+ddout;
    delt=dout_new(2:end)-dout_new(1:end-1);
    for k=1:length(delt)
        if(delt(k)<delta_min)
            delt(k)=delta_min;
        end
        dout_new(k+1)=dout_new(k)+delt(k);
    end
    dout_new=dout_new/dout_new(end);
    %     dout_new=dout_new/s;
    R_out_new=interp1(dist,cx,dout_new);
    Z_out_new=interp1(dist,cy,dout_new);
    
    cost_new=zeros(1,long-1);
    for k=2:long-1
        vx=Rin(k+1)-Rin(k);
        vy=Zin(k+1)-Zin(k);
        vx1=R_out_new(k)-Rin(k);
        vy1=Z_out_new(k)-Zin(k);
        qual=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2));%+...
            %(vx2*vy-vy2*vx)/(sqrt(vx^2+vy^2)*sqrt(vx2^2+vy2^2));
        l1=dout_new(k+1)-dout_new(k);
        l2=dout_new(k)-dout_new(k-1);
        cost_new(k)=qual*sens*orthoptim+((l1/l2+l2/l1)-1)*(1-orthoptim);
    end
    
%         figure(4)
%         clf
%         hold on
%         plot([Rin;R_out_new],[Zin;Z_out_new],'b-')
%         plot([R_out_new],[Z_out_new],'ro')
%         plot([Rin;R_out_new]',[Zin;Z_out_new]','b-')
%         figure(5)
%         plot(cost_new,'b.-')
%         hold on
%         plot(mean(cost_new)*ones(size(cost_new)),'g-')
%         plot(min(cost_new)*ones(size(cost_new)),'r-')
%     
    
    dout=dout_new;
    Rout=R_out_new;
    Zout=Z_out_new;
    cost=cost_new;
    
    cost_er=[cost_er,min(cost(2:end))];
    n=n+1;
%         figure(6)
%         plot(cost_er,'k.-')
end

% figure(4)
% clf
% hold on
% plot([Rin;Rout],[Zin;Zout],'b-')
% plot([Rout],[Zout],'ro')
% plot([Rin;Rout]',[Zin;Zout]','b-')
% figure(5)
% plot(cost,'b.-')
% hold on
% plot(mean(cost)*ones(size(cost)),'g-')
% plot(min(cost)*ones(size(cost)),'r-')
% 
