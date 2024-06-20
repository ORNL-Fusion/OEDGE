function c=contour_better(x,y,f,val,color)

if(nargin<5)
    color='b';
end

% c_raw=contour(x(1,:),y(:,1),f,[val val],color);
c_raw=contourc(x(1,:),y(:,1),f,[val val]);

[m,p]=size(c_raw);

if(p==0)
    c.num=0;
    return
end

k=1;
c.num=0;
while(k<p)
    k_=k;
    k=k+c_raw(2,k)+1;
    c.num=c.num+1;
    c.arc(c.num).x=c_raw(1,k_+1:k-1);
    c.arc(c.num).y=c_raw(2,k_+1:k-1);
%     %determining trigo
%     cx=mean(c.arc(c.num).x);
%     cy=mean(c.arc(c.num).y);
%     A=0;
%     for k=1:length(c.arc(c.num).x)-1
%         vx=c.arc(c.num).x(k)-cx;
%         vy=c.arc(c.num).y(k)-cy;
%         vx1=c.arc(c.num).x(k+1)-cx;
%         vy1=c.arc(c.num).y(k+1)-cy;
%         A=A+(vx*vy1-vx1*vy);
%     end
%     if(A<0)
%         c.arc(c.num).x=c.arc(c.num).x(end:-1:1);
%         c.arc(c.num).y=c.arc(c.num).y(end:-1:1);
%     end
end

end