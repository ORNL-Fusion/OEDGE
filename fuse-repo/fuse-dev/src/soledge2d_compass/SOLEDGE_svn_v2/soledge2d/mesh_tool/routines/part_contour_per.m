function cout = part_contour_per(cin,p1,p2,dir)

d=sqrt((cin.x-p1.x).^2+(cin.y-p1.y).^2);
a=find(d==min(d));

if(dir==1)
    if(a<length(d))
        v1x=cin.x(a+1)-cin.x(a);
        v1y=cin.y(a+1)-cin.y(a);
        vx=p1.x-cin.x(a);
        vy=p1.y-cin.y(a);
        if(v1x*vx+v1y*vy>0)
            cleft.x=circshift(cin.x(1:end-1)',-a)';
            cleft.y=circshift(cin.y(1:end-1)',-a)';
        else
            cleft.x=circshift(cin.x(1:end-1)',-a+1)';
            cleft.y=circshift(cin.y(1:end-1)',-a+1)';
        end
        cleft.x=[cleft.x,cleft.x(1)];
        cleft.y=[cleft.y,cleft.y(1)];
    else
        cleft.x=cin.x;
        cleft.y=cin.y;
    end
    
    
    %check if p2 in [p1, cleft(1)]
    d=sqrt((cleft.x-p2.x).^2+(cleft.y-p2.y).^2);
    a2=find(d==min(d));
    if((a2(1)==length(d))||(a2(1)==length(d)-1)||(a2(1)==1))
        vx=p2.x-p1.x;
        vy=p2.y-p2.x;
        v1x=cleft.x(end)-cleft.x(end-1);
        v1y=cleft.y(end)-cleft.y(end-1);
        if(vx*v1x+vy*v1y>0)
            cout.x=[];
            cout.y=[];
        else
            cout.x=cleft.x(1:end-1)
            cout.y=cleft.y(1:end-1)
        end
    else
        v1x=cleft.x(a2+1)-cleft.x(a2);
        v1y=cleft.y(a2+1)-cleft.y(a2);
        vx=p2.x-cleft.x(a2);
        vy=p2.y-cleft.y(a2);
        if(v1x*vx+v1y*vy<0)
            cout.x=cleft.x(1:a2-1);
            cout.y=cleft.y(1:a2-1);
        else
            cout.x=cleft.x(1:a2);
            cout.y=cleft.y(1:a2);
        end
    end
    
    cout.x=[p1.x,cout.x,p2.x];
    cout.y=[p1.y,cout.y,p2.y];
else
    
    if(a<length(d))
        v1x=cin.x(a+1)-cin.x(a);
        v1y=cin.y(a+1)-cin.y(a);
        vx=p1.x-cin.x(a);
        vy=p1.y-cin.y(a);
        if(v1x*vx+v1y*vy>0)
            cleft.x=circshift(cin.x(1:end-1)',-a)';
            cleft.y=circshift(cin.y(1:end-1)',-a)';
        else
            cleft.x=circshift(cin.x(1:end-1)',-a+1)';
            cleft.y=circshift(cin.y(1:end-1)',-a+1)';
        end
        cleft.x=[cleft.x,cleft.x(1)];
        cleft.y=[cleft.y,cleft.y(1)];
    else
        cleft.x=cin.x;
        cleft.y=cin.y;
    end
    cright.x=cleft.x(end-1:-1:1);
    cright.y=cleft.y(end-1:-1:1);
    cright.x=[cright.x,cright.x(1)];
    cright.y=[cright.y,cright.y(1)];
    
    d=sqrt((cright.x-p2.x).^2+(cright.y-p2.y).^2);
    a2=find(d==min(d));
    if((a2(1)==length(d))||(a2(1)==length(d)-1)||(a2(1)==1))
        vx=p2.x-p1.x;
        vy=p2.y-p1.y;
        v1x=cright.x(end)-cright.x(end-1);
        v1y=cright.y(end)-cright.y(end-1);
        if(vx*v1x+vy*v1y>0)
            cout.x=[];
            cout.y=[];
        else
            cout.x=cright.x(1:end-1)
            cout.y=cright.y(1:end-1)
        end
    else
        v1x=cright.x(a2+1)-cright.x(a2);
        v1y=cright.y(a2+1)-cright.y(a2);
        vx=p2.x-cright.x(a2);
        vy=p2.y-cright.y(a2);
        if(v1x*vx+v1y*vy<0)
            cout.x=cright.x(1:a2-1);
            cout.y=cright.y(1:a2-1);
        else
            cout.x=cright.x(1:a2);
            cout.y=cright.y(1:a2);
        end
    end
    
    cout.x=[p1.x,cout.x,p2.x];
    cout.y=[p1.y,cout.y,p2.y];
    
end



end