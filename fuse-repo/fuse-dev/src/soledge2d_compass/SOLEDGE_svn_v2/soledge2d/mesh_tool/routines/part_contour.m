function C=part_contour(Cin,p1)

% input ordered contour
% output pieces of contour delimited by p1

d=sqrt((Cin.x-p1.x).^2+(Cin.y-p1.y).^2);
ind=find(d==min(d));
vx=p1.x-Cin.x(ind);
vy=p1.y-Cin.y(ind);

if(ind>1)
    vmx=p1.x-Cin.x(ind-1);
    vmy=p1.y-Cin.y(ind-1);
end
if(ind<length(d))
    vpx=p1.x-Cin.x(ind+1);
    vpy=p1.y-Cin.y(ind+1);
end

if(ind==1)
    if(vx*vpx+vy*vpy>0)
        C.arc(1).x=p1.x;
        C.arc(1).y=p1.y;
        C.arc(2).x=Cin.x;
        C.arc(2).y=Cin.y;
    else
        C.arc(1).x=[Cin.x(1),p1.x];
        C.arc(1).y=[Cin.y(1),p1.y];
        C.arc(2).x=[p1.x,Cin.x(2:end)];
        C.arc(2).y=[p1.y,Cin.y(2:end)];
    end
else
    if(ind==length(d))
        if(vx*vmx+vy*vmy>0)
            C.arc(2).x=p1.x;
            C.arc(2).y=p1.y;
            C.arc(1).x=Cin.x;
            C.arc(1).y=Cin.y;
        else
            C.arc(1).x=[Cin.x(1:end-1),p1.x];
            C.arc(1).y=[Cin.y(1:end-1),p1.y];
            C.arc(2).x=[p1.x,Cin.x(end)];
            C.arc(2).y=[p1.y,Cin.y(end)];
        end
    else %general case
        if(vx*vpx+vy*vpy>0) %p1 between ind and ind-1 probably
            C.arc(1).x=[Cin.x(1:ind-1),p1.x];
            C.arc(1).y=[Cin.y(1:ind-1),p1.y];
            C.arc(2).x=[p1.x,Cin.x(ind:end)];
            C.arc(2).y=[p1.y,Cin.y(ind:end)];
        else
            C.arc(1).x=[Cin.x(1:ind),p1.x];
            C.arc(1).y=[Cin.y(1:ind),p1.y];
            C.arc(2).x=[p1.x,Cin.x(ind+1:end)];
            C.arc(2).y=[p1.y,Cin.y(ind+1:end)];
        end
    end
end

C.ind=ind;

end