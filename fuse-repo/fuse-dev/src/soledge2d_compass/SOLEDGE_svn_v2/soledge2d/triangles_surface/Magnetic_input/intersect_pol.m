function [x,y]=intersect_pol(P1x,P1y,P2x,P2y,polx,poly)

tol=1e-4;
err=1;
while(err>tol)
    a=inpolygon(P1x,P1y,polx,poly);
    b=inpolygon(P2x,P2y,polx,poly);
    P3x=(P1x+P2x)/2;
    P3y=(P1y+P2y)/2;
    c=inpolygon(P3x,P3y,polx,poly);
    if(a~=c)
        P2x=P3x;
        P2y=P3y;
    else
        P1x=P3x;
        P1y=P3y;
    end
    err=sqrt((P1x-P2x)^2+(P1y-P2y)^2);
end

x=(P1x+P2x)/2;
y=(P1y+P2y)/2;
return

end