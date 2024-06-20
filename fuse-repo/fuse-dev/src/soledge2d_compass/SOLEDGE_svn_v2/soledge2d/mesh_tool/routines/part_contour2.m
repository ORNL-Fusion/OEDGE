function C = part_contour2(Cin,p1,p2)

C1=part_contour(Cin,p1);

d=sqrt((Cin.x-p2.x).^2+(Cin.y-p2.y).^2);
ind=find(d==min(d));
if(ind>C1.ind)
    Cin2.x=C1.arc(2).x;
    Cin2.y=C1.arc(2).y;
    C2=part_contour(Cin2,p2);
    C.x=C2.arc(1).x;
    C.y=C2.arc(1).y;
else
    if(ind==C1.ind)
        C.x=[p1.x, p2.x];
        C.y=[p1.y, p2.y];
    else
        Cin2.x=C1.arc(1).x;
        Cin2.y=C1.arc(1).y;
        C2=part_contour(Cin2,p2);
        C.x=C2.arc(2).x;
        C.y=C2.arc(2).y;
    end
end

end

