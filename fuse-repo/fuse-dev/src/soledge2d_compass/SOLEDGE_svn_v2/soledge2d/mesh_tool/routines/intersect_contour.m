function x=intersect_contour(c1,c2)
    
    x.num=0;
    for i=1:c1.num
        for j=1:c2.num
            [xb,yb]=polyxpoly(c1.arc(i).x,c1.arc(i).y,c2.arc(j).x,c2.arc(j).y,'unique');
            [m,p]=size(xb);
            for k=1:m
                x.num=x.num+1;
                x.x(x.num)=xb(k);
                x.y(x.num)=yb(k);
                x.arc1(x.num)=i;
                x.arc2(x.num)=j;
            end
        end
    end

end