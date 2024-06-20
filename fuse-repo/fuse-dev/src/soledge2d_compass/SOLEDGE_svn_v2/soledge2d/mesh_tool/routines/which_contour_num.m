function n = which_contour_num(cin,p1)

    d=zeros(1,cin.num);
    for k=1:cin.num
        dist=sqrt((cin.arc(k).x-p1.x).^2+(cin.arc(k).y-p1.y).^2);
        d(k)=min(dist);
    end
    
    n=find(d==min(d));

end