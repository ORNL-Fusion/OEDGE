for k=1:X_points.num
    for n=1:4
       c1.num=1;
       c1.arc(1).x=X_points.cut(k).arc(n).R;
       c1.arc(1).y=X_points.cut(k).arc(n).Z;
       c2=contour_better(r2D,z2D,flux2D,X_points.cut(k).psilim(n));
       X=intersect_contour(c1,c2);
       X_points.cut(k).arc(n).psiR=c2.arc(X.arc2).x;
       X_points.cut(k).arc(n).psiZ=c2.arc(X.arc2).y;
    end
end