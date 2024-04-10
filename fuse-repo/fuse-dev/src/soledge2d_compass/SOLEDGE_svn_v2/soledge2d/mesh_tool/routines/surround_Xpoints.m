
r=2*(r2D(1,2)-r2D(1,1));
theta=linspace(0,2*pi,100);

for k=1:X_points.num
    r=X_points.R(k)/100;
    c1=contour_better(r2D,z2D,flux2D,X_points.psi(k));
    c2.num=1;
    c2.arc(1).x=X_points.R(k)+r*cos(theta);
    c2.arc(1).y=X_points.Z(k)+r*sin(theta);
    I=intersect_contour(c1,c2);
    arg=zeros(1,4);
    for n=1:4
       cplx=(I.x(n)-X_points.R(k))+complex(0,1)*(I.y(n)-X_points.Z(k));
       arg(n)=angle(cplx);
    end
    [a,b]=sort(arg);
    for n=1:4
        X_points.branch(k).R(n)=I.x(b(n));
        X_points.branch(k).Z(n)=I.y(b(n));
        X_points.branch(k).theta(n)=arg(b(n));
    end
    
    
    B=circshift([1:4]',-1);
    for n=1:4
        X_points.cut(k).theta(n)=(X_points.branch(k).theta(n)+...
            X_points.branch(k).theta(B(n)))/2;
        if(abs(X_points.cut(k).theta(n)-X_points.branch(k).theta(n))>pi/2)
            X_points.cut(k).theta(n)=X_points.cut(k).theta(n)+pi;
        end % choose the small angle
        X_points.cut(k).R(n)=X_points.R(k)+r*cos(X_points.cut(k).theta(n));
        X_points.cut(k).Z(n)=X_points.Z(k)+r*sin(X_points.cut(k).theta(n));
        X_points.cut(k).psi(n)=interp2(r2D,z2D,flux2D,X_points.cut(k).R(n),X_points.cut(k).Z(n));
    end
    
    if(X_points.cut(k).psi(1)>X_points.cut(k).psi(2)) %reorder
       R_=X_points.branch(k).R(4);
       Z_=X_points.branch(k).Z(4);
       theta_=X_points.branch(k).theta(4);
       for i=3:-1:1
           X_points.branch(k).R(i+1)=X_points.branch(k).R(i);
           X_points.branch(k).Z(i+1)=X_points.branch(k).Z(i);
           X_points.branch(k).theta(i+1)=X_points.branch(k).theta(i);
       end
       X_points.branch(k).R(1)=R_;
       X_points.branch(k).Z(1)=Z_;
       X_points.branch(k).theta(1)=theta_;
       R_=X_points.cut(k).R(4);
       Z_=X_points.cut(k).Z(4);
       theta_=X_points.cut(k).psi(4);
       for i=3:-1:1
           X_points.cut(k).R(i+1)=X_points.cut(k).R(i);
           X_points.cut(k).Z(i+1)=X_points.cut(k).Z(i);
           X_points.cut(k).psi(i+1)=X_points.cut(k).psi(i);
       end
       X_points.cut(k).R(1)=R_;
       X_points.cut(k).Z(1)=Z_;
       X_points.cut(k).psi(1)=theta_;
    end
    
    text(X_points.branch(k).R(1),X_points.branch(k).Z(1),'1')
    text(X_points.branch(k).R(2),X_points.branch(k).Z(2),'2')
    text(X_points.branch(k).R(3),X_points.branch(k).Z(3),'3')
    text(X_points.branch(k).R(4),X_points.branch(k).Z(4),'4')
    text(X_points.cut(k).R(1),X_points.cut(k).Z(1),'A')
    text(X_points.cut(k).R(2),X_points.cut(k).Z(2),'B')
    text(X_points.cut(k).R(3),X_points.cut(k).Z(3),'C')
    text(X_points.cut(k).R(4),X_points.cut(k).Z(4),'D')
end