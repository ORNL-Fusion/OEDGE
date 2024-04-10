global eirene

ntri_wall=eirene.wall.num;
tri_left=[2:ntri_wall];

tri_seq=[1];

R=zeros(ntriangles_,3);
Z=zeros(ntriangles_,3);
for n=1:ntriangles_
    R(n,1)=R_(triangles_(n).p1)*100;
    R(n,2)=R_(triangles_(n).p2)*100;
    R(n,3)=R_(triangles_(n).p3)*100;
    Z(n,1)=Z_(triangles_(n).p1)*100;
    Z(n,2)=Z_(triangles_(n).p2)*100;
    Z(n,3)=Z_(triangles_(n).p3)*100;
end

err=1e-8;
while(length(tri_left)>0)
    ntri_end=triangle(tri_seq(end)).ntri;
    if(triangle(tri_seq(end)).side==1)
        p2R=R(ntri_end,2);
        p2Z=Z(ntri_end,2);
    end
    if(triangle(tri_seq(end)).side==2)
        p2R=R(ntri_end,3);
        p2Z=Z(ntri_end,3);
    end
    if(triangle(tri_seq(end)).side==3)
        p2R=R(ntri_end,1);
        p2Z=Z(ntri_end,1);
    end
    for k=1:length(tri_left)
        ntri_=triangle(tri_left(k)).ntri;
        if(triangle(tri_left(k)).side==1)
            p1R=R(ntri_,1);
            p1Z=Z(ntri_,1);
        end
        if(triangle(tri_left(k)).side==2)
            p1R=R(ntri_,2);
            p1Z=Z(ntri_,2);
        end
        if(triangle(tri_left(k)).side==3)
            p1R=R(ntri_,3);
            p1Z=Z(ntri_,3);
        end
        d=sqrt((p2R-p1R)^2+(p2Z-p1Z)^2);
        if(d<err)
            break
        end
    end
    tri_seq=[tri_seq,tri_left(k)];
    tri_left=[tri_left([1:k-1]),tri_left([k+1:end])];
end

eirene.wall.sequence=tri_seq;
eirene.wall.typemat=ones(size(eirene.wall.sequence));
eirene.wall.ispump=false(size(eirene.wall.sequence));