global eirene

for n=1:ntriangles
    %segment1
    [ntri1,np1]=find(tri_knots==triangles(n).p1);
    [ntri2,np2]=find(tri_knots==triangles(n).p2);
    ntric=intersect(ntri1,ntri2);
    ntrineigh1=setdiff(ntric,n);
    if(isempty(ntrineigh1))
        if(triangles(n).area==1)
            triangles(n).neigh1=0;
            triangles(n).typeneigh1=0;
            triangles(n).BC1=1; % core BC       
        else
            triangles(n).neigh1=0;
            triangles(n).typeneigh1=0;
            triangles(n).BC1=3; % wall BC
        end
    else
        triangles(n).BC1=0;
        triangles(n).neigh1=ntrineigh1;
        seg=[np1(find(ntri1==ntrineigh1)) np2(find(ntri2==ntrineigh1))];
        if((sum(seg==[1 2])==2)||(sum(seg==[2 1])==2))
            triangles(n).typeneigh1=1;
        else
            if((sum(seg==[2 3])==2)||(sum(seg==[3 2])==2))
                triangles(n).typeneigh1=2;
            else
                triangles(n).typeneigh1=3;
            end
        end
    end
    %segment2
    [ntri1,np1]=find(tri_knots==triangles(n).p2);
    [ntri2,np2]=find(tri_knots==triangles(n).p3);
    ntric=intersect(ntri1,ntri2);
    ntrineigh1=setdiff(ntric,n);
    if(isempty(ntrineigh1))
        if(triangles(n).area==1)
            triangles(n).neigh2=0;
            triangles(n).typeneigh2=0;
            triangles(n).BC2=1; % core BC
        else
            triangles(n).neigh2=0;
            triangles(n).typeneigh2=0;
            triangles(n).BC2=3; % wall BC
        end
    else
        triangles(n).BC2=0;
        triangles(n).neigh2=ntrineigh1;
        seg=[np1(find(ntri1==ntrineigh1)) np2(find(ntri2==ntrineigh1))];
        if((sum(seg==[1 2])==2)||(sum(seg==[2 1])==2))
            triangles(n).typeneigh2=1;
        else
            if((sum(seg==[2 3])==2)||(sum(seg==[3 2])==2))
                triangles(n).typeneigh2=2;
            else
                triangles(n).typeneigh2=3;
            end
        end
    end
    %segment3
    [ntri1,np1]=find(tri_knots==triangles(n).p3);
    [ntri2,np2]=find(tri_knots==triangles(n).p1);
    ntric=intersect(ntri1,ntri2);
    ntrineigh1=setdiff(ntric,n);
    if(isempty(ntrineigh1))
        if(triangles(n).area==1)
            triangles(n).neigh3=0;
            triangles(n).typeneigh3=0;
            triangles(n).BC3=1; % core BC
        else
            triangles(n).neigh3=0;
            triangles(n).typeneigh3=0;
            triangles(n).BC3=3; % wall BC
        end
    else
        triangles(n).BC3=0;
        triangles(n).neigh3=ntrineigh1;
        seg=[np1(find(ntri1==ntrineigh1)) np2(find(ntri2==ntrineigh1))];
        if((sum(seg==[1 2])==2)||(sum(seg==[2 1])==2))
            triangles(n).typeneigh3=1;
        else
            if((sum(seg==[2 3])==2)||(sum(seg==[3 2])==2))
                triangles(n).typeneigh3=2;
            else
                triangles(n).typeneigh3=3;
            end
        end
    end
end