global zones
global eirene

ntrirem=0;

ntrinew=1;
ntrinew_table=zeros(1,ntriangles);
remknot=zeros(1,eirene.nknots);
nknotnew_table=zeros(1,eirene.nknots);
for n=1:ntriangles
    if(triangles(n).step==0)
        ntrirem=ntrirem+1;
        %remove knots only in this triangle
        nknot1=triangles(n).p1;
        [m,p]=find(tri_knots==nknot1);
        remknot(nknot1)=1;
        for n2=1:length(m)
            if(triangles(m(n2)).step==1)
                remknot(nknot1)=0;
            end
        end
        nknot2=triangles(n).p2;
        [m,p]=find(tri_knots==nknot2);
        remknot(nknot2)=1;
        for n2=1:length(m)
            if(triangles(m(n2)).step==1)
                remknot(nknot2)=0;
            end
        end
        nknot3=triangles(n).p3;
        [m,p]=find(tri_knots==nknot3);
        remknot(nknot3)=1;
        for n2=1:length(m)
            if(triangles(m(n2)).step==1)
                remknot(nknot3)=0;
            end
        end
    else
        ntrinew_table(n)=ntrinew;
        ntrinew=ntrinew+1;
    end
end

disp(['triangles removed ',num2str(ntrirem)])

ntriangles_=ntrinew-1;

nknotnew=1;
for n=1:eirene.nknots
    if(remknot(n)==0)
        nknotnew_table(n)=nknotnew;
        nknotnew=nknotnew+1;
    end
end

nknots_=nknotnew-1;

for n=1:ntriangles
    if(ntrinew_table(n)~=0)
        triangles_(ntrinew_table(n)).k=triangles(n).k;
        triangles_(ntrinew_table(n)).i=triangles(n).i;
        triangles_(ntrinew_table(n)).j=triangles(n).j;
        triangles_(ntrinew_table(n)).p1=nknotnew_table(triangles(n).p1);
        triangles_(ntrinew_table(n)).p2=nknotnew_table(triangles(n).p2);
        triangles_(ntrinew_table(n)).p3=nknotnew_table(triangles(n).p3);
        triangles_(ntrinew_table(n)).area=triangles(n).area;
        triangles_(ntrinew_table(n)).BC1=triangles(n).BC1;
        triangles_(ntrinew_table(n)).BC2=triangles(n).BC2;
        triangles_(ntrinew_table(n)).BC3=triangles(n).BC3;
        if(triangles(n).neigh1>0)
            triangles_(ntrinew_table(n)).neigh1=ntrinew_table(triangles(n).neigh1);
        else
            triangles_(ntrinew_table(n)).neigh1=triangles(n).neigh1;
        end
        if(triangles(n).neigh2>0)
            triangles_(ntrinew_table(n)).neigh2=ntrinew_table(triangles(n).neigh2);
        else
            triangles_(ntrinew_table(n)).neigh2=triangles(n).neigh2;
        end
        if(triangles(n).neigh3>0)
            triangles_(ntrinew_table(n)).neigh3=ntrinew_table(triangles(n).neigh3);
        else
            triangles_(ntrinew_table(n)).neigh3=triangles(n).neigh3;
        end
        triangles_(ntrinew_table(n)).typeneigh1=triangles(n).typeneigh1;
        triangles_(ntrinew_table(n)).typeneigh2=triangles(n).typeneigh2;
        triangles_(ntrinew_table(n)).typeneigh3=triangles(n).typeneigh3;
        triangles_(ntrinew_table(n)).step=triangles(n).step;
    end
end

R_=zeros(1,nknots_);
Z_=zeros(1,nknots_);

for n=1:eirene.nknots
    if(nknotnew_table(n)~=0)
        R_(nknotnew_table(n))=eirene.R(n);
        Z_(nknotnew_table(n))=eirene.Z(n);
    end
end

for k=1:nzones
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(zones.zone(k).knotA(i,j)~=0)
                zones.zone(k).knotA_(i,j)=nknotnew_table(zones.zone(k).knotA(i,j));
            else
                zones.zone(k).knotA_(i,j)=0;
            end
        end
    end
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(zones.zone(k).knotB(i,j)~=0)
                zones.zone(k).knotB_(i,j)=nknotnew_table(zones.zone(k).knotB(i,j));
            else
                zones.zone(k).knotB_(i,j)=0;
            end
        end
    end
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(zones.zone(k).knotC(i,j)~=0)
                zones.zone(k).knotC_(i,j)=nknotnew_table(zones.zone(k).knotC(i,j));
            else
                zones.zone(k).knotC_(i,j)=0;
            end
        end
    end
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(zones.zone(k).knotD(i,j)~=0)
                zones.zone(k).knotD_(i,j)=nknotnew_table(zones.zone(k).knotD(i,j));
            else
                zones.zone(k).knotD_(i,j)=0;
            end
        end
    end
    
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(zones.zone(k).knotE(i,j)~=0)
                zones.zone(k).knotE_(i,j)=nknotnew_table(zones.zone(k).knotE(i,j));
            else
                zones.zone(k).knotE_(i,j)=0;
            end
        end
    end
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(zones.zone(k).knotF(i,j)~=0)
                zones.zone(k).knotF_(i,j)=nknotnew_table(zones.zone(k).knotF(i,j));
            else
                zones.zone(k).knotF_(i,j)=0;
            end
        end
    end
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(zones.zone(k).knotG(i,j)~=0)
                zones.zone(k).knotG_(i,j)=nknotnew_table(zones.zone(k).knotG(i,j));
            else
                zones.zone(k).knotG_(i,j)=0;
            end
        end
    end
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(zones.zone(k).knotH(i,j)~=0)
                zones.zone(k).knotH_(i,j)=nknotnew_table(zones.zone(k).knotH(i,j));
            else
                zones.zone(k).knotH_(i,j)=0;
            end
        end
    end
end

tri_knots_=zeros(ntriangles_,3);
for n=1:ntriangles_
    tri_knots_(n,1:3)=[triangles_(n).p1,triangles_(n).p2,triangles_(n).p3];
end
