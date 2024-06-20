global eirene;
ntriangles_=eirene.ntriangles_;
nknots_=eirene.nknots_;
R_=eirene.R_;
Z_=eirene.Z_;
knots_interp=eirene.knots_interp;
triangle=eirene.triangle;
triangles_=eirene.triangles_;

ntri_wall=eirene.wall.num;

hdf5write('triangles.h5','/Ntriangles',int64(ntriangles_));
hdf5write('triangles.h5','/Nknots',int64(nknots_),'WriteMode','append');

hdf5write('triangles.h5','/knots/R',100*R_,'WriteMode','append');
hdf5write('triangles.h5','/knots/Z',100*Z_,'WriteMode','append');

pass=zeros(1,nknots_);
neir=zeros(1,nknots_);
nsol=zeros(1,nknots_);
eir=zeros(nknots_,3);
sol=zeros(nknots_,8,3);

for n=1:nknots_
    pass(n)=knots_interp(n).pass;
    nsol(n)=knots_interp(n).nsol;
    neir(n)=knots_interp(n).neir;
    for n2=1:knots_interp(n).nsol
       sol(n,n2,1)=knots_interp(n).sol(n2,2); 
       sol(n,n2,2)=knots_interp(n).sol(n2,3); 
       sol(n,n2,3)=knots_interp(n).sol(n2,1); 
    end
    for n2=1:knots_interp(n).neir
       eir(n,n2)=knots_interp(n).eir(n2); 
    end
end

hdf5write('triangles.h5','/knots/pass',int64(pass),'WriteMode','append')
hdf5write('triangles.h5','/knots/neir',int64(neir),'WriteMode','append')
hdf5write('triangles.h5','/knots/nsol',int64(nsol),'WriteMode','append')
hdf5write('triangles.h5','/knots/sol',int64(sol),'WriteMode','append')
hdf5write('triangles.h5','/knots/eir',int64(eir),'WriteMode','append')

hdf5write('triangles.h5','/wall/Ntriangles',int64(eirene.wall.num),'WriteMode','append')

N_triwall=eirene.wall.num;
back_interp=zeros(N_triwall,4);
for k=1:N_triwall
    back_interp(k,1)=eirene.wall.tri(k);
    back_interp(k,2)=eirene.wall.tri_k(k);
    back_interp(k,3)=eirene.wall.tri_i(k);
    back_interp(k,4)=eirene.wall.tri_j(k);
end
hdf5write('triangles.h5','/wall/back_interp',int64(back_interp),'WriteMode','append');

North=zeros(ntri_wall,3);
for n=1:ntri_wall
    [m,p]=size(triangle(n).list_n);
    if(p>=1)
        North(n,1)=triangle(n).list_n(1,1);
        North(n,2)=triangle(n).list_n(2,1);
        North(n,3)=triangle(n).list_n(3,1);
    end
end
hdf5write('triangles.h5','/wall/north_s2d_to_use',int64(North),'WriteMode','append');
North=zeros(ntri_wall,1);
for n=1:ntri_wall
    [m,p]=size(triangle(n).list_n);
    if(p>=1)
        North(n,1)=triangle(n).weight_n(1,1);
    end
end
hdf5write('triangles.h5','/wall/north_s2d_weight',North,'WriteMode','append');

South=zeros(ntri_wall,3);
for n=1:ntri_wall
    [m,p]=size(triangle(n).list_s);
    if(p>=1)
        South(n,1)=triangle(n).list_s(1,1);
        South(n,2)=triangle(n).list_s(2,1);
        South(n,3)=triangle(n).list_s(3,1);
    end
end
hdf5write('triangles.h5','/wall/south_s2d_to_use',int64(South),'WriteMode','append');
South=zeros(ntri_wall,1);
for n=1:ntri_wall
    [m,p]=size(triangle(n).list_s);
    if(p>=1)
        South(n,1)=triangle(n).weight_s(1,1);
    end
end
hdf5write('triangles.h5','/wall/south_s2d_weight',South,'WriteMode','append');

East=zeros(ntri_wall,3);
for n=1:ntri_wall
    [m,p]=size(triangle(n).list_e);
    if(p>=1)
        East(n,1)=triangle(n).list_e(1,1);
        East(n,2)=triangle(n).list_e(2,1);
        East(n,3)=triangle(n).list_e(3,1);
    end
end
hdf5write('triangles.h5','/wall/east_s2d_to_use',int64(East),'WriteMode','append');
East=zeros(ntri_wall,1);
for n=1:ntri_wall
    [m,p]=size(triangle(n).list_e);
    if(p>=1)
        East(n,1)=triangle(n).weight_e(1,1);
    end
end
hdf5write('triangles.h5','/wall/east_s2d_weight',East,'WriteMode','append');

West=zeros(ntri_wall,3);
for n=1:ntri_wall
    [m,p]=size(triangle(n).list_w);
    if(p>=1)
        West(n,1)=triangle(n).list_w(1,1);
        West(n,2)=triangle(n).list_w(2,1);
        West(n,3)=triangle(n).list_w(3,1);
    end
end
hdf5write('triangles.h5','/wall/west_s2d_to_use',int64(West),'WriteMode','append');
West=zeros(ntri_wall,1);
for n=1:ntri_wall
    [m,p]=size(triangle(n).list_w);
    if(p>=1)
        West(n,1)=triangle(n).weight_w(1,1);
    end
end
hdf5write('triangles.h5','/wall/west_s2d_weight',West,'WriteMode','append');

tri_knots=zeros(ntriangles_,3);
for k=1:ntriangles_
    tri_knots(k,1)=triangles_(k).p1;
    tri_knots(k,2)=triangles_(k).p2;
    tri_knots(k,3)=triangles_(k).p3;
end
hdf5write('triangles.h5','/triangles/tri_knots',int64(tri_knots),'WriteMode','append')

type_face=zeros(ntriangles_,3);
for k=1:ntriangles_
    type_face(k,1)=triangles_(k).BC1;
    type_face(k,2)=triangles_(k).BC2;
    type_face(k,3)=triangles_(k).BC3;
end
hdf5write('triangles.h5','/triangles/type_face',int64(type_face),'WriteMode','append')

H=zeros(ntriangles_,3);
for k=1:ntriangles_
    H(k,1)=triangles_(k).k;
    H(k,2)=triangles_(k).i;
    H(k,3)=triangles_(k).j;
end
hdf5write('triangles.h5','/triangles/back_interp',int64(H),'WriteMode','append')
