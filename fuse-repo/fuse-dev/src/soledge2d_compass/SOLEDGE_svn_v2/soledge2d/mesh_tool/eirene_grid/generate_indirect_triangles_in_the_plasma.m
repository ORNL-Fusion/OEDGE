% global Rwall
% global Zwall
global zones
global eirene

ntriangles=0;

tri_knots=[];

%attributing triangle numbers and vertices
for k=1:nzones
    Nx=zones.zone(k).Nx;
    Nz=zones.zone(k).Nz;
    for i=1:Nx
        for j=1:Nz
            if(((sum(zones.zone(k).iscrossed(i,j,:))==0)&&(zones.zone(k).inplasma(i,j)==1))||((zones.zone(k).isaligned(i,j)==1)&&(inpolygon(zones.zone(k).R(i,j),zones.zone(k).Z(i,j),eirene.Rwall,eirene.Zwall))))
                %premier triangle
                nt1=zones.zone(k).knotA(i,j);
                nt2=zones.zone(k).knotB(i,j);
                nt3=zones.zone(k).knotC(i,j);
                ntriangles=ntriangles+1;
                
                triangles(ntriangles).k=k;
                triangles(ntriangles).i=i;
                triangles(ntriangles).j=j;
                triangles(ntriangles).p1=nt1;
                triangles(ntriangles).p2=nt3;
                triangles(ntriangles).p3=nt2;
                
                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                
                %deuxieme triangle
                nt1=zones.zone(k).knotA(i,j);
                nt2=zones.zone(k).knotC(i,j);
                nt3=zones.zone(k).knotD(i,j);
                ntriangles=ntriangles+1;
                
                triangles(ntriangles).k=k;
                triangles(ntriangles).i=i;
                triangles(ntriangles).j=j;
                triangles(ntriangles).p1=nt1;
                triangles(ntriangles).p2=nt3;
                triangles(ntriangles).p3=nt2;
                
                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
            end
        end
    end
end

eirene.triangles=triangles;
eirene.tri_knots=tri_knots;
eirene.ntriangles=ntriangles;