cd ../Mesh/
for k=1:nzones
    name=strcat('chi_',num2str(k,'%.3d'));
    name=strcat(name,'.txt');
    zone(k).chi=load(name);
end
cd ../Magnetic_input/

ntriangles=0;


for k=1:nzones
    for i=1:zone(k).Nx
        for j=1:zone(k).Nz
%             if((zone(k).chi(i,j)==1)&&(zone(k).iscrossed(i,j)==0))
%                 zone(k).knotA(i,j)=0;
%                 zone(k).knotB(i,j)=0;
%                 zone(k).knotC(i,j)=0;
%                 zone(k).knotD(i,j)=0;
%                 zone(k).knotE(i,j)=0;
%                 zone(k).knotF(i,j)=0;
%                 zone(k).knotG(i,j)=0;
%                 zone(k).knotH(i,j)=0;
%             end
            if(zone(k).iscrossed(i,j)==1)
                PCr=zone(k).R(i,j);
                PCz=zone(k).Z(i,j);
                %cas 1 et 2
                if(zone(k).knotE(i,j)~=0)
                    %cas 1
                    if(zone(k).knotH(i,j)~=0)
                        triR=[R(zone(k).knotE(i,j));...
                            R(zone(k).knotH(i,j));...
                            R(zone(k).knotA(i,j))];
                        triZ=[Z(zone(k).knotE(i,j));...
                            Z(zone(k).knotH(i,j));...
                            Z(zone(k).knotA(i,j))];
                        %cas 1a
                        if(inpolygon(PCr,PCz,triR,triZ))
                            zone(k).knotB(i,j)=0;
                            zone(k).knotF(i,j)=0;
                            zone(k).knotC(i,j)=0;
                            zone(k).knotG(i,j)=0;
                            zone(k).knotD(i,j)=0;
                            % on crée 1 triangle
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotA(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotE(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotH(i,j);
                        else %cas 1b
                            zone(k).knotA(i,j)=0;
                            % on crée 3 triangles
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotE(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotB(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotH(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotB(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotC(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotH(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotH(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotC(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotD(i,j);
                        end
                    end
                    %cas 2
                    if(zone(k).knotF(i,j)~=0)
                        triR=[R(zone(k).knotB(i,j));...
                            R(zone(k).knotE(i,j));...
                            R(zone(k).knotF(i,j))];
                        triZ=[Z(zone(k).knotB(i,j));...
                            Z(zone(k).knotE(i,j));...
                            Z(zone(k).knotF(i,j))];
                        %cas 2a
                        if(inpolygon(PCr,PCz,triR,triZ))
                            zone(k).knotC(i,j)=0;
                            zone(k).knotG(i,j)=0;
                            zone(k).knotD(i,j)=0;
                            zone(k).knotH(i,j)=0;
                            zone(k).knotA(i,j)=0;
                            % on crée 1 triangle
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotB(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotF(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotE(i,j);
                        else %cas 2b
                            zone(k).knotB(i,j)=0;
                            % on crée 3 triangles
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotE(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotF(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotA(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotA(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotF(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotD(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotF(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotC(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotD(i,j);
                        end
                    end
                    %cas 3
                    if(zone(k).knotG(i,j)~=0)
                        triR=[R(zone(k).knotB(i,j));...
                            R(zone(k).knotC(i,j));...
                            R(zone(k).knotG(i,j));...
                            R(zone(k).knotE(i,j))];
                        triZ=[Z(zone(k).knotB(i,j));...
                            Z(zone(k).knotC(i,j));...
                            Z(zone(k).knotG(i,j));...
                            Z(zone(k).knotE(i,j))];
                        %cas 3a
                        if(inpolygon(PCr,PCz,triR,triZ))
                            zone(k).knotA(i,j)=0;
                            zone(k).knotH(i,j)=0;
                            zone(k).knotD(i,j)=0;
                            % on crée 2 triangles
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotB(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotG(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotE(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotB(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotC(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotG(i,j);
                        else %cas 3b
                            zone(k).knotB(i,j)=0;
                            zone(k).knotF(i,j)=0;
                            zone(k).knotC(i,j)=0;
                            % on crée 2 triangles
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotE(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotD(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotA(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotE(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotG(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotD(i,j);
                        end
                    end
                end
                if(zone(k).knotG(i,j)~=0)
                    %cas4
                    if(zone(k).knotF(i,j)~=0)
                        triR=[R(zone(k).knotF(i,j));...
                            R(zone(k).knotC(i,j));...
                            R(zone(k).knotG(i,j))];
                        triZ=[Z(zone(k).knotF(i,j));...
                            Z(zone(k).knotC(i,j));...
                            Z(zone(k).knotG(i,j))];
                        %cas 4a
                        if(inpolygon(PCr,PCz,triR,triZ))
                            zone(k).knotB(i,j)=0;
                            zone(k).knotE(i,j)=0;
                            zone(k).knotA(i,j)=0;
                            zone(k).knotH(i,j)=0;
                            zone(k).knotD(i,j)=0;
                            % on crée 1 triangle
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotF(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotC(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotG(i,j);
                        else %cas 4b
                            zone(k).knotC(i,j)=0;
                            % on crée 3 triangles
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotB(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotF(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotA(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotA(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotF(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotD(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotF(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotG(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotD(i,j);
                        end
                    end
                    %cas 5
                    if(zone(k).knotH(i,j)~=0)
                        triR=[R(zone(k).knotH(i,j));...
                            R(zone(k).knotG(i,j));...
                            R(zone(k).knotD(i,j))];
                        triZ=[Z(zone(k).knotH(i,j));...
                            Z(zone(k).knotG(i,j));...
                            Z(zone(k).knotD(i,j))];
                        %cas 5a
                        if(inpolygon(PCr,PCz,triR,triZ))
                            zone(k).knotD(i,j)=0;
                            zone(k).knotE(i,j)=0;
                            zone(k).knotB(i,j)=0;
                            zone(k).knotF(i,j)=0;
                            zone(k).knotC(i,j)=0;
                            % on crée 1 triangle
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotH(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotG(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotD(i,j);
                        else %cas 5b
                            zone(k).knotD(i,j)=0;
                            % on crée 3 triangles
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotA(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotB(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotH(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotB(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotC(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotH(i,j);
                            ntriangles=ntriangles+1;
                            triangles(ntriangles).p(1)=zone(k).knotH(i,j);
                            triangles(ntriangles).p(2)=zone(k).knotC(i,j);
                            triangles(ntriangles).p(3)=zone(k).knotG(i,j);
                        end
                    end
                end
                %cas6
                if((zone(k).knotF(i,j)~=0)&&(zone(k).knotH(i,j)~=0))
                    triR=[R(zone(k).knotB(i,j));...
                        R(zone(k).knotF(i,j));...
                        R(zone(k).knotH(i,j));...
                        R(zone(k).knotA(i,j))];
                    triZ=[Z(zone(k).knotB(i,j));...
                        Z(zone(k).knotF(i,j));...
                        Z(zone(k).knotH(i,j));...
                        Z(zone(k).knotA(i,j))];
                    %cas 3a
                    if(inpolygon(PCr,PCz,triR,triZ))
                        zone(k).knotC(i,j)=0;
                        zone(k).knotG(i,j)=0;
                        zone(k).knotD(i,j)=0;
                        % on crée 2 triangles
                        ntriangles=ntriangles+1;
                        triangles(ntriangles).p(1)=zone(k).knotB(i,j);
                        triangles(ntriangles).p(2)=zone(k).knotH(i,j);
                        triangles(ntriangles).p(3)=zone(k).knotA(i,j);
                        ntriangles=ntriangles+1;
                        triangles(ntriangles).p(1)=zone(k).knotB(i,j);
                        triangles(ntriangles).p(2)=zone(k).knotF(i,j);
                        triangles(ntriangles).p(3)=zone(k).knotH(i,j);
                    else %cas 3b
                        zone(k).knotB(i,j)=0;
                        zone(k).knotE(i,j)=0;
                        zone(k).knotA(i,j)=0;
                        % on crée 2 triangles
                        ntriangles=ntriangles+1;
                        triangles(ntriangles).p(1)=zone(k).knotF(i,j);
                        triangles(ntriangles).p(2)=zone(k).knotD(i,j);
                        triangles(ntriangles).p(3)=zone(k).knotH(i,j);
                        ntriangles=ntriangles+1;
                        triangles(ntriangles).p(1)=zone(k).knotF(i,j);
                        triangles(ntriangles).p(2)=zone(k).knotC(i,j);
                        triangles(ntriangles).p(3)=zone(k).knotD(i,j);
                    end
                end
            end
        end
    end
end

figure
hold on
for i=1:2
    plot(R([triangles(i).p(1) triangles(i).p(2)]),Z([triangles(i).p(1) triangles(i).p(2)]),'k-')
    plot(R([triangles(i).p(1) triangles(i).p(3)]),Z([triangles(i).p(1) triangles(i).p(3)]),'k-')
    plot(R([triangles(i).p(2) triangles(i).p(3)]),Z([triangles(i).p(2) triangles(i).p(3)]),'k-')
end

