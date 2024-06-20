global zones
global eirene

for k=1:zones.num
    for i=1:zones.zone(k).Nx
        for j=1:zones.zone(k).Nz
            if(~(zones.zone(k).isaligned(i,j)==1))
                if((sum(zones.zone(k).iscrossed(i,j,:))~=0)&&(zones.zone(k).inplasma(i,j)==1))
                    %cas 1 et 2
                    if(zones.zone(k).knotE(i,j)~=0)
                        %cas 1
                        if(zones.zone(k).knotH(i,j)~=0)
                            %cas 1a
                            if(zones.zone(k).knotA(i,j)~=0)
                                % on crée 1 triangle
                                ntriangles=ntriangles+1;
                                triangles(ntriangles).k=k;
                                triangles(ntriangles).i=i;
                                triangles(ntriangles).j=j;
                                triangles(ntriangles).p1=zones.zone(k).knotA(i,j);
                                triangles(ntriangles).p2=zones.zone(k).knotE(i,j);
                                triangles(ntriangles).p3=zones.zone(k).knotH(i,j);
                                triangles(ntriangles).area=0; % 0 on the wall
                                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                            else %cas 1b
                                %                             on crée 3 triangles
                                if(zones.zone(k).knotB(i,j)~=0)
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotE(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotB(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotH(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                                if((zones.zone(k).knotC(i,j)~=0)&&(zones.zone(k).knotB(i,j)~=0))
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotB(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotC(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotH(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                                if((zones.zone(k).knotC(i,j)~=0)&&(zones.zone(k).knotD(i,j)~=0))
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotH(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotC(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotD(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                            end
                        end
                        %cas 2
                        if(zones.zone(k).knotF(i,j)~=0)
                            %cas 2a
                            if(zones.zone(k).knotB(i,j)~=0)
                                % on crée 1 triangle
                                ntriangles=ntriangles+1;
                                triangles(ntriangles).k=k;
                                triangles(ntriangles).i=i;
                                triangles(ntriangles).j=j;
                                triangles(ntriangles).p1=zones.zone(k).knotB(i,j);
                                triangles(ntriangles).p2=zones.zone(k).knotF(i,j);
                                triangles(ntriangles).p3=zones.zone(k).knotE(i,j);
                                triangles(ntriangles).area=0; % 0 on the wall
                                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                            else %cas 2b
                                % on crée 3 triangles
                                if(zones.zone(k).knotA(i,j)~=0)
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotE(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotF(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotA(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                                if((zones.zone(k).knotA(i,j)~=0)&&(zones.zone(k).knotD(i,j)~=0))
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotA(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotF(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotD(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                                if((zones.zone(k).knotC(i,j)~=0)&&(zones.zone(k).knotD(i,j)~=0))
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotF(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotC(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotD(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                            end
                        end
                        %cas 3
                        if(zones.zone(k).knotG(i,j)~=0)
                            %cas 3a
                            if(zones.zone(k).knotB(i,j)~=0)
                                % on crée 2 triangles
                                ntriangles=ntriangles+1;
                                triangles(ntriangles).k=k;
                                triangles(ntriangles).i=i;
                                triangles(ntriangles).j=j;
                                triangles(ntriangles).p1=zones.zone(k).knotB(i,j);
                                triangles(ntriangles).p2=zones.zone(k).knotG(i,j);
                                triangles(ntriangles).p3=zones.zone(k).knotE(i,j);
                                triangles(ntriangles).area=0; % 0 on the wall
                                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                if(zones.zone(k).knotC(i,j)~=0)
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotB(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotC(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotG(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                            else %cas 3b
                                % on crée 2 triangles
                                if((zones.zone(k).knotA(i,j)~=0)&&(zones.zone(k).knotD(i,j)~=0))
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotE(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotD(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotA(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                                if(zones.zone(k).knotD(i,j)~=0)
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotE(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotG(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotD(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                            end
                        end
                    end
                    if(zones.zone(k).knotG(i,j)~=0)
                        %cas4
                        if(zones.zone(k).knotF(i,j)~=0)
                            %cas 4a
                            if(zones.zone(k).knotC(i,j)~=0)
                                % on crée 1 triangle
                                ntriangles=ntriangles+1;
                                triangles(ntriangles).k=k;
                                triangles(ntriangles).i=i;
                                triangles(ntriangles).j=j;
                                triangles(ntriangles).p1=zones.zone(k).knotF(i,j);
                                triangles(ntriangles).p2=zones.zone(k).knotC(i,j);
                                triangles(ntriangles).p3=zones.zone(k).knotG(i,j);
                                triangles(ntriangles).area=0; % 0 on the wall
                                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                            else %cas 4b
                                % on crée 3 triangles
                                if((zones.zone(k).knotA(i,j)~=0)&&(zones.zone(k).knotB(i,j)~=0))
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotB(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotF(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotA(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                                if((zones.zone(k).knotA(i,j)~=0)&&(zones.zone(k).knotD(i,j)~=0))
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotA(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotF(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotD(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                                if(zones.zone(k).knotD(i,j)~=0)
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotF(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotG(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotD(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                            end
                        end
                        %cas 5
                        if(zones.zone(k).knotH(i,j)~=0)
                            %cas 5a
                            if(zones.zone(k).knotD(i,j)~=0)
                                % on crée 1 triangle
                                ntriangles=ntriangles+1;
                                triangles(ntriangles).k=k;
                                triangles(ntriangles).i=i;
                                triangles(ntriangles).j=j;
                                triangles(ntriangles).p1=zones.zone(k).knotH(i,j);
                                triangles(ntriangles).p2=zones.zone(k).knotG(i,j);
                                triangles(ntriangles).p3=zones.zone(k).knotD(i,j);
                                triangles(ntriangles).area=0; % 0 on the wall
                                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                            else %cas 5b
                                % on crée 3 triangles
                                if((zones.zone(k).knotA(i,j)~=0)&&(zones.zone(k).knotB(i,j)~=0))
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotA(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotB(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotH(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                                if((zones.zone(k).knotC(i,j)~=0)&&(zones.zone(k).knotB(i,j)~=0))
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotB(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotC(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotH(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                                if(zones.zone(k).knotC(i,j)~=0)
                                    ntriangles=ntriangles+1;
                                    triangles(ntriangles).k=k;
                                    triangles(ntriangles).i=i;
                                    triangles(ntriangles).j=j;
                                    triangles(ntriangles).p1=zones.zone(k).knotH(i,j);
                                    triangles(ntriangles).p2=zones.zone(k).knotC(i,j);
                                    triangles(ntriangles).p3=zones.zone(k).knotG(i,j);
                                    triangles(ntriangles).area=0; % 0 on the wall
                                    tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                                end
                            end
                        end
                    end
                    %cas6
                    if((zones.zone(k).knotF(i,j)~=0)&&(zones.zone(k).knotH(i,j)~=0))
                        %cas 3a
                        if(zones.zone(k).knotB(i,j)~=0)
                            % on crée 2 triangles
                            if(zones.zone(k).knotA(i,j)~=0)
                                ntriangles=ntriangles+1;
                                triangles(ntriangles).k=k;
                                triangles(ntriangles).i=i;
                                triangles(ntriangles).j=j;
                                triangles(ntriangles).p1=zones.zone(k).knotB(i,j);
                                triangles(ntriangles).p2=zones.zone(k).knotH(i,j);
                                triangles(ntriangles).p3=zones.zone(k).knotA(i,j);
                                triangles(ntriangles).area=0; % 0 on the wall
                                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                            end
                            if(zones.zone(k).knotB(i,j)~=0)
                                ntriangles=ntriangles+1;
                                triangles(ntriangles).k=k;
                                triangles(ntriangles).i=i;
                                triangles(ntriangles).j=j;
                                triangles(ntriangles).p1=zones.zone(k).knotB(i,j);
                                triangles(ntriangles).p2=zones.zone(k).knotF(i,j);
                                triangles(ntriangles).p3=zones.zone(k).knotH(i,j);
                                triangles(ntriangles).area=0; % 0 on the wall
                                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                            end
                        else %cas 3b
                            % on crée 2 triangles
                            if(zones.zone(k).knotD(i,j)~=0)
                                ntriangles=ntriangles+1;
                                triangles(ntriangles).k=k;
                                triangles(ntriangles).i=i;
                                triangles(ntriangles).j=j;
                                triangles(ntriangles).p1=zones.zone(k).knotF(i,j);
                                triangles(ntriangles).p2=zones.zone(k).knotD(i,j);
                                triangles(ntriangles).p3=zones.zone(k).knotH(i,j);
                                triangles(ntriangles).area=0; % 0 on the wall
                                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                            end
                            if((zones.zone(k).knotC(i,j)~=0)&&(zones.zone(k).knotD(i,j)~=0))
                                ntriangles=ntriangles+1;
                                triangles(ntriangles).k=k;
                                triangles(ntriangles).i=i;
                                triangles(ntriangles).j=j;
                                triangles(ntriangles).p1=zones.zone(k).knotF(i,j);
                                triangles(ntriangles).p2=zones.zone(k).knotC(i,j);
                                triangles(ntriangles).p3=zones.zone(k).knotD(i,j);
                                triangles(ntriangles).area=0; % 0 on the wall
                                tri_knots=[tri_knots;triangles(ntriangles).p1,triangles(ntriangles).p2,triangles(ntriangles).p3];
                            end
                        end
                    end
                end
            end
        end
    end
end

eirene.triangles=triangles;
eirene.tri_knots=tri_knots;
eirene.ntriangles=ntriangles;