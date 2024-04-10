for n=1:nknots_
    [m,p]=find(tri_knots_==n);
    ntri_in=length(m);
    sol_mesh=[];
    s=0;
    for n2=1:ntri_in
        sol_mesh_=[triangles_(m(n2)).k triangles_(m(n2)).i triangles_(m(n2)).j];
        in=0;
        for n3=1:s
            if(sum(sol_mesh(n3,:)==sol_mesh_)==3)
                in=1;
            end
        end
        if(in==0)
            sol_mesh=[sol_mesh;sol_mesh_];
            s=s+1;
        end
    end
    % selection of the soledge quadrangle in the plasma
    sol_mesh2=[];
    s2=0;
    for n2=1:s
        k=sol_mesh(n2,1);
        i=sol_mesh(n2,2);
        j=sol_mesh(n2,3);
        if(k~=0)
            if(zones.zone(k).chi(i,j)==0)
                sol_mesh2=[sol_mesh2;sol_mesh(n2,:)];
                s2=s2+1;
            end
        end
    end
    if(s2>=4) % The knot is connected with more than 4 soledge quadrangle in the plasma => first pass candidate
        knots_interp(n).pass=1;
        knots_interp(n).nsol=s2;
        knots_interp(n).sol=sol_mesh2;
        knots_interp(n).neir=0;
        knots_interp(n).eir=[];
    else
        if(s2==3) % The knot is connected with 3 soledge quadrangle in the plasma => second pass candidate
            knots_interp(n).pass=2;
            knots_interp(n).nsol=s2;
            knots_interp(n).sol=sol_mesh2;
            knots_interp(n).neir=0;
            knots_interp(n).eir=[];
        else
            if(s2>=1)
                knots_interp(n).pass=3;
                knots_interp(n).nsol=s2;
                knots_interp(n).sol=sol_mesh2;
                knots_interp(n).neir=0;
                knots_interp(n).eir=[];
            else
                knots_interp(n).pass=4; % no soledge point (corner in penalized quadrangle typically)
                knots_interp(n).nsol=0;
                knots_interp(n).sol=[];
                knots_interp(n).neir=0;
                knots_interp(n).eir=[];
            end
        end
    end
end



