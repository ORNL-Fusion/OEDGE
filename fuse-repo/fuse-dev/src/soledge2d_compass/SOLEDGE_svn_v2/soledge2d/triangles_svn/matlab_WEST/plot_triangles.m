figure(1)
hold on
axis equal
for i=1:16
    Nx=zone(i).Nx;
    Nz=zone(i).Nz;
    for k=1:Nx
        for j=1:Nz
            %premier triangle
            nt1=zone(i).ntrinum(k,j);
            ind=find(tri==nt1);
            R1=R(ind);
            Z1=Z(ind);
            nt2=zone(i).ntrinum(k+1,j);
            ind=find(tri==nt2);
            R2=R(ind);
            Z2=Z(ind);
            nt3=zone(i).ntrinum(k,j+1);
            ind=find(tri==nt3);
            R3=R(ind);
            Z3=Z(ind);
            plot([R1,R2],[Z1,Z2],'k-');
            plot([R1,R3],[Z1,Z3],'k-');
            plot([R3,R2],[Z3,Z2],'k-');
            %deuxieme triangle
            nt1=zone(i).ntrinum(k+1,j+1);
            ind=find(tri==nt1);
            R1=R(ind);
            Z1=Z(ind);
            nt2=zone(i).ntrinum(k+1,j);
            ind=find(tri==nt2);
            R2=R(ind);
            Z2=Z(ind);
            nt3=zone(i).ntrinum(k,j+1);
            ind=find(tri==nt3);
            R3=R(ind);
            Z3=Z(ind);
            plot([R1,R2],[Z1,Z2],'k-');
            plot([R1,R3],[Z1,Z3],'k-');
            plot([R3,R2],[Z3,Z2],'k-');
        end
    end
                drawnow
end