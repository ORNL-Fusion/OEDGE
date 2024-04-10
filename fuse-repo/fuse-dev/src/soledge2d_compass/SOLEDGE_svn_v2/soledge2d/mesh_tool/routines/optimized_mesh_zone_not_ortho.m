function optimized_mesh_zone_not_ortho(nz,refR,refZ,direction)

global zones;
global megazone;
global r2D;
global z2D;
global flux2D;
global Pmegazone;
global Rwall;
global Zwall;

long=length(refR);
dist=zeros(1,long);
for k=2:long
    dist(k)=dist(k-1)+sqrt((refR(k)-refR(k-1))^2+(refZ(k)-refZ(k-1))^2);
end
dist=dist/dist(end);

vx=zones.zone(nz).south.R(2)-zones.zone(nz).south.R(1);
vy=zones.zone(nz).south.Z(2)-zones.zone(nz).south.Z(1);
vx1=zones.zone(nz).west.R(2)-zones.zone(nz).west.R(1);
vy1=zones.zone(nz).west.Z(2)-zones.zone(nz).west.Z(1);
qual=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2));
sens=sign(qual);

psilist=megazone.mz(zones.zone(nz).mz).refpoints.psi;
pmz=zones.zone(nz).pmz;

if(direction==1) %south to north
    %first one
    Rm=zones.zone(nz).south.R;
    Zm=zones.zone(nz).south.Z;
    dm=zeros(size(Rm));
    for j=2:length(Rm)
        dm(j)=dm(j-1)+sqrt((Rm(j)-Rm(j-1))^2+(Zm(j)-Zm(j-1))^2);
    end
    dm=dm/dm(end);
    for j=1:length(dist)
        zones.zone(nz).gridR(1,j)=interp1(dm,Rm,dist(j));
        zones.zone(nz).gridZ(1,j)=interp1(dm,Zm,dist(j));
    end
    zones.zone(nz).gridR(1,1)=zones.zone(nz).south.R(1);
    zones.zone(nz).gridZ(1,1)=zones.zone(nz).south.Z(1);
    zones.zone(nz).gridR(1,length(dist))=zones.zone(nz).south.R(end);
    zones.zone(nz).gridZ(1,length(dist))=zones.zone(nz).south.Z(end);
    for i=2:length(psilist)-1
        c1=contour_better(r2D,z2D,flux2D,psilist(i));
        cE.num=1;
        cE.arc(1).x=zones.zone(nz).east.R;
        cE.arc(1).y=zones.zone(nz).east.Z;
        cW.num=1;
        cW.arc(1).x=zones.zone(nz).west.R;
        cW.arc(1).y=zones.zone(nz).west.Z;
        XE=intersect_contour(c1,cE);
        XW=intersect_contour(c1,cW);
        cnum=XE.arc1(1);
        cin.x=c1.arc(cnum).x;
        cin.y=c1.arc(cnum).y;
        p1.x=XW.x(1);
        p1.y=XW.y(1);
        p2.x=XE.x(1);
        p2.y=XE.y(1);
        clear cout
        cout(1)=part_contour_per(cin,p1,p2,1);
        cout(2)=part_contour_per(cin,p1,p2,-1);
        %choose arc
        Areas=zeros(1,2);
        for na=1:2
            Areas(na)=polyarea([zones.zone(nz).gridR(i-1,:),cout(na).x(end:-1:1)],...
                [zones.zone(nz).gridZ(i-1,:),cout(na).y(end:-1:1)]);
        end
        na=find(Areas==min(Areas));
        Rm=cout(na(1)).x;
        Zm=cout(na(1)).y;
        
        if(~Pmegazone.mz(pmz).isaligned)
            Rin=zones.zone(nz).gridR(i-1,:);
            Zin=zones.zone(nz).gridZ(i-1,:);
            din=zeros(size(Rin));
            for j=2:length(Rin);
                din(j)=din(j-1)+sqrt((Rin(j)-Rin(j-1))^2+(Zin(j)-Zin(j-1))^2);
            end
            din=din/din(end);
            [dout,Rout,Zout]=optim_mesh(din,Rin,Zin,Rm,Zm,sens);
            zones.zone(nz).gridR(i,:)=Rout;
            zones.zone(nz).gridZ(i,:)=Zout;
        else
            if((psilist(i)<=Pmegazone.mz(pmz).align_psimax)&&(psilist(i)>=Pmegazone.mz(pmz).align_psimin))
                %align --> two parts
                Rin1=zones.zone(nz).gridR(i-1,1:length(Pmegazone.mz(pmz).subrefpoints(1).R));
                Zin1=zones.zone(nz).gridZ(i-1,1:length(Pmegazone.mz(pmz).subrefpoints(1).R));
                din1=zeros(size(Rin1));
                for j=2:length(Rin1);
                    din1(j)=din1(j-1)+sqrt((Rin1(j)-Rin1(j-1))^2+(Zin1(j)-Zin1(j-1))^2);
                end
                din1=din1/din1(end);
                c1.num=1;
                c1.arc(1).x=Rwall;
                c1.arc(1).y=Zwall;
                c2.num=1;
                c2.arc(1).x=Rm;
                c2.arc(1).y=Zm;
                X=intersect_contour(c1,c2);
                cin.x=Rm;
                cin.y=Zm;
                if(X.num>1)
                    d=zeros(X.num,1);
                    if(inpolygon(c2.arc(1).x(1),c2.arc(1).y(1),Rwall,Zwall))
                        for k=1:X.num
                            d(k)=sqrt((X.x(k)- c2.arc(1).x(1))^2+(X.y(k)-c2.arc(1).y(1))^2);
                        end
                        a=find(d==min(d));
                        X.num=1;
                        X.x(1)=X.x(a);
                        X.y(1)=X.y(a);
                    else
                        for k=1:X.num
                            d(k)=sqrt((X.x(k)- c2.arc(1).x(end))^2+(X.y(k)-c2.arc(1).y(end))^2);
                        end
                        a=find(d==min(d));
                        X.num=1;
                        X.x(1)=X.x(a);
                        X.y(1)=X.y(a);
                    end
                end
                if(X.num==1)
                    p1.x=X.x(1);
                    p1.y=X.y(1);
                end
                cout=part_contour(cin,p1);
                Rm1=cout.arc(1).x;
                Zm1=cout.arc(1).y;
                Rm2=cout.arc(2).x;
                Zm2=cout.arc(2).y;
                [dout,Rout1,Zout1]=optim_mesh(din1,Rin1,Zin1,Rm1,Zm1,sens);
                Rin2=zones.zone(nz).gridR(i-1,length(Pmegazone.mz(pmz).subrefpoints(1).R):end);
                Zin2=zones.zone(nz).gridZ(i-1,length(Pmegazone.mz(pmz).subrefpoints(1).R):end);
                din2=zeros(size(Rin2));
                for j=2:length(Rin2);
                    din2(j)=din2(j-1)+sqrt((Rin2(j)-Rin2(j-1))^2+(Zin2(j)-Zin2(j-1))^2);
                end
                din2=din2/din2(end);
                [dout,Rout2,Zout2]=optim_mesh(din2,Rin2,Zin2,Rm2,Zm2,sens);
                Rout=[Rout1(1:end-1),Rout2];
                Zout=[Zout1(1:end-1),Zout2];
                zones.zone(nz).gridR(i,:)=Rout;
                zones.zone(nz).gridZ(i,:)=Zout;
            else
                %normal
                Rin=zones.zone(nz).gridR(i-1,:);
                Zin=zones.zone(nz).gridZ(i-1,:);
                din=zeros(size(Rin));
                for j=2:length(Rin);
                    din(j)=din(j-1)+sqrt((Rin(j)-Rin(j-1))^2+(Zin(j)-Zin(j-1))^2);
                end
                din=din/din(end);
                [dout,Rout,Zout]=optim_mesh(din,Rin,Zin,Rm,Zm,sens);
                zones.zone(nz).gridR(i,:)=Rout;
                zones.zone(nz).gridZ(i,:)=Zout;
            end
        end
        
    end
    %last one
    i=length(psilist);
    Rm=zones.zone(nz).north.R;
    Zm=zones.zone(nz).north.Z;
    if(~Pmegazone.mz(pmz).isaligned)
        Rin=zones.zone(nz).gridR(i-1,:);
        Zin=zones.zone(nz).gridZ(i-1,:);
        din=zeros(size(Rin));
        for j=2:length(Rin);
            din(j)=din(j-1)+sqrt((Rin(j)-Rin(j-1))^2+(Zin(j)-Zin(j-1))^2);
        end
        din=din/din(end);
        [dout,Rout,Zout]=optim_mesh(din,Rin,Zin,Rm,Zm,sens);
        zones.zone(nz).gridR(i,:)=Rout;
        zones.zone(nz).gridZ(i,:)=Zout;
    else
        if((psilist(i)<=Pmegazone.mz(pmz).align_psimax)&&(psilist(i)>=Pmegazone.mz(pmz).align_psimin))
            %align --> two parts
            Rin1=zones.zone(nz).gridR(i-1,1:length(Pmegazone.mz(pmz).subrefpoints(1).R));
            Zin1=zones.zone(nz).gridZ(i-1,1:length(Pmegazone.mz(pmz).subrefpoints(1).R));
            din1=zeros(size(Rin1));
            for j=2:length(Rin1);
                din1(j)=din1(j-1)+sqrt((Rin1(j)-Rin1(j-1))^2+(Zin1(j)-Zin1(j-1))^2);
            end
            din1=din1/din1(end);
            c1.num=1;
            c1.arc(1).x=Rwall;
            c1.arc(1).y=Zwall;
            c2.num=1;
            c2.arc(1).x=Rm;
            c2.arc(1).y=Zm;
            X=intersect_contour(c1,c2);
            cin.x=Rm;
            cin.y=Zm;
            if(X.num>1)
                d=zeros(X.num,1);
                if(inpolygon(c2.arc(1).x(1),c2.arc(1).y(1),Rwall,Zwall))
                    for k=1:X.num
                        d(k)=sqrt((X.x(k)- c2.arc(1).x(1))^2+(X.y(k)-c2.arc(1).y(1))^2);
                    end
                    a=find(d==min(d));
                    X.num=1;
                    X.x(1)=X.x(a);
                    X.y(1)=X.y(a);
                else
                    for k=1:X.num
                        d(k)=sqrt((X.x(k)- c2.arc(1).x(end))^2+(X.y(k)-c2.arc(1).y(end))^2);
                    end
                    a=find(d==min(d));
                    X.num=1;
                    X.x(1)=X.x(a);
                    X.y(1)=X.y(a);
                end
            end
            if(X.num==1)
                p1.x=X.x(1);
                p1.y=X.y(1);
            end
            cout=part_contour(cin,p1);
            Rm1=cout.arc(1).x;
            Zm1=cout.arc(1).y;
            Rm2=cout.arc(2).x;
            Zm2=cout.arc(2).y;
            [dout,Rout1,Zout1]=optim_mesh(din1,Rin1,Zin1,Rm1,Zm1,sens);
            Rin2=zones.zone(nz).gridR(i-1,length(Pmegazone.mz(pmz).subrefpoints(1).R):end);
            Zin2=zones.zone(nz).gridZ(i-1,length(Pmegazone.mz(pmz).subrefpoints(1).R):end);
            din2=zeros(size(Rin2));
            for j=2:length(Rin2);
                din2(j)=din2(j-1)+sqrt((Rin2(j)-Rin2(j-1))^2+(Zin2(j)-Zin2(j-1))^2);
            end
            din2=din2/din2(end);
            [dout,Rout2,Zout2]=optim_mesh(din2,Rin2,Zin2,Rm2,Zm2,sens);
            Rout=[Rout1(1:end-1),Rout2];
            Zout=[Zout1(1:end-1),Zout2];
            zones.zone(nz).gridR(i,:)=Rout;
            zones.zone(nz).gridZ(i,:)=Zout;
        else
            %normal
            Rin=zones.zone(nz).gridR(i-1,:);
            Zin=zones.zone(nz).gridZ(i-1,:);
            din=zeros(size(Rin));
            for j=2:length(Rin);
                din(j)=din(j-1)+sqrt((Rin(j)-Rin(j-1))^2+(Zin(j)-Zin(j-1))^2);
            end
            din=din/din(end);
            [dout,Rout,Zout]=optim_mesh(din,Rin,Zin,Rm,Zm,sens);
            zones.zone(nz).gridR(i,:)=Rout;
            zones.zone(nz).gridZ(i,:)=Zout;
        end
    end
    
    zones.zone(nz).gridR(length(psilist),1)=zones.zone(nz).north.R(1);
    zones.zone(nz).gridZ(length(psilist),1)=zones.zone(nz).north.Z(1);
    zones.zone(nz).gridR(length(psilist),length(dist))=zones.zone(nz).north.R(end);
    zones.zone(nz).gridZ(length(psilist),length(dist))=zones.zone(nz).north.Z(end);
    
else %north to south
    %last one
    Rm=zones.zone(nz).north.R;
    Zm=zones.zone(nz).north.Z;
    dm=zeros(size(Rm));
    for j=2:length(Rm)
        dm(j)=dm(j-1)+sqrt((Rm(j)-Rm(j-1))^2+(Zm(j)-Zm(j-1))^2);
    end
    dm=dm/dm(end);
    for j=1:length(dist)
        zones.zone(nz).gridR(length(psilist),j)=interp1(dm,Rm,dist(j));
        zones.zone(nz).gridZ(length(psilist),j)=interp1(dm,Zm,dist(j));
    end
    for i=length(psilist)-1:-1:2
        c1=contour_better(r2D,z2D,flux2D,psilist(i));
        cE.num=1;
        cE.arc(1).x=zones.zone(nz).east.R;
        cE.arc(1).y=zones.zone(nz).east.Z;
        cW.num=1;
        cW.arc(1).x=zones.zone(nz).west.R;
        cW.arc(1).y=zones.zone(nz).west.Z;
        XE=intersect_contour(c1,cE);
        XW=intersect_contour(c1,cW);
        cnum=XE.arc1(1);
        cin.x=c1.arc(cnum).x;
        cin.y=c1.arc(cnum).y;
        p1.x=XW.x(1);
        p1.y=XW.y(1);
        p2.x=XE.x(1);
        p2.y=XE.y(1);
        clear cout
        cout(1)=part_contour_per(cin,p1,p2,1);
        cout(2)=part_contour_per(cin,p1,p2,-1);
        %choose arc
        Areas=zeros(1,2);
        for na=1:2
            Areas(na)=polyarea([zones.zone(nz).gridR(i+1,:),cout(na).x(end:-1:1)],...
                [zones.zone(nz).gridZ(i+1,:),cout(na).y(end:-1:1)]);
        end
        na=find(Areas==min(Areas));
        Rm=cout(na(1)).x;
        Zm=cout(na(1)).y;
        
        if(~Pmegazone.mz(pmz).isaligned)
            Rin=zones.zone(nz).gridR(i+1,:);
            Zin=zones.zone(nz).gridZ(i+1,:);
            din=zeros(size(Rin));
            for j=2:length(Rin);
                din(j)=din(j-1)+sqrt((Rin(j)-Rin(j-1))^2+(Zin(j)-Zin(j-1))^2);
            end
            din=din/din(end);
            [dout,Rout,Zout]=optim_mesh(din,Rin,Zin,Rm,Zm,-sens);
            zones.zone(nz).gridR(i,:)=Rout;
            zones.zone(nz).gridZ(i,:)=Zout;
        else
            if((psilist(i)<=Pmegazone.mz(pmz).align_psimax)&&(psilist(i)>=Pmegazone.mz(pmz).align_psimin))
                %align --> two parts
                Rin1=zones.zone(nz).gridR(i+1,1:length(Pmegazone.mz(pmz).subrefpoints(1).R));
                Zin1=zones.zone(nz).gridZ(i+1,1:length(Pmegazone.mz(pmz).subrefpoints(1).R));
                din1=zeros(size(Rin1));
                for j=2:length(Rin1);
                    din1(j)=din1(j-1)+sqrt((Rin1(j)-Rin1(j-1))^2+(Zin1(j)-Zin1(j-1))^2);
                end
                din1=din1/din1(end);
                c1.num=1;
                c1.arc(1).x=Rwall;
                c1.arc(1).y=Zwall;
                c2.num=1;
                c2.arc(1).x=Rm;
                c2.arc(1).y=Zm;
                X=intersect_contour(c1,c2);
                cin.x=Rm;
                cin.y=Zm;
                if(X.num>1)
                    d=zeros(X.num,1);
                    if(inpolygon(c2.arc(1).x(1),c2.arc(1).y(1),Rwall,Zwall))
                        for k=1:X.num
                            d(k)=sqrt((X.x(k)- c2.arc(1).x(1))^2+(X.y(k)-c2.arc(1).y(1))^2);
                        end
                        a=find(d==min(d));
                        X.num=1;
                        X.x(1)=X.x(a);
                        X.y(1)=X.y(a);
                    else
                        for k=1:X.num
                            d(k)=sqrt((X.x(k)- c2.arc(1).x(end))^2+(X.y(k)-c2.arc(1).y(end))^2);
                        end
                        a=find(d==min(d));
                        X.num=1;
                        X.x(1)=X.x(a);
                        X.y(1)=X.y(a);
                    end
                end
                if(X.num==1)
                    p1.x=X.x(1);
                    p1.y=X.y(1);
                end
                cout=part_contour(cin,p1);
                Rm1=cout.arc(1).x;
                Zm1=cout.arc(1).y;
                Rm2=cout.arc(2).x;
                Zm2=cout.arc(2).y;
                [dout,Rout1,Zout1]=optim_mesh(din1,Rin1,Zin1,Rm1,Zm1,-sens);
                Rin2=zones.zone(nz).gridR(i+1,length(Pmegazone.mz(pmz).subrefpoints(1).R):end);
                Zin2=zones.zone(nz).gridZ(i+1,length(Pmegazone.mz(pmz).subrefpoints(1).R):end);
                din2=zeros(size(Rin2));
                for j=2:length(Rin2);
                    din2(j)=din2(j-1)+sqrt((Rin2(j)-Rin2(j-1))^2+(Zin2(j)-Zin2(j-1))^2);
                end
                din2=din2/din2(end);
                [dout,Rout2,Zout2]=optim_mesh(din2,Rin2,Zin2,Rm2,Zm2,-sens);
                Rout=[Rout1(1:end-1),Rout2];
                Zout=[Zout1(1:end-1),Zout2];
                zones.zone(nz).gridR(i,:)=Rout;
                zones.zone(nz).gridZ(i,:)=Zout;
            else
                %normal
                Rin=zones.zone(nz).gridR(i+1,:);
                Zin=zones.zone(nz).gridZ(i+1,:);
                din=zeros(size(Rin));
                for j=2:length(Rin);
                    din(j)=din(j-1)+sqrt((Rin(j)-Rin(j-1))^2+(Zin(j)-Zin(j-1))^2);
                end
                din=din/din(end);
                [dout,Rout,Zout]=optim_mesh(din,Rin,Zin,Rm,Zm,-sens);
                zones.zone(nz).gridR(i,:)=Rout;
                zones.zone(nz).gridZ(i,:)=Zout;
            end
        end
    end
    %first one
    i=1;
    Rm=zones.zone(nz).south.R;
    Zm=zones.zone(nz).south.Z;
    if(~Pmegazone.mz(pmz).isaligned)
        Rin=zones.zone(nz).gridR(i+1,:);
        Zin=zones.zone(nz).gridZ(i+1,:);
        din=zeros(size(Rin));
        for j=2:length(Rin);
            din(j)=din(j-1)+sqrt((Rin(j)-Rin(j-1))^2+(Zin(j)-Zin(j-1))^2);
        end
        din=din/din(end);
        [dout,Rout,Zout]=optim_mesh(din,Rin,Zin,Rm,Zm,-sens);
        zones.zone(nz).gridR(i,:)=Rout;
        zones.zone(nz).gridZ(i,:)=Zout;
    else
        if((psilist(i)<=Pmegazone.mz(pmz).align_psimax)&&(psilist(i)>=Pmegazone.mz(pmz).align_psimin))
            %align --> two parts
            Rin1=zones.zone(nz).gridR(i+1,1:length(Pmegazone.mz(pmz).subrefpoints(1).R));
            Zin1=zones.zone(nz).gridZ(i+1,1:length(Pmegazone.mz(pmz).subrefpoints(1).R));
            din1=zeros(size(Rin1));
            for j=2:length(Rin1);
                din1(j)=din1(j-1)+sqrt((Rin1(j)-Rin1(j-1))^2+(Zin1(j)-Zin1(j-1))^2);
            end
            din1=din1/din1(end);
            c1.num=1;
            c1.arc(1).x=Rwall;
            c1.arc(1).y=Zwall;
            c2.num=1;
            c2.arc(1).x=Rm;
            c2.arc(1).y=Zm;
            X=intersect_contour(c1,c2);
            cin.x=Rm;
            cin.y=Zm;
            if(X.num>1)
                d=zeros(X.num,1);
                if(inpolygon(c2.arc(1).x(1),c2.arc(1).y(1),Rwall,Zwall))
                    for k=1:X.num
                        d(k)=sqrt((X.x(k)- c2.arc(1).x(1))^2+(X.y(k)-c2.arc(1).y(1))^2);
                    end
                    a=find(d==min(d));
                    X.num=1;
                    X.x(1)=X.x(a);
                    X.y(1)=X.y(a);
                else
                    for k=1:X.num
                        d(k)=sqrt((X.x(k)- c2.arc(1).x(end))^2+(X.y(k)-c2.arc(1).y(end))^2);
                    end
                    a=find(d==min(d));
                    X.num=1;
                    X.x(1)=X.x(a);
                    X.y(1)=X.y(a);
                end
            end
            if(X.num==1)
                p1.x=X.x(1);
                p1.y=X.y(1);
            end
            cout=part_contour(cin,p1);
            Rm1=cout.arc(1).x;
            Zm1=cout.arc(1).y;
            Rm2=cout.arc(2).x;
            Zm2=cout.arc(2).y;
            [dout,Rout1,Zout1]=optim_mesh(din1,Rin1,Zin1,Rm1,Zm1,-sens);
            Rin2=zones.zone(nz).gridR(i+1,length(Pmegazone.mz(pmz).subrefpoints(1).R):end);
            Zin2=zones.zone(nz).gridZ(i+1,length(Pmegazone.mz(pmz).subrefpoints(1).R):end);
            din2=zeros(size(Rin2));
            for j=2:length(Rin2);
                din2(j)=din2(j-1)+sqrt((Rin2(j)-Rin2(j-1))^2+(Zin2(j)-Zin2(j-1))^2);
            end
            din2=din2/din2(end);
            [dout,Rout2,Zout2]=optim_mesh(din2,Rin2,Zin2,Rm2,Zm2,-sens);
            Rout=[Rout1(1:end-1),Rout2];
            Zout=[Zout1(1:end-1),Zout2];
            zones.zone(nz).gridR(i,:)=Rout;
            zones.zone(nz).gridZ(i,:)=Zout;
        else
            %normal
            Rin=zones.zone(nz).gridR(i+1,:);
            Zin=zones.zone(nz).gridZ(i+1,:);
            din=zeros(size(Rin));
            for j=2:length(Rin);
                din(j)=din(j-1)+sqrt((Rin(j)-Rin(j-1))^2+(Zin(j)-Zin(j-1))^2);
            end
            din=din/din(end);
            [dout,Rout,Zout]=optim_mesh(din,Rin,Zin,Rm,Zm,-sens);
            zones.zone(nz).gridR(i,:)=Rout;
            zones.zone(nz).gridZ(i,:)=Zout;
        end
    end
    
    zones.zone(nz).gridR(1,1)=zones.zone(nz).south.R(1);
    zones.zone(nz).gridZ(1,1)=zones.zone(nz).south.Z(1);
    zones.zone(nz).gridR(1,length(dist))=zones.zone(nz).south.R(end);
    zones.zone(nz).gridZ(1,length(dist))=zones.zone(nz).south.Z(end);
    
end

