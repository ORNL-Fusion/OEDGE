global zones

nzones=zones.num;

Nxv=[];
Nzv=[];
for k=1:zones.num
    Nxv=[Nxv,zones.zone(k).Nx];
    Nzv=[Nzv,zones.zone(k).Nz];
end
Nxm=max(Nxv);
Nzm=max(Nzv);

d2=ones(Nxm,Nzm,nzones)*10000;
% treating pass 3 points
for n=1:nknots_
    if(knots_interp(n).pass==3)
        if(knots_interp(n).nsol==2) %one needs to provide one more point
            % looking for closest eirene points
            d=zeros(1,nknots_);
            for n2=1:nknots_
                if(knots_interp(n2).pass==1)
                    d(n2)=sqrt((R_(n)-R_(n2))^2+(Z_(n)-Z_(n2))^2);
                else
                    d(n2)=10000;
                end
            end
            d(n)=10000;
            [deir,p]=min(d);
            knots_interp(n).neir=1;
            knots_interp(n).eir=p;
        end
    end
end

for n=1:nknots_
    if(knots_interp(n).pass==3)
        if(knots_interp(n).nsol==1)
            d=zeros(1,nknots_);
            for n2=1:nknots_
                if(knots_interp(n2).pass==1)
                    d(n2)=sqrt((R_(n)-R_(n2))^2+(Z_(n)-Z_(n2))^2);
                else
                    d(n2)=10000;
                end
            end
            d(n)=10000;
            [m,p]=min(d);
            d(p)=10000;
            [m2,p2]=min(d);
            dmins=[m,m2];
            % looking for closest soledge points
            for k=1:nzones
                for i=1:zones.zone(k).Nx
                    for j=1:zones.zone(k).Nz
                        if(zones.zone(k).chi(i,j)==0)
                            d2(i,j,k)=sqrt((R_(n)-zones.zone(k).gridRc(i,j))^2+(Z_(n)-zones.zone(k).gridZc(i,j))^2);
                        else
                            d2(i,j,k)=10000;
                        end
                    end
                end
            end
            k=knots_interp(n).sol(1,1);
            i=knots_interp(n).sol(1,2);
            j=knots_interp(n).sol(1,3);
            d2(i,j,k)=10000;
            [i1,j1,k1]=find(d2==min(min(min(d2))));
            j1_=mod(j1,Nzm);
            if(j1_~=0)
                k1=(j1-j1_)/Nzm+1;
                j1=j1_;
            else
                k1=j1/Nzm;
                j1=Nzm;
            end
            dsol=min(min(min(d2)));
            dmins=[dmins,dsol];
            %taking the two closest points
            [A,B]=sort(dmins);
            switch B(1)
                case 1
                    knots_interp(n).neir=knots_interp(n).neir+1;
                    knots_interp(n).eir=[knots_interp(n).eir,p];
                case 2 %impossible
                    knots_interp(n).neir=knots_interp(n).neir+1;
                    knots_interp(n).eir=[knots_interp(n).eir,p2];
                case 3
                    knots_interp(n).nsol=knots_interp(n).nsol+1;
                    knots_interp(n).sol=[knots_interp(n).sol;k1,i1,j1];
            end
            switch B(2)
                case 1
                    knots_interp(n).neir=knots_interp(n).neir+1;
                    knots_interp(n).eir=[knots_interp(n).eir,p];
                case 2
                    knots_interp(n).neir=knots_interp(n).neir+1;
                    knots_interp(n).eir=[knots_interp(n).eir,p2];
                case 3
                    knots_interp(n).nsol=knots_interp(n).nsol+1;
                    knots_interp(n).sol=[knots_interp(n).sol;k1,i1,j1];
            end
        end
    end
end

% treating pass 4 points
for n=1:nknots_
    if(knots_interp(n).pass==4)
        %one needs to provide three points
        d=zeros(1,nknots_);
        for n2=1:nknots_
            if(knots_interp(n2).pass~=4)
                d(n2)=sqrt((R_(n)-R_(n2))^2+(Z_(n)-Z_(n2))^2);
            else
                d(n2)=10000;
            end
        end
        d(n)=10000;
        [m,p]=min(d);
        d(p)=10000;
        [m2,p2]=min(d);
        dmins=[m,m2];
        % looking for closest soledge points
        for k=1:nzones
            for i=1:zones.zone(k).Nx
                for j=1:zones.zone(k).Nz
                    if(zones.zone(k).chi(i,j)==0)
                        d2(i,j,k)=sqrt((R_(n)-zones.zone(k).gridRc(i,j))^2+(Z_(n)-zones.zone(k).gridZc(i,j))^2);
                    else
                        d2(i,j,k)=10000;
                    end
                end
            end
        end
        [i1,j1,k1]=find(d2==min(min(min(d2))));
        j1_=mod(j1,Nzm);
        if(j1_~=0)
            k1=(j1-j1_)/Nzm+1;
            j1=j1_;
        else
            k1=j1/Nzm;
            j1=Nzm;
        end
        dsol=min(min(min(d2)));
        dmins=[dmins,dsol];
        d2(i1,j1,k1)=10000;
        [i2,j2,k2]=find(d2==min(min(min(d2))));
        j2_=mod(j2,Nzm);
        if(j2_~=0)
            k2=(j2-j2_)/Nzm+1;
            j2=j2_;
        else
            k2=j2/Nzm;
            j2=Nzm;
        end
        dsol=min(min(min(d2)));
        dmins=[dmins,dsol];
        %taking the three closest points
        [A,B]=sort(dmins);
        switch B(1)
            case 1
                knots_interp(n).neir=knots_interp(n).neir+1;
                knots_interp(n).eir=[knots_interp(n).eir,p];
            case 2 %impossible
                knots_interp(n).neir=knots_interp(n).neir+1;
                knots_interp(n).eir=[knots_interp(n).eir,p2];
            case 3
                knots_interp(n).nsol=knots_interp(n).nsol+1;
                knots_interp(n).sol=[knots_interp(n).sol;k1,i1,j1];
            case 4 %impossible
                knots_interp(n).nsol=knots_interp(n).nsol+1;
                knots_interp(n).sol=[knots_interp(n).sol;k2,i2,j2];
        end
        switch B(2)
            case 1
                knots_interp(n).neir=knots_interp(n).neir+1;
                knots_interp(n).eir=[knots_interp(n).eir,p];
            case 2
                knots_interp(n).neir=knots_interp(n).neir+1;
                knots_interp(n).eir=[knots_interp(n).eir,p2];
            case 3
                knots_interp(n).nsol=knots_interp(n).nsol+1;
                knots_interp(n).sol=[knots_interp(n).sol;k1,i1,j1];
            case 4
                knots_interp(n).nsol=knots_interp(n).nsol+1;
                knots_interp(n).sol=[knots_interp(n).sol;k2,i2,j2];
        end
        switch B(3)
            case 1
                knots_interp(n).neir=knots_interp(n).neir+1;
                knots_interp(n).eir=[knots_interp(n).eir,p];
            case 2
                knots_interp(n).neir=knots_interp(n).neir+1;
                knots_interp(n).eir=[knots_interp(n).eir,p2];
            case 3
                knots_interp(n).nsol=knots_interp(n).nsol+1;
                knots_interp(n).sol=[knots_interp(n).sol;k1,i1,j1];
            case 4
                knots_interp(n).nsol=knots_interp(n).nsol+1;
                knots_interp(n).sol=[knots_interp(n).sol;k2,i2,j2];
        end
    end
end