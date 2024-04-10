clear cote;

Nseg=zone_elements.Nseg; % total number of segments
NsegX=zone_elements.Nseg-zone_elements.Nseg2; % segments from X points (two sided)

ncote=0;

for k=1:X_points.num
    
    for n=1:4
        
        %         psi_list=zeros(1,1+X_points.num);
        psi_list=zeros(1,2);
        psi_list(1)=X_points.cut(k).psilim(n);
        psi_list(2)=X_points.psi(k);
        ind=2;
        for k1=1:X_points.num
            if(k1~=k)
                %if psi_X < psi < psilim  ou bien  psilim < psi < psi_X
                if((X_points.psi(k1)-X_points.psi(k))*(X_points.psi(k1)-X_points.cut(k).psilim(n))<0)
                    ind=ind+1;
                    psi_list(ind)=X_points.psi(k1);
                end
            end
        end
        psi_list_red=psi_list(1:ind);
        psi_list_red=sort(psi_list_red);
        
        %first find intersections
        c1.num=1;
        c1.arc(1).x=X_points.cut(k).arc(n).R;
        c1.arc(1).y=X_points.cut(k).arc(n).Z;
        for k2=1:length(psi_list_red)
            if(psi_list_red(k2)==X_points.psi(k)) %special treatment here
                
                c2=contour_better(r2D,z2D,flux2D,psi_list_red(k2));
                for nb=1:4 %for each leg
                    X.num=1;
                    X.x(1)=X_points.branch(k).R(nb);
                    X.y(1)=X_points.branch(k).Z(nb);
                    p1.x=X.x(1);
                    p1.y=X.y(1);
                    number=which_contour_num(c2,p1);
                    cprog.x=c2.arc(number).x;
                    cprog.y=c2.arc(number).y;
                    d=sqrt((cprog.x-X.x(1)).^2+(cprog.y-X.y(1)).^2);
                    a=find(d==min(d));
                    if(a<length(d))
                        v1x=cprog.x(a+1)-cprog.x(a);
                        v1y=cprog.y(a+1)-cprog.y(a);
                        vx=X.x(1)-cprog.x(a);
                        vy=X.y(1)-cprog.y(a);
                        if(v1x*vx+v1y*vy>0)
                            cleft.x=circshift(cprog.x(1:end-1)',-a)';
                            cleft.y=circshift(cprog.y(1:end-1)',-a)';
                        else
                            cleft.x=circshift(cprog.x(1:end-1)',-a+1)';
                            cleft.y=circshift(cprog.y(1:end-1)',-a+1)';
                        end
                        cleft.x=[cleft.x,cleft.x(1)];
                        cleft.y=[cleft.y,cleft.y(1)];
                    else
                        cleft.x=cprog.x;
                        cleft.y=cprog.y;
                    end
                    
                    %should we turn left or right ?
                    distX=sqrt((cleft.x(1)-X_points.R(k))^2+(cleft.y(1)-X_points.Z(k))^2);
                    distX2=sqrt((cleft.x(2)-X_points.R(k))^2+(cleft.y(2)-X_points.Z(k))^2);
                    if(distX2<distX) %cleft NOK : cleft -> cright
                        sens=-1;
                        cleft.x=cleft.x(end:-1:1);
                        cleft.y=cleft.y(end:-1:1);
                    else
                        sens=1;
                    end
                    
                    
                    % then same as below
                    
                    dist=[];
                    type=[];
                    fine=[];
                    num=[];
                    for ns=1:Nseg
                        if(segments(ns).psimin==psi_list_red(k2))
                                                        if(segments(ns).psimin_numc==number)
                            X.num=1;
                            X.x(1)=segments(ns).R(1);
                            X.y(1)=segments(ns).Z(1);
                                                        else
                                                            X.num=0;
                                                        end
                        else
                            if(segments(ns).psimax==psi_list_red(k2))
                                                                if(segments(ns).psimax_numc==number)
                                X.num=1;
                                X.x(1)=segments(ns).R(end);
                                X.y(1)=segments(ns).Z(end);
                                                                else
                                                                    X.num=0;
                                                                end
                            else
                                X.num=0;
                            end
                        end
                        if(X.num==1) %intersect
                            d=sqrt((cleft.x-X.x(1)).^2+(cleft.y-X.y(1)).^2);
                            a1=find(d==min(d));
                            cleftp.x=circshift(cleft.x(1:end-1)',-1)';
                            cleftp.y=circshift(cleft.y(1:end-1)',-1)';
                            cleftp.x=[cleftp.x,cleftp.x(1)];
                            cleftp.y=[cleftp.y,cleftp.y(1)];
                            dist=[dist,a1];
                            num=[num,ns*ones(size(a1))];
                            if(ns<=NsegX)
                                type=[type,1*ones(size(a1))]; % Xpoint
                            else
                                type=[type,2*ones(size(a1))]; % lim
                            end
                            vx=X.x(1)-cleft.x(a1);
                            vy=X.y(1)-cleft.y(a1);
                            vx1=cleftp.x(a1)-cleft.x(a1);
                            vy1=cleftp.y(a1)-cleft.y(a1);
                            fine=[fine,min(d)*sign(vx.*vx1+vy.*vy1)];
                        end
                    end %end for nsegment
                    pasmax=max(abs(fine))+1;
                    [a,b]=sort(dist*pasmax+fine);
                    stopit=false;
                    ind=1;
                    segpa=[];
                    while(~stopit)
                        p1.psi=psi_list_red(k2);
                        p2.psi=psi_list_red(k2);
                        if(ind==1)%start from X point surrounding
                            p1.x=X_points.branch(k).R(nb);
                            p1.y=X_points.branch(k).Z(nb);
                            p1.theta=k;
                            p1.branch=-1;
                            p1.type=1;
                        else
                            if(segments(num(b(ind-1))).psimin==psi_list_red(k2))
                                p1.x=segments(num(b(ind-1))).R(1);
                                p1.y=segments(num(b(ind-1))).Z(1);
                            else
                                p1.x=segments(num(b(ind-1))).R(end);
                                p1.y=segments(num(b(ind-1))).Z(end);
                            end
                            if(segments(num(b(ind-1))).nX~=k)
                                p1.theta=segments(num(b(ind-1))).nX;
                                p1.branch=segments(num(b(ind-1))).nB;
                            else % modif
                                p1.theta=segments(num(b(ind-1))).nX;
                                p1.branch=-1;
                            end
                            p1.type=segments(num(b(ind-1))).type;
                        end
                        if(ind>length(b))
                            p2.x=X_points.R(k);
                            p2.y=X_points.Z(k);
                            p2.theta=k;
                            p2.branch=-1;
                            p2.type=1;
                        else
                            if(segments(num(b(ind))).psimin==psi_list_red(k2))
                                p2.x=segments(num(b(ind))).R(1);
                                p2.y=segments(num(b(ind))).Z(1);
                            else
                                p2.x=segments(num(b(ind))).R(end);
                                p2.y=segments(num(b(ind))).Z(end);
                            end
                            if(segments(num(b(ind))).nX~=k)
                                p2.theta=segments(num(b(ind))).nX;
                                p2.branch=segments(num(b(ind))).nB;
                            else %modif
                                p2.theta=segments(num(b(ind))).nX;
                                p2.branch=-1;
                            end
                            p2.type=segments(num(b(ind))).type;
                            segpa=[segpa,num(b(ind))];
                        end
                        if(sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)==0)
                            ind=ind+1;
                        else
                            
                            %     if(ncote==47)
                            %         disp('oc')
                            %     end
                            ncote=ncote+1;
                            cote(ncote).p1.R=p1.x;
                            cote(ncote).p1.Z=p1.y;
                            cote(ncote).p2.R=p2.x;
                            cote(ncote).p2.Z=p2.y;
                            cout=part_contour_per(cprog,p1,p2,sens);
                            cote(ncote).p1.coord=[p1.psi,p1.theta,p1.branch,p1.type];
                            cote(ncote).p2.coord=[p2.psi,p2.theta,p2.branch,p2.type];
                            if(ind>length(b))
                                %remove and close
                                dbr=zeros(1,4);
                                for nb2=1:4
                                    if(nb2~=nb)
                                        dbr(nb2)=min(sqrt((cout.x-X_points.branch(k).R(nb2)).^2+...
                                            (cout.y-X_points.branch(k).Z(nb2)).^2));
                                    else
                                        dbr(nb2)=1e10;
                                    end
                                end
                                p2.x=X_points.branch(k).R(find(dbr==min(dbr)));
                                p2.y=X_points.branch(k).Z(find(dbr==min(dbr)));
                                dist2=sqrt((cout.x-p2.x).^2+(cout.y-p2.y).^2);
                                a=find(dist2==min(dist2));
                                cout.x=[cout.x(1:a-1),p2.x];
                                cout.y=[cout.y(1:a-1),p2.y];
                            end
                            if(sens==1)
                                plot(cout.x,cout.y,'g-','LineWidth',3);
                            else
                                plot(cout.x,cout.y,'m-','LineWidth',3);
                            end
                            cote(ncote).R=cout.x;
                            cote(ncote).Z=cout.y;
                            ind=ind+1;
                            if(ind>length(b))
                                stopit=true;
                            else
                                if((type(b(ind))==2)||(sum(find(segpa==num(b(ind))))~=0)) %lim
                                    stopit=true;
                                end
                                if((p2.type==1)&&(p2.theta==k))
                                    stopit=true;
                                end
                            end
                        end
                    end
                    
                end
                
            else % end of special treatment
                %                 if(k==1&&n==4)
                %                     disp('f')
                %                 end
                c2=contour_better(r2D,z2D,flux2D,psi_list_red(k2));
                X=intersect_contour(c1,c2); %should intersect once and just once
                cprog.x=c2.arc(X.arc2(1)).x;
                cprog.y=c2.arc(X.arc2(1)).y;
                number=X.arc2(1);
                %all contour close normally if extrapolation correct
                %look left
                d=sqrt((cprog.x-X.x(1)).^2+(cprog.y-X.y(1)).^2);
                a=find(d==min(d));
                if(a<length(d))
                    v1x=cprog.x(a+1)-cprog.x(a);
                    v1y=cprog.y(a+1)-cprog.y(a);
                    vx=X.x(1)-cprog.x(a);
                    vy=X.y(1)-cprog.y(a);
                    if(v1x*vx+v1y*vy>0)
                        cleft.x=circshift(cprog.x(1:end-1)',-a)';
                        cleft.y=circshift(cprog.y(1:end-1)',-a)';
                    else
                        cleft.x=circshift(cprog.x(1:end-1)',-a+1)';
                        cleft.y=circshift(cprog.y(1:end-1)',-a+1)';
                    end
                    cleft.x=[cleft.x,cleft.x(1)];
                    cleft.y=[cleft.y,cleft.y(1)];
                else
                    cleft.x=cprog.x;
                    cleft.y=cprog.y;
                end
                dist=[];
                type=[];
                fine=[];
                num=[];
                for ns=1:Nseg
                    if(segments(ns).psimin==psi_list_red(k2))
                        if(segments(ns).psimin_numc==number)
                            X.num=1;
                            X.x(1)=segments(ns).R(1);
                            X.y(1)=segments(ns).Z(1);
                        else
                            X.num=0;
                        end
                    else
                        if(segments(ns).psimax==psi_list_red(k2))
                            if(segments(ns).psimax_numc==number)
                                X.num=1;
                                X.x(1)=segments(ns).R(end);
                                X.y(1)=segments(ns).Z(end);
                            else
                                X.num=0;
                            end
                        else
                            X.num=0;
                        end
                    end
                    if(X.num==1) %intersect
                        d=sqrt((cleft.x-X.x(1)).^2+(cleft.y-X.y(1)).^2);
                        a1=find(d==min(d));
                        cleftp.x=circshift(cleft.x(1:end-1)',-1)';
                        cleftp.y=circshift(cleft.y(1:end-1)',-1)';
                        cleftp.x=[cleftp.x,cleftp.x(1)];
                        cleftp.y=[cleftp.y,cleftp.y(1)];
                        dist=[dist,a1];
                        num=[num,ns*ones(size(a1))];
                        if(ns<=NsegX)
                            type=[type,1*ones(size(a1))]; % Xpoint
                        else
                            type=[type,2*ones(size(a1))]; % lim
                        end
                        vx=X.x(1)-cleft.x(a1);
                        vy=X.y(1)-cleft.y(a1);
                        vx1=cleftp.x(a1)-cleft.x(a1);
                        vy1=cleftp.y(a1)-cleft.y(a1);
                        fine=[fine,min(d)*sign(vx.*vx1+vy.*vy1)];
                    end
                end %end for nsegment
                pasmax=max(abs(fine))+1;
                [a,b]=sort(dist*pasmax+fine);
                stopit=false;
                ind=1;
                segpa=[];
                if(num(b(1))~=num(b(end)))
                    b=[b(end),b];
                end
                if(length(num)==1)
                    num=[num num];
                    b=[1 2];
                    type=[type type];
                end
                if(ind+1>length(b))
                    stopit=true;
                end
                while(~stopit)
                    p1.psi=psi_list_red(k2);
                    p2.psi=psi_list_red(k2);
                    if(segments(num(b(ind))).psimin==psi_list_red(k2))
                        p1.x=segments(num(b(ind))).R(1);
                        p1.y=segments(num(b(ind))).Z(1);
                    else
                        p1.x=segments(num(b(ind))).R(end);
                        p1.y=segments(num(b(ind))).Z(end);
                    end
                    p1.theta=segments(num(b(ind))).nX;
                    p1.branch=segments(num(b(ind))).nB;
                    p1.type=segments(num(b(ind))).type;
                    if(segments(num(b(ind+1))).psimin==psi_list_red(k2))
                        p2.x=segments(num(b(ind+1))).R(1);
                        p2.y=segments(num(b(ind+1))).Z(1);
                    else
                        p2.x=segments(num(b(ind+1))).R(end);
                        p2.y=segments(num(b(ind+1))).Z(end);
                    end
                    p2.theta=segments(num(b(ind+1))).nX;
                    p2.branch=segments(num(b(ind+1))).nB;
                    p2.type=segments(num(b(ind+1))).type;
                    if((sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)==0)&&(psi_list_red(k2)~=psicore))
                        %                     if(sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)==0)
                        ind=ind+1;
                    else
                        segpa=[segpa,num(b(ind))];
                        if(sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)<1e-6)
                            cout.x=[p1.x,cleft.x(1:end-1),p2.x];
                            cout.y=[p1.y,cleft.y(1:end-1),p2.y];
                        else
                            cout=part_contour_per(cprog,p1,p2,1);
                        end
                        OK=true;
                        for k3=1:X_points.num
                            for n3=1:4
                                c11.num=1;
                                c11.arc(1).x=[X_points.R(k3) X_points.branch(k3).R(n3)];
                                c11.arc(1).y=[X_points.Z(k3) X_points.branch(k3).Z(n3)];
                                c21.num=1;
                                c21.arc(1).x=cout.x;
                                c21.arc(1).y=cout.y;
                                X=intersect_contour(c11,c21);
                                if(X.num>0)
                                    OK=false
                                end
                            end
                        end
                        if(OK)
                            
                            %     if(ncote==47)
                            %         disp('oc')
                            %     end
                            ncote=ncote+1;
                            cote(ncote).p1.R=p1.x;
                            cote(ncote).p1.Z=p1.y;
                            cote(ncote).p2.R=p2.x;
                            cote(ncote).p2.Z=p2.y;
                            cote(ncote).p1.coord=[p1.psi,p1.theta,p1.branch,p1.type];
                            cote(ncote).p2.coord=[p2.psi,p2.theta,p2.branch,p2.type];
                            plot(cout.x,cout.y,'g-','LineWidth',3);
                            cote(ncote).R=cout.x;
                            cote(ncote).Z=cout.y;
                        end
                        ind=ind+1;
                    end
                    if((type(b(ind))==2)||(sum(find(segpa==num(b(ind))))~=0)||(ind+1>length(b))) %lim
                        stopit=true;
                    end
%                     if((p2.type==1)&&(p2.theta==k))
%                                     stopit=true;
%                                 end
                end
                %------------------------------------------------------------
                %look right
                cright.x=cleft.x(end:-1:1);
                cright.y=cleft.y(end:-1:1);
                dist=[];
                type=[];
                fine=[];
                num=[];
                for ns=1:Nseg
                    if(segments(ns).psimin==psi_list_red(k2))
                        if(segments(ns).psimin_numc==number)
                            X.num=1;
                            X.x(1)=segments(ns).R(1);
                            X.y(1)=segments(ns).Z(1);
                        else
                            X.num=0;
                        end
                    else
                        if(segments(ns).psimax==psi_list_red(k2))
                            if(segments(ns).psimax_numc==number)
                                X.num=1;
                                X.x(1)=segments(ns).R(end);
                                X.y(1)=segments(ns).Z(end);
                            else
                                X.num=0;
                            end
                        else
                            X.num=0;
                        end
                    end
                    if(X.num==1) %intersect
                        d=sqrt((cright.x-X.x(1)).^2+(cright.y-X.y(1)).^2);
                        a1=find(d==min(d));
                        crightp.x=circshift(cright.x(1:end-1)',-1)';
                        crightp.y=circshift(cright.y(1:end-1)',-1)';
                        crightp.x=[crightp.x,crightp.x(1)];
                        crightp.y=[crightp.y,crightp.y(1)];
                        dist=[dist,a1];
                        num=[num,ns*ones(size(a1))];
                        if(ns<=NsegX)
                            type=[type,1*ones(size(a1))]; % Xpoint
                        else
                            type=[type,2*ones(size(a1))]; % lim
                        end
                        vx=X.x(1)-cright.x(a1);
                        vy=X.y(1)-cright.y(a1);
                        vx1=crightp.x(a1)-cright.x(a1);
                        vy1=crightp.y(a1)-cright.y(a1);
                        fine=[fine,min(d)*sign(vx.*vx1+vy.*vy1)];
                    end
                end %end for nsegment
                pasmax=max(abs(fine))+1;
                [a,b]=sort(dist*pasmax+fine);
                if(num(b(1))~=num(b(end)))
                    b=[b(end),b];
                end
                if(length(num)==1)
                    num=[num num];
                    b=[1 2];
                    type=[type type];
                end
                stopit=false;
                ind=1;
                segpa=[];
                if(ind+1>length(b))
                    if(ind~=1)
                        stopit=true;
                    end
                end
                while(~stopit)
                    p1.psi=psi_list_red(k2);
                    p2.psi=psi_list_red(k2);
                    if(segments(num(b(ind))).psimin==psi_list_red(k2))
                        p1.x=segments(num(b(ind))).R(1);
                        p1.y=segments(num(b(ind))).Z(1);
                    else
                        p1.x=segments(num(b(ind))).R(end);
                        p1.y=segments(num(b(ind))).Z(end);
                    end
                    p1.theta=segments(num(b(ind))).nX;
                    p1.branch=segments(num(b(ind))).nB;
                    p1.type=segments(num(b(ind))).type;
                    if(segments(num(b(ind+1))).psimin==psi_list_red(k2))
                        p2.x=segments(num(b(ind+1))).R(1);
                        p2.y=segments(num(b(ind+1))).Z(1);
                    else
                        p2.x=segments(num(b(ind+1))).R(end);
                        p2.y=segments(num(b(ind+1))).Z(end);
                    end
                    p2.theta=segments(num(b(ind+1))).nX;
                    p2.branch=segments(num(b(ind+1))).nB;
                    p2.type=segments(num(b(ind+1))).type;
                    if((sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)==0)&&(psi_list_red(k2)~=psicore))
                        %                        if(sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)==0)
                        ind=ind+1;
                    else
                        segpa=[segpa,num(b(ind))];
                        if(sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)<1e-6)
                            cout.x=[p1.x,cright.x(2:end),p2.x];
                            cout.y=[p1.y,cright.y(2:end),p2.y];
                        else
                            cout=part_contour_per(cprog,p1,p2,-1);
                        end
                        OK=true;
                        for k3=1:X_points.num
                            for n3=1:4
                                c11.num=1;
                                c11.arc(1).x=[X_points.R(k3) X_points.branch(k3).R(n3)];
                                c11.arc(1).y=[X_points.Z(k3) X_points.branch(k3).Z(n3)];
                                c21.num=1;
                                c21.arc(1).x=cout.x;
                                c21.arc(1).y=cout.y;
                                X=intersect_contour(c11,c21);
                                if(X.num>0)
                                    OK=false
                                end
                            end
                        end
                        if(OK)
                            %               if(ncote==47)
                            %         disp('oc')
                            %     end
                            ncote=ncote+1;
                            cote(ncote).p1.R=p1.x;
                            cote(ncote).p1.Z=p1.y;
                            cote(ncote).p2.R=p2.x;
                            cote(ncote).p2.Z=p2.y;
                            cote(ncote).p1.coord=[p1.psi,p1.theta,p1.branch,p1.type];
                            cote(ncote).p2.coord=[p2.psi,p2.theta,p2.branch,p2.type];
                            plot(cout.x,cout.y,'m-','LineWidth',3);
                            cote(ncote).R=cout.x;
                            cote(ncote).Z=cout.y;
                        end
                        ind=ind+1;
                    end
                    if((type(b(ind))==2)||(sum(find(segpa==num(b(ind))))~=0)||(ind+1>length(b))) %lim
                        stopit=true;
                    end
%                     if((p2.type==1)&&(p2.theta==k))
%                                     stopit=true;
%                                 end
                end
                
            end % end if not separatrice
        end
        
    end
    
end
