h2 = waitbar(0,'Please wait');

Nsteps=X_points.num*4+frontiers.num;

Npt=X_points.num;
Nseg=0;

for k=1:X_points.num
    
    for n=1:4
        %for each arc, find intersections of interest
        c1.num=1;
        c1.arc(1).x=X_points.cut(k).arc(n).R;
        c1.arc(1).y=X_points.cut(k).arc(n).Z;
        psi_list=zeros(1,1+X_points.num);
        psi_list(1)=X_points.cut(k).psilim(n);
        psi_list(2)=X_points.psi(k);
        ind=2;
        for k1=1:X_points.num
            if(k1~=k)
                if((X_points.psi(k1)-X_points.psi(k))*(X_points.psi(k1)-X_points.cut(k).psilim(n))<0)
                    ind=ind+1;
                    psi_list(ind)=X_points.psi(k1);
                end
            end
        end
        psi_list_red=psi_list(1:ind);
        psi_list_red=sort(psi_list_red);
        for k1=1:length(psi_list_red)
            if(psi_list_red(k1)~=X_points.psi(k))
                c2=contour_better(r2D,z2D,flux2D,psi_list_red(k1));
                X=intersect_contour(c1,c2);
                X_points.cut(k).arc(n).Pt_lim(k1).R=X.x(1);
                X_points.cut(k).arc(n).Pt_lim(k1).Z=X.y(1);
                X_points.cut(k).arc(n).Pt_lim(k1).psi=psi_list_red(k1);
                X_points.cut(k).arc(n).Pt_lim(k1).psi_numc=X.arc2(1);
                Npt=Npt+1;
                X_points.cut(k).arc(n).Pt_lim(k1).num=Npt;
            else
                X_points.cut(k).arc(n).Pt_lim(k1).R=X_points.R(k);
                X_points.cut(k).arc(n).Pt_lim(k1).Z=X_points.Z(k);
                X_points.cut(k).arc(n).Pt_lim(k1).psi=psi_list_red(k1);
                c2=contour_better(r2D,z2D,flux2D,psi_list_red(k1));
                p1.x=X_points.cut(k).arc(n).R(1);
                p1.y=X_points.cut(k).arc(n).Z(1);
                X_points.cut(k).arc(n).Pt_lim(k1).psi_numc=which_contour_num(c2,p1);
                X_points.cut(k).arc(n).Pt_lim(k1).num=k;
            end
        end
        for k1=1:length(psi_list_red)-1
            p1.x=X_points.cut(k).arc(n).Pt_lim(k1).R;
            p1.y=X_points.cut(k).arc(n).Pt_lim(k1).Z;
            p2.x=X_points.cut(k).arc(n).Pt_lim(k1+1).R;
            p2.y=X_points.cut(k).arc(n).Pt_lim(k1+1).Z;
            cin.x=c1.arc(1).x;
            cin.y=c1.arc(1).y;
            cout=part_contour2(cin,p1,p2);
            X_points.cut(k).arc(n).seg_lim(k1).R=cout.x;
            X_points.cut(k).arc(n).seg_lim(k1).Z=cout.y;
            Nseg=Nseg+1;
            X_points.cut(k).arc(n).seg_lim(k1).num=Nseg;
        end
        X_points.cut(k).arc(n).nPt_lim=length(psi_list_red);
        Nstep=(k-1)*4+n;
        waitbar(Nstep/Nsteps)
        
    end
    
end

Nseg2=0;

for k=1:frontiers.num
    
    c1.num=1;
    c1.arc(1).x=frontiers.lim(k).R';
    c1.arc(1).y=frontiers.lim(k).Z';
    psi_list=zeros(1,2+X_points.num);
    psi_list(1)=frontiers.lim(k).psimin;
    psi_list(2)=frontiers.lim(k).psimax;
    ind=2;
    for k1=1:X_points.num
        if((X_points.psi(k1)>frontiers.lim(k).psimin)&&(X_points.psi(k1)<frontiers.lim(k).psimax))
            ind=ind+1;
            psi_list(ind)=X_points.psi(k1);
        end
    end
    psi_list_red=psi_list(1:ind);
    psi_list_red=sort(psi_list_red);
    for k1=1:length(psi_list_red)
            c2=contour_better(r2D,z2D,flux2D,psi_list_red(k1));
            X=intersect_contour(c1,c2);
            frontiers.lim(k).Pt_lim(k1).R=X.x(1);
            frontiers.lim(k).Pt_lim(k1).Z=X.y(1);
            frontiers.lim(k).Pt_lim(k1).psi=psi_list_red(k1);
            frontiers.lim(k).Pt_lim(k1).psi_numc=X.arc2(1);
            Npt=Npt+1;
            frontiers.lim(k).Pt_lim(k1).num=Npt;
    end
    for k1=1:length(psi_list_red)-1
        p1.x=frontiers.lim(k).Pt_lim(k1).R;
        p1.y=frontiers.lim(k).Pt_lim(k1).Z;
        p2.x=frontiers.lim(k).Pt_lim(k1+1).R;
        p2.y=frontiers.lim(k).Pt_lim(k1+1).Z;
        cin.x=c1.arc(1).x;
        cin.y=c1.arc(1).y;
        cout=part_contour2(cin,p1,p2);
        frontiers.lim(k).seg_lim(k1).R=cout.x;
        frontiers.lim(k).seg_lim(k1).Z=cout.y;
        Nseg=Nseg+1;
        Nseg2=Nseg2+1;
        frontiers.lim(k).seg_lim(k1).num=Nseg;
    end
    frontiers.lim(k).nPt_lim=length(psi_list_red);
    Nstep=k+X_points.num*4;
    waitbar(Nstep/Nsteps)
    
end

zone_elements.Nseg=Nseg;
zone_elements.Nseg2=Nseg2;
zone_elements.Npt=Npt;

close(h2)