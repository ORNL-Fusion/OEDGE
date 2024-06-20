Nseg=zone_elements.Nseg; % total number of segments
NsegX=zone_elements.Nseg-zone_elements.Nseg2; % segments from X points (two sided)
Nseglim=zone_elements.Nseg2; % segments on limits (one sided)

SegCnt=zeros(1,Nseg);

clear segments

for k=1:X_points.num
    for n=1:4
        for ns=1:X_points.cut(k).arc(n).nPt_lim-1
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).psimin=...
                min([X_points.cut(k).arc(n).Pt_lim(ns).psi X_points.cut(k).arc(n).Pt_lim(ns+1).psi]);
            wh=find([X_points.cut(k).arc(n).Pt_lim(ns).psi X_points.cut(k).arc(n).Pt_lim(ns+1).psi]==...
                min([X_points.cut(k).arc(n).Pt_lim(ns).psi X_points.cut(k).arc(n).Pt_lim(ns+1).psi]));
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).psimin_numc=X_points.cut(k).arc(n).Pt_lim(ns-1+wh).psi_numc;
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).psimax=...
                max([X_points.cut(k).arc(n).Pt_lim(ns).psi X_points.cut(k).arc(n).Pt_lim(ns+1).psi]);
            wh=find([X_points.cut(k).arc(n).Pt_lim(ns).psi X_points.cut(k).arc(n).Pt_lim(ns+1).psi]==...
                max([X_points.cut(k).arc(n).Pt_lim(ns).psi X_points.cut(k).arc(n).Pt_lim(ns+1).psi]));
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).psimax_numc=X_points.cut(k).arc(n).Pt_lim(ns-1+wh).psi_numc;
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).type=1;
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).Xtype=X_points.cut(k).arc(n).type;
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).nX=k;
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).nB=n;
            psi1=interp2(r2D,z2D,flux2D,X_points.cut(k).arc(n).seg_lim(ns).R(1),X_points.cut(k).arc(n).seg_lim(ns).Z(1));
            psi2=interp2(r2D,z2D,flux2D,X_points.cut(k).arc(n).seg_lim(ns).R(end),X_points.cut(k).arc(n).seg_lim(ns).Z(end));
            if(psi1<psi2)
                segments(X_points.cut(k).arc(n).seg_lim(ns).num).R=X_points.cut(k).arc(n).seg_lim(ns).R;
                segments(X_points.cut(k).arc(n).seg_lim(ns).num).Z=X_points.cut(k).arc(n).seg_lim(ns).Z;
            else
                segments(X_points.cut(k).arc(n).seg_lim(ns).num).R=X_points.cut(k).arc(n).seg_lim(ns).R(end:-1:1);
                segments(X_points.cut(k).arc(n).seg_lim(ns).num).Z=X_points.cut(k).arc(n).seg_lim(ns).Z(end:-1:1);
            end
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).p1.coord=[psi1,k,segments(X_points.cut(k).arc(n).seg_lim(ns).num).nB,1];
            segments(X_points.cut(k).arc(n).seg_lim(ns).num).p2.coord=[psi2,k,segments(X_points.cut(k).arc(n).seg_lim(ns).num).nB,1]
        end
    end
end
for k=1:frontiers.num
    for ns=1:frontiers.lim(k).nPt_lim-1
        segments(frontiers.lim(k).seg_lim(ns).num).psimin=...
            min([frontiers.lim(k).Pt_lim(ns).psi frontiers.lim(k).Pt_lim(ns+1).psi]);
        wh=find([frontiers.lim(k).Pt_lim(ns).psi frontiers.lim(k).Pt_lim(ns+1).psi]==...
                min([frontiers.lim(k).Pt_lim(ns).psi frontiers.lim(k).Pt_lim(ns+1).psi]));
        segments(frontiers.lim(k).seg_lim(ns).num).psimin_numc=frontiers.lim(k).Pt_lim(ns-1+wh).psi_numc;
        segments(frontiers.lim(k).seg_lim(ns).num).psimax=...
            max([frontiers.lim(k).Pt_lim(ns).psi frontiers.lim(k).Pt_lim(ns+1).psi]);
        wh=find([frontiers.lim(k).Pt_lim(ns).psi frontiers.lim(k).Pt_lim(ns+1).psi]==...
                max([frontiers.lim(k).Pt_lim(ns).psi frontiers.lim(k).Pt_lim(ns+1).psi]));
        segments(frontiers.lim(k).seg_lim(ns).num).psimax_numc=frontiers.lim(k).Pt_lim(ns-1+wh).psi_numc;
        segments(frontiers.lim(k).seg_lim(ns).num).type=2;
        segments(frontiers.lim(k).seg_lim(ns).num).Xtype=0;
        segments(frontiers.lim(k).seg_lim(ns).num).nX=k;
        segments(frontiers.lim(k).seg_lim(ns).num).nB=0;
        segments(frontiers.lim(k).seg_lim(ns).num).R=frontiers.lim(k).seg_lim(ns).R;
        segments(frontiers.lim(k).seg_lim(ns).num).Z=frontiers.lim(k).seg_lim(ns).Z;
        segments(frontiers.lim(k).seg_lim(ns).num).p1.coord=[segments(frontiers.lim(k).seg_lim(ns).num).psimin,k,segments(frontiers.lim(k).seg_lim(ns).num).nB,2];
        segments(frontiers.lim(k).seg_lim(ns).num).p2.coord=[segments(frontiers.lim(k).seg_lim(ns).num).psimax,k,segments(frontiers.lim(k).seg_lim(ns).num).nB,2]
    end
end