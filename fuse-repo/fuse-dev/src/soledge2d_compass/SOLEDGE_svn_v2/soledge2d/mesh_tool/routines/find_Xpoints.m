function X_points = find_Xpoints(r2D,z2D,flux2D)
clear X_points
X_points.num=0;

% compute first derivatives
dp_dR=(flux2D(2:end-1,3:end)-flux2D(2:end-1,1:end-2))./(r2D(2:end-1,3:end)-r2D(2:end-1,1:end-2));
dp_dZ=(flux2D(3:end,2:end-1)-flux2D(1:end-2,2:end-1))./(z2D(3:end,2:end-1)-z2D(1:end-2,2:end-1));

rred=r2D(2:end-1,2:end-1);
zred=z2D(2:end-1,2:end-1);
d2p_dR2=(dp_dR(2:end-1,3:end)-dp_dR(2:end-1,1:end-2))./(rred(2:end-1,3:end)-rred(2:end-1,1:end-2));
d2p_dZ2=(dp_dZ(3:end,2:end-1)-dp_dZ(1:end-2,2:end-1))./(zred(3:end,2:end-1)-zred(1:end-2,2:end-1));

% plot contour lines where the first derivatives are zero
figure
cx = contour_better(r2D(2:end-1,2:end-1),z2D(2:end-1,2:end-1),dp_dR,[0 0],'r');
hold on
cy = contour_better(r2D(2:end-1,2:end-1),z2D(2:end-1,2:end-1),dp_dZ,[0 0],'g');

x_cxcy=intersect_contour(cx,cy);
for k=1:x_cxcy.num
    d2p_dR2_=interp2(rred(2:end-1,2:end-1),zred(2:end-1,2:end-1),d2p_dR2,x_cxcy.x(k),x_cxcy.y(k));
    d2p_dZ2_=interp2(rred(2:end-1,2:end-1),zred(2:end-1,2:end-1),d2p_dZ2,x_cxcy.x(k),x_cxcy.y(k));
    %threshold for merdiques cases
    th=abs(max(max(flux2D))-min(min(flux2D)))/abs(max(max(r2D))-min(min(r2D))).^2*0.6;
    if((abs(d2p_dR2_)<th)||(abs(d2p_dZ2_)<th)) %manual
        hf=figure(1)
        rguess=x_cxcy.x(k);
        zguess=x_cxcy.y(k);
        psiguess=interp2(r2D,z2D,flux2D,rguess,zguess);
        contour(r2D,z2D,flux2D,[psiguess psiguess])
        hold on
        plot(rguess,zguess,'ks')
        hold off
        choice = questdlg('X or O point?','X or O point',...
            'O','X','X');
        close(hf);
        switch choice
            case 'X'
                X_points.num=X_points.num+1;
                X_points.R(X_points.num)=x_cxcy.x(k);
                X_points.Z(X_points.num)=x_cxcy.y(k);
                X_points.psi(X_points.num)=interp2(r2D,z2D,flux2D,X_points.R(X_points.num),X_points.Z(X_points.num));
                plot(x_cxcy.x(k),x_cxcy.y(k),'bo')
            case 'O'
        end
    else %auto
        if(d2p_dR2_*d2p_dZ2_<0)
            X_points.num=X_points.num+1;
            X_points.R(X_points.num)=x_cxcy.x(k);
            X_points.Z(X_points.num)=x_cxcy.y(k);
            X_points.psi(X_points.num)=interp2(r2D,z2D,flux2D,X_points.R(X_points.num),X_points.Z(X_points.num));
            plot(x_cxcy.x(k),x_cxcy.y(k),'bo')
        end
    end
end

psi=zeros(1,X_points.num);
for k=1:X_points.num
    psi(k)=X_points.psi(k);
end
[a,b]=sort(psi);
for k=1:X_points.num
    X_points.index(k)=find(a==X_points.psi(k));
    X_points.cut(k).psilim=zeros(1,4);
end

end