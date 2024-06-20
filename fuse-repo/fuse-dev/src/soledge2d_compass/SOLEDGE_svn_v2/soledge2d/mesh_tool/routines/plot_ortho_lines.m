
h =waitbar(0,'Please wait: Xpoint number 1');

Nsteps=4*X_points.num;

for k=1:X_points.num
   
    for n=1:4
       
        %start point
        R_=X_points.cut(k).R(n);
        Z_=X_points.cut(k).Z(n);
        psi_=X_points.cut(k).psi(n);
        psiX=X_points.psi(k);
        if(psi_<psiX) %psi descendant --> psimin=psicore
            Rv=[X_points.R(k), R_];
            Zv=[X_points.Z(k), Z_];
        else %psi montant --> psimax = psiout
            Rv=[X_points.R(k), R_];
            Zv=[X_points.Z(k), Z_];
        end
        
        X_points.cut(k).arc(n).R=Rv;
        X_points.cut(k).arc(n).Z=Zv;
        X_points.cut(k).arc(n).type=1;
    end
    
    psilim=zeros(1,4);
    psilim(1)=X_points.cut(k).psilim(1);
    psilim(2)=psiout;
    psilim(3)=X_points.cut(k).psilim(3);
    psilim(4)=psiout;
    for n=1:4
        [R,Z]=follow_grad(X_points.cut(k).arc(n).R(end),X_points.cut(k).arc(n).Z(end),psilim(n),r2D,z2D,flux2D);
        X_points.cut(k).arc(n).R=[X_points.cut(k).arc(n).R,R(2:end)];
        X_points.cut(k).arc(n).Z=[X_points.cut(k).arc(n).Z,Z(2:end)];
        Nstep=(k-1)*4+n;
        waitbar(Nstep/Nsteps,h,['Please wait: Xpoint number ',num2str(k)])
    end
end 

delete(h)