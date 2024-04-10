function [R,Z]=follow_grad(Rstart,Zstart,psi_objective,r2D,z2D,flux2D,step)

if(nargin<7)
    step=3e-3;
end

% compute first derivatives
dp_dR=(flux2D(2:end-1,3:end)-flux2D(2:end-1,1:end-2))./(r2D(2:end-1,3:end)-r2D(2:end-1,1:end-2));
dp_dZ=(flux2D(3:end,2:end-1)-flux2D(1:end-2,2:end-1))./(z2D(3:end,2:end-1)-z2D(1:end-2,2:end-1));

rred=r2D(2:end-1,2:end-1);
zred=z2D(2:end-1,2:end-1);


psistart=interp2(r2D,z2D,flux2D,Rstart,Zstart);
if(psistart>psi_objective)
    signe=1;
else
    signe=-1;
end

R=Rstart;
Z=Zstart;
psi=psistart;

while(signe*(psi-psi_objective)>0)
    dpdR_=interp2(rred,zred,dp_dR,R(end),Z(end));
    dpdZ_=interp2(rred,zred,dp_dZ,R(end),Z(end));
    Rnew=R(end)-signe*dpdR_/sqrt(dpdR_^2+dpdZ_^2)*step;
    Znew=Z(end)-signe*dpdZ_/sqrt(dpdR_^2+dpdZ_^2)*step;
    R=[R,Rnew];
    Z=[Z,Znew];
    psiold=psi;
    psi=interp2(r2D,z2D,flux2D,R(end),Z(end));
    if(signe*(psi-psiold)>0)
        break %change sign in decent
    end
end


end