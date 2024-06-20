global raise_value;

fluxmin=min(min(flux2D));
fluxmax=max(max(flux2D));

Rmin=min(min(r2D));
Zmin=min(min(z2D));
Rmax=max(max(r2D));
Zmax=max(max(z2D));

dpsi=(fluxmax-fluxmin)*raise_value;

%50/200 ==> +25% de flux

% extrapolate radial 
for i=1:50
    r2D=[r2D,Rmax*ones(size(r2D(:,end)))+i/50*((Rmax-Rmin)*extrapol_val)];
    z2D=[z2D,z2D(:,end)];
    flux2D=[flux2D,2*flux2D(:,end)-flux2D(:,end-1)+dpsi*i];
    Br2D=[Br2D,2*Br2D(:,end)-Br2D(:,end-1)];
    Bz2D=[Bz2D,2*Bz2D(:,end)-Bz2D(:,end-1)];
    Bphi2D=[Bphi2D,2*Bphi2D(:,end)-Bphi2D(:,end-1)];
end
for i=1:50
    r2D=[Rmin*ones(size(r2D(:,1)))-i/50*((Rmax-Rmin)*extrapol_val),r2D];
    z2D=[z2D(:,1),z2D];
    flux2D=[2*flux2D(:,1)-flux2D(:,2)+dpsi*i,flux2D];
    Br2D=[2*Br2D(:,1)-Br2D(:,2),Br2D];
    Bz2D=[2*Bz2D(:,1)-Bz2D(:,2),Bz2D];
    Bphi2D=[2*Bphi2D(:,1)-Bphi2D(:,2),Bphi2D];
end


% extrapolate vertical
for i=1:50
    r2D=[r2D;r2D(end,:)];
    z2D=[z2D;Zmax*ones(size(z2D(end,:)))+i/50*((Zmax-Zmin)*extrapol_val)];
    flux2D=[flux2D;2*flux2D(end,:)-flux2D(end-1,:)+dpsi*i];
    Br2D=[Br2D;2*Br2D(end,:)-Br2D(end-1,:)];
    Bz2D=[Bz2D;2*Bz2D(end,:)-Bz2D(end-1,:)];
    Bphi2D=[Bphi2D;2*Bphi2D(end,:)-Bphi2D(end-1,:)];
end
for i=1:50
    r2D=[r2D(1,:);r2D];
    z2D=[Zmin*ones(size(z2D(1,:)))-i/50*((Zmax-Zmin)*extrapol_val);z2D];
    flux2D=[2*flux2D(1,:)-flux2D(2,:)+dpsi*i;flux2D];
    Br2D=[2*Br2D(1,:)-Br2D(2,:);Br2D];
    Bz2D=[2*Bz2D(1,:)-Bz2D(2,:);Bz2D];
    Bphi2D=[2*Bphi2D(1,:)-Bphi2D(2,:);Bphi2D];
end
