n0=  1.6666667E+20;
 T0=   135.373690640864     ;
 c0=   80537.4225381667     ;
 tau0=  1.560314474728064E-004;
 R0=2;

nl=1;

load ps_files/Metrics/G_002.txt
load ps_files/Metrics/G_005.txt
load ps_files/Metrics/G_012.txt
load ps_files/Metrics/G_018.txt
load ps_files/Metrics/G_015.txt


G=[G_002(nl+1,2:end-1),G_005(nl+1,2:end-1),G_012(nl+1,2:end-1)...
    ,G_018(nl+1,2:end-1),G_015(nl+1,2:end-1)];


[m,p]=size(G);

load ps_files/Mesh/chi_002.txt
load ps_files/Mesh/chi_005.txt
load ps_files/Mesh/chi_012.txt
load ps_files/Mesh/chi_018.txt
load ps_files/Mesh/chi_015.txt
load ps_files/Mesh/meshzb_002.txt
load ps_files/Mesh/meshzb_005.txt
load ps_files/Mesh/meshzb_012.txt
load ps_files/Mesh/meshzb_018.txt
load ps_files/Mesh/meshzb_015.txt
chi=[chi_002,chi_005,chi_012,chi_018,chi_015];

Pts=find(chi(nl,:)==0);

B=zeros(1,p);
for i=1:p
    if(chi(nl,i)==1)
        B(i)=NaN;
    end
end
dz=[meshzb_002(1,2:end)-meshzb_002(1,1:end-1),...
    meshzb_005(1,2:end)-meshzb_005(1,1:end-1),...
    meshzb_012(1,2:end)-meshzb_012(1,1:end-1),...
    meshzb_018(1,2:end)-meshzb_018(1,1:end-1),...
    meshzb_015(1,2:end)-meshzb_015(1,1:end-1)];

s=0;
for i=1:p-1
    s=[s,s(i)+2/(G(i)+G(i+1))*(dz(i))];
end
s=s*2*pi*R0;

load ps_files/Magnetic_input/R_002.txt
load ps_files/Magnetic_input/R_005.txt
load ps_files/Magnetic_input/R_012.txt
load ps_files/Magnetic_input/R_018.txt
load ps_files/Magnetic_input/R_015.txt

Rsep=[R_002(nl,:),R_005(nl,:),R_012(nl,:),R_018(nl,:),R_015(nl,:)];

load ps_files/Magnetic_input/Z_002.txt
load ps_files/Magnetic_input/Z_005.txt
load ps_files/Magnetic_input/Z_012.txt
load ps_files/Magnetic_input/Z_018.txt
load ps_files/Magnetic_input/Z_015.txt

Zsep=[Z_002(nl,:),Z_005(nl,:),Z_012(nl,:),Z_018(nl,:),Z_015(nl,:)];

Rsep=Rsep(Pts);
Zsep=Zsep(Pts);

s=s(Pts);
s=s-s(1);

figure
nsep=griddata(R/100,Z/100,density,Rsep,Zsep);
plot(s,nsep*n0,'b.-')
title('density (m^{-3})')
xlabel('connection length (m)')

figure
Msep=griddata(R/100,Z/100,M,Rsep,Zsep);
plot(s,Msep,'b.-')
title('Mach number')
xlabel('connection length (m)')

figure
Tesep=griddata(R/100,Z/100,Te,Rsep,Zsep);
hold on
plot(s,Tesep*T0,'b.-')
Tisep=griddata(R/100,Z/100,Ti,Rsep,Zsep);
plot(s,Tisep*T0,'r.-')
legend('Te','Ti')
title('Temperatures (eV)')
xlabel('connection length (m)')

figure
Gsep=griddata(R/100,Z/100,Gamma,Rsep,Zsep);
Pisep=nsep.*(Tesep+Tisep)+Gsep.^2./nsep;
plot(s,Pisep*T0*1.6e-19*n0,'b.-')
title('Total plasma pressure (Pa)')
xlabel('connection length (m)')

figure
SNsep=griddata(R/100,Z/100,SN*n0/tau0,Rsep,Zsep);
plot(s,SNsep,'b.-')
title('Ionization source (m^{-3}s^{-1})')
xlabel('connection length (m)')

figure
PNsep=griddata(R/100,Z/100,PN,Rsep,Zsep);
plot(s,PNsep*n0,'b.-')
title('Neutral pressure (Pa)')
xlabel('connection length (m)')