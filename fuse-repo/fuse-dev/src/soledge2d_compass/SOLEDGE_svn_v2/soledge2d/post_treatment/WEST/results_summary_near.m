close all

n0=  1.6666667E+20;
 T0=   114.015765857941     ;
 c0=   73911.6753502835     ;
 tau0=  1.700187494709652E-004;

load ps_files_near/WEST_paroi2.txt
load ps_files_near/psi_WEST2

cd(case_dir);

cd Results
load R.txt
load Z.txt

load density.txt
load Gamma.txt
load Ti.txt
load Te.txt
load NN.txt
load SN.txt
load TN.txt

PN=NN*1.6e-19.*TN; %  Pa

M=Gamma./(density.*sqrt(Te+Ti));

cd ..

figure;
patch(R',Z',density'*n0);
axis equal tight
shading flat
hold on
plot(WEST_paroi2(:,1)*100,WEST_paroi2(:,2)*100,'k-','LineWidth',2)
contour(r*100,z*100,psi,[psisep1 psisep2],'k--')
title('electron density (m^{-3})')
caxis([min(min(density*n0)) max(max(density*n0))])
colorbar

figure;
patch(R',Z',Gamma');
axis equal tight
shading flat
hold on
plot(WEST_paroi2(:,1)*100,WEST_paroi2(:,2)*100,'k-','LineWidth',2)
contour(r*100,z*100,psi,[psisep1 psisep2],'k--')
title('Gamma')
caxis([min(min(Gamma)) max(max(Gamma))])
colorbar

figure;
patch(R',Z',M');
axis equal tight
shading flat
hold on
plot(WEST_paroi2(:,1)*100,WEST_paroi2(:,2)*100,'k-','LineWidth',2)
contour(r*100,z*100,psi,[psisep1 psisep2],'k--')
title('Mach number')
caxis([-1 1])
colorbar

figure;
patch(R',Z',Te'*T0);
axis equal tight
shading flat
hold on
plot(WEST_paroi2(:,1)*100,WEST_paroi2(:,2)*100,'k-','LineWidth',2)
contour(r*100,z*100,psi,[psisep1 psisep2],'k--')
title('Te (eV)')
caxis([min(min(Te*T0)) max(max(Te*T0))])
colorbar

figure;
patch(R',Z',Ti'*T0);
axis equal tight
shading flat
hold on
plot(WEST_paroi2(:,1)*100,WEST_paroi2(:,2)*100,'k-','LineWidth',2)
contour(r*100,z*100,psi,[psisep1 psisep2],'k--')
title('Ti (eV)')
caxis([min(min(Ti*T0)) max(max(Ti*T0))])
colorbar

figure;
patch(R',Z',log10(NN'*n0));
axis equal tight
shading flat
hold on
plot(WEST_paroi2(:,1)*100,WEST_paroi2(:,2)*100,'k-','LineWidth',2)
contour(r*100,z*100,psi,[psisep1 psisep2],'k--')
title('Neutral density (m^{-3} log)')
caxis([max(max(log10(NN*n0)))-3 max(max(log10(NN*n0)))])
colorbar

figure;
patch(R',Z',log10(abs(PN'*n0)));
axis equal tight
shading flat
hold on
plot(WEST_paroi2(:,1)*100,WEST_paroi2(:,2)*100,'k-','LineWidth',2)
contour(r*100,z*100,psi,[psisep1 psisep2],'k--')
title('Neutral pressure (Pa)')
caxis([max(max(log10(abs(PN*n0))))-3 max(max(log10(abs(PN*n0))))])
% caxis([0 max(max(PN*n0))])
colorbar

figure;
patch(R',Z',log10(SN'*n0/tau0));
axis equal tight
shading flat
hold on
plot(WEST_paroi2(:,1)*100,WEST_paroi2(:,2)*100,'k-','LineWidth',2)
contour(r*100,z*100,psi,[psisep1 psisep2],'k--')
title('Ionization source all (log)')
caxis([max(max(log10(SN*n0/tau0)))-3 max(max(log10(SN*n0/tau0)))])
colorbar

load soledge2D.energy_fluxes
figure
plot(soledge2D(:,4),soledge2D(:,8),'b-')
title('energy fluxes (W/m^2)')

load soledge2D.particle_fluxes_wall
figure
plot(soledge2D(:,4),soledge2D(:,5),'b-')
title('particle fluxes (m^{-2})')

cd ..