n0=  1.6666667E+20;
 T0=   135.373690640864     ;
 c0=   80537.4225381667     ;
 tau0=  1.560314474728064E-004;
 R0=2;

Rperp_out=linspace(286.5,297.7,100)/100;
Zperp_out=linspace(0,0,100);
Rsep1_out=293.32/100;
Rsep2_out=295.548/100;

figure
nperp_out=griddata(R/100,Z/100,density,Rperp_out,Zperp_out);
plot(Rperp_out,nperp_out*n0,'b.-')
hold on
plot([Rsep1_out Rsep1_out],[0 max(nperp_out*n0)*1.1],'k--')
plot([Rsep2_out Rsep2_out],[0 max(nperp_out*n0)*1.1],'r--')
title('Inner density (m^{-3})')
ylim([0 max(nperp_out)*n0*1.1])
xlabel('radial coordinate (m)')

figure
Mperp_out=griddata(R/100,Z/100,M,Rperp_out,Zperp_out);
plot(Rperp_out,Mperp_out,'b.-')
hold on
plot([Rsep1_out Rsep1_out],[1.1*min(Mperp_out) max(Mperp_out)*1.1],'k--')
plot([Rsep2_out Rsep2_out],[1.1*min(Mperp_out) max(Mperp_out)*1.1],'r--')
title('Inner Mach number')
ylim([1.1*min(Mperp_out) max(Mperp_out)*1.1]);
xlabel('radial coordinate (m)')

figure
Teperp_out=griddata(R/100,Z/100,Te,Rperp_out,Zperp_out);
hold on
plot(Rperp_out,Teperp_out*T0,'b.-')
Tiperp_out=griddata(R/100,Z/100,Ti,Rperp_out,Zperp_out);
plot(Rperp_out,Tiperp_out*T0,'r.-')
legend('Te','Ti')
hold on
plot([Rsep1_out Rsep1_out],[0 max(Tiperp_out*T0)*1.1],'k--')
plot([Rsep2_out Rsep2_out],[0 max(Tiperp_out*T0)*1.1],'r--')
ylim([0 max(Tiperp_out)*1.1*T0]);
title('Inner Temperatures (eV)')
xlabel('radial coordinate (m)')


figure
SNperp_out=griddata(R/100,Z/100,SN*n0/tau0,Rperp_out,Zperp_out);
plot(Rperp_out,SNperp_out,'b.-')
hold on
plot([Rsep1_out Rsep1_out],[min(SNperp_out)*1.1 max(SNperp_out)*1.1],'k--')
plot([Rsep2_out Rsep2_out],[min(SNperp_out)*1.1 max(SNperp_out)*1.1],'r--')
plot([min(Rperp_out) max(Rperp_out)],[0 0],'k-')
ylim([min(SNperp_out)*1.1 max(SNperp_out)*1.1])
title('Inner Ionization source (m^{-3}s^{-1})')
xlabel('radial coordinate (m)')

figure
PNperp_out=griddata(R/100,Z/100,PN,Rperp_out,Zperp_out);
plot(Rperp_out,PNperp_out*n0,'b.-')
hold on
plot([Rsep1_out Rsep1_out],[0 max(PNperp_out*n0)*1.1],'k--')
plot([Rsep2_out Rsep2_out],[0 max(PNperp_out*n0)*1.1],'r--')
ylim([0 max(PNperp_out*n0)*1.1])
title('Inner Neutral pressure (Pa)')
xlabel('radial coordinate (m)')