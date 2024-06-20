n0=  1.6666667E+20;
 T0=   135.373690640864     ;
 c0=   80537.4225381667     ;
 tau0=  1.560314474728064E-004;
 R0=2;

Rperp_in=linspace(221.44,180.015,100)/100;
Zperp_in=linspace(0,0,100);
Rsep1_in=211.51/100;
Rsep2_in=208.3675/100;

figure
nperp_in=griddata(R/100,Z/100,density,Rperp_in,Zperp_in);
plot(Rperp_in,nperp_in*n0,'b.-')
hold on
plot([Rsep1_in Rsep1_in],[0 max(nperp_in*n0)*1.1],'k--')
plot([Rsep2_in Rsep2_in],[0 max(nperp_in*n0)*1.1],'r--')
title('Inner density (m^{-3})')
ylim([0 max(nperp_in)*n0*1.1])
xlabel('radial coordinate (m)')

figure
Mperp_in=griddata(R/100,Z/100,M,Rperp_in,Zperp_in);
plot(Rperp_in,Mperp_in,'b.-')
hold on
plot([Rsep1_in Rsep1_in],[1.1*min(Mperp_in) max(Mperp_in)*1.1],'k--')
plot([Rsep2_in Rsep2_in],[1.1*min(Mperp_in) max(Mperp_in)*1.1],'r--')
title('Inner Mach number')
ylim([1.1*min(Mperp_in) max(Mperp_in)*1.1]);
xlabel('radial coordinate (m)')

figure
Teperp_in=griddata(R/100,Z/100,Te,Rperp_in,Zperp_in);
hold on
plot(Rperp_in,Teperp_in*T0,'b.-')
Tiperp_in=griddata(R/100,Z/100,Ti,Rperp_in,Zperp_in);
plot(Rperp_in,Tiperp_in*T0,'r.-')
legend('Te','Ti')
hold on
plot([Rsep1_in Rsep1_in],[0 max(Tiperp_in*T0)*1.1],'k--')
plot([Rsep2_in Rsep2_in],[0 max(Tiperp_in*T0)*1.1],'r--')
ylim([0 max(Tiperp_in)*1.1*T0]);
title('Inner Temperatures (eV)')
xlabel('radial coordinate (m)')


figure
SNperp_in=griddata(R/100,Z/100,SN*n0/tau0,Rperp_in,Zperp_in);
plot(Rperp_in,SNperp_in,'b.-')
hold on
plot([Rsep1_in Rsep1_in],[min(SNperp_in)*1.1 max(SNperp_in)*1.1],'k--')
plot([Rsep2_in Rsep2_in],[min(SNperp_in)*1.1 max(SNperp_in)*1.1],'r--')
plot([min(Rperp_in) max(Rperp_in)],[0 0],'k-')
ylim([min(SNperp_in)*1.1 max(SNperp_in)*1.1])
title('Inner Ionization source (m^{-3}s^{-1})')
xlabel('radial coordinate (m)')

figure
PNperp_in=griddata(R/100,Z/100,PN,Rperp_in,Zperp_in);
plot(Rperp_in,PNperp_in*n0,'b.-')
hold on
plot([Rsep1_in Rsep1_in],[0 max(PNperp_in*n0)*1.1],'k--')
plot([Rsep2_in Rsep2_in],[0 max(PNperp_in*n0)*1.1],'r--')
ylim([0 max(PNperp_in*n0)*1.1])
title('Inner Neutral pressure (Pa)')
xlabel('radial coordinate (m)')