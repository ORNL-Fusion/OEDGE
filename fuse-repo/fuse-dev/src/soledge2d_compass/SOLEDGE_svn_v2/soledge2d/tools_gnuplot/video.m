load ../Results/R.txt
load ../Results/Z.txt

close all
for k=1:18
%    name=['S_',num2str(k,'%.3d'),'.h5'];
%    s=h5read(name,'/Sabs');
%    name=['S_',num2str(k,'%.3d'),'.h5'];
%    F=h5read(name,'/Flux_Ue');
%    plot(s,F,'b-')
%    drawnow
    name=['S_',num2str(k,'%.3d'),'.h5'];    
%     G=h5read(name,'/Gamma');
%    
%     n=h5read(name,'/density');
%   
%     Te=h5read(name,'/Te');
%     
%     Ti=h5read(name,'/Ti');
%     patch(R',Z',(G./(n.*sqrt(Te+Ti)))');shading flat;colorbar
%     caxis([-1 1])

    Sn=h5read(name,'/Sn');
    patch(R',Z',log10(Sn'*1.66e20/1.560314474728064E-004));shading flat;colorbar;
    axis equal;
    caxis([21 24]);
    drawnow
end