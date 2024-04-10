hdf5write('mesh.h5','/Nzones',zones.num);
for k=1:zones.num
    if(~isempty(zones.zone(k).gridR))
        hdf5write('mesh.h5',['/zone',num2str(k),'/Rcorner'],zones.zone(k).gridR,'WriteMode','append');
        hdf5write('mesh.h5',['/zone',num2str(k),'/Zcorner'],zones.zone(k).gridZ,'WriteMode','append');
    end
end