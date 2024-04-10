nz=[13 4];

R_=[];
Z_=[];
for k=1:length(nz)
    R=hdf5read('mesh.h5',['/zone',num2str(nz(k)),'/Rcorner']);
    Z=hdf5read('mesh.h5',['/zone',num2str(nz(k)),'/Zcorner']);
    R_=[R_(1:end-1,:);R];
    Z_=[Z_(1:end-1,:);Z];
    
end

figure
hold on
plot(R_,Z_,'b-')
plot(R_',Z_','b-')

R=R_;
Z=Z_;
figure
hold on
[m,p]=size(R);
for i=1:m-1
    for j=1:p-1
        vx=R(i,j+1)-R(i,j);
        vy=Z(i,j+1)-Z(i,j);
        vx1=R(i+1,j)-R(i,j);
        vy1=Z(i+1,j)-Z(i,j);
        qual=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2));
        patch([R(i,j), R(i+1,j), ...
            R(i+1,j+1), R(i,j+1)],...
            [Z(i,j), Z(i+1,j), ...
            Z(i+1,j+1), Z(i,j+1)], ...
            [qual, qual, qual, qual])
    end
end
colorbar
