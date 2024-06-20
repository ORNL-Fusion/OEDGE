 nz=[8];
% nz=[14 6];
% nz=[13 4];

R_=[];
Z_=[];
for k=1:length(nz)
    R=hdf5read('mesh.h5',['/zone',num2str(nz(k)),'/Rcorner']);
    Z=hdf5read('mesh.h5',['/zone',num2str(nz(k)),'/Zcorner']);
    R_=[R_(1:end-1,:);R];
    Z_=[Z_(1:end-1,:);Z];    
end

figure(1)
subplot(1,2,1)
hold on
plot(R_,Z_,'b-')
plot(R_',Z_','b-')
subplot(1,2,2)
hold on
R=R_;
Z=Z_
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

R=R_;
Z=Z_;
[m,p]=size(R);

vx=R(1,2)-R(1,1);
vy=Z(1,2)-Z(1,1);
vx1=R(2,1)-R(1,1);
vy1=Z(2,1)-Z(1,1);
qual=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2));
sens=sign(qual);

R_=R(1,:);
Z_=Z(1,:);
for i=2:m
    din=zeros(size(R(i-1,:)));
    for j=2:p
        din(j)=din(j-1)+sqrt((R_(i-1,j)-R_(i-1,j-1))^2+(Z_(i-1,j)-Z_(i-1,j-1))^2);
    end
    din=din/din(end);
    [dout,Rout,Zout]=optim_mesh(din,R_(i-1,:),Z_(i-1,:),R(i,:),Z(i,:),sens);
    R_=[R_;Rout];
    Z_=[Z_;Zout];
end

figure(2)
subplot(1,2,1)
hold on
plot(R_,Z_,'k-')
plot(R_',Z_','k-')
subplot(1,2,2)
hold on
R=R_;
Z=Z_
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


nz=[9];
% nz=[11];
% nz=[12];

R_=[];
Z_=[];
for k=1:length(nz)
    R=hdf5read('mesh.h5',['/zone',num2str(nz(k)),'/Rcorner']);
    Z=hdf5read('mesh.h5',['/zone',num2str(nz(k)),'/Zcorner']);
    R_=[R_(1:end-1,:);R];
    Z_=[Z_(1:end-1,:);Z];    
end

figure(1)
subplot(1,2,1)
hold on
plot(R_,Z_,'b-')
plot(R_',Z_','b-')
subplot(1,2,2)
hold on
R=R_;
Z=Z_
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
            [qual, qual, qual, qual]*sens)
    end
end
colorbar

R=R_(end:-1:1,:);
Z=Z_(end:-1:1,:);
[m,p]=size(R);

vx=R(1,2)-R(1,1);
vy=Z(1,2)-Z(1,1);
vx1=R(2,1)-R(1,1);
vy1=Z(2,1)-Z(1,1);
qual=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2));
sens=sign(qual);

R_=R(1,:);
Z_=Z(1,:);
for i=2:m
    din=zeros(size(R(i-1,:)));
    for j=2:p
        din(j)=din(j-1)+sqrt((R_(i-1,j)-R_(i-1,j-1))^2+(Z_(i-1,j)-Z_(i-1,j-1))^2);
    end
    din=din/din(end);
    [dout,Rout,Zout]=optim_mesh(din,R_(i-1,:),Z_(i-1,:),R(i,:),Z(i,:),sens);
    R_=[R_;Rout];
    Z_=[Z_;Zout];
end

R_=R_(end:-1:1,:);
Z_=Z_(end:-1:1,:);

figure(2)
subplot(1,2,1)
hold on
plot(R_,Z_,'k-')
plot(R_',Z_','k-')
subplot(1,2,2)
hold on
R=R_;
Z=Z_;
vx=R(1,2)-R(1,1);
vy=Z(1,2)-Z(1,1);
vx1=R(2,1)-R(1,1);
vy1=Z(2,1)-Z(1,1);
qual=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2));
sens=sign(qual);

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
            [qual, qual, qual, qual]*sens)
    end
end
colorbar