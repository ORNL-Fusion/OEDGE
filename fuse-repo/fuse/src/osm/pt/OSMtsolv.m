%=========================================================
%  OSMtsolv
%  Etude des sorties du code transp_parallel.f90
%
%=========================================================


clear all

% Little or big endian ??
endian = 'l';
%endian = 'b';

Nz = 200
dtdiag = 10

% File location
folder = ['Result/'];

% File names
fileN = [folder,'N.dat'];
fileV = [folder,'V.dat'];
fileTe = [folder,'Te.dat'];
fileTi = [folder,'Ti.dat'];
filez = [folder,'zmesh.dat'];
fileB = [folder,'Bamp.dat'];

% Determining number of time steps saved
fidsize = fopen(fileN,'r',endian);
usize=fread(fidsize,inf,'uint32');
fclose(fidsize);
Nt=size(usize,1)/(2*Nz+4)
tt = [0:Nt-1]*dtdiag;
clear usize

% Calling my general loading routine
N = mybinload(fileN,Nt,[Nz],'float64',endian);


V = mybinload(fileV,Nt,[Nz],'float64',endian);
Te = mybinload(fileTe,Nt,[Nz],'float64',endian);
Ti = mybinload(fileTi,Nt,[Nz],'float64',endian);
zz = mybinload(filez,1,[Nz],'float64',endian);
Bamp = mybinload(fileB,1,[Nz],'float64',endian);

period = 2*zz(Nz)-zz(Nz-1)

V = V./(sqrt((Te+Ti)/(3.345E-27/1.602E-19)));

Gamma = V;
%  for it = 1:Nt
%     Gamma(it,:)=Gamma(it,:)./Bamp;
%  end

% Analytical solution
%  Vanal = 2*(zz-0.5)./(1+sqrt(1-(2*(zz-0.5)).^2));
%  Nanal = (zz-0.5)./Vanal;
%  Nanal = Nanal*(1.E25/(2*Vanal(Nz)))/Nanal(1);


% 1D plots
figure(1)
plot(zz,N(1:floor(Nt/10):Nt,:))%,'+')
%plot(zz,Nanal,'o')
%  plot(zz,N(1:floor(Nt/10):Nt,Nz:-1:1),'o')
%  hold off

figure(2)
plot(zz,Gamma(1:floor(Nt/10):Nt,:))%,'+')
%plot(zz,Vanal,'o')
%  plot(zz,-V(1:floor(Nt/10):Nt,Nz:-1:1),'o')
%  hold off

figure(3)
plot(zz,Te(1:floor(Nt/10):Nt,:))%,'+')

figure(4)
plot(zz,Ti(1:floor(Nt/10):Nt,:))%,'+')


% 2D plotsfigure(2)
plot(zz,Gamma(1:floor(Nt/10):Nt,:))%,'+')
figure(5)
pcolor(zz,tt,N)
shading flat
colorbar

figure(6)
pcolor(zz,tt,Gamma)
shading flat
colorbar

figure(7)
pcolor(zz,tt,Te)
shading flat
colorbar

figure(8)
pcolor(zz,tt,Ti)
shading flat
colorbar

% time plots
figure(9)
plot(tt,N(:,1:floor(Nz/10):Nz)')%,'+')

figure(10)
plot(tt,V(:,1:floor(Nz/10):Nz)')%,'+')

figure(11)
plot(tt,Te(:,1:floor(Nz/10):Nz)')%,'+')

figure(12)
plot(tt,Ti(:,1:floor(Nz/10):Nz)')%,'+')


% conservation checks

% we need the volume of each cell
Bamp_edge = zeros(1,Nz+1);
zz_edge = zeros(1,Nz+1);
cellvol = zeros(Nt,Nz);
for iz = 2:Nz
   Bamp_edge(iz) = 0.5*(Bamp(iz-1)+Bamp(iz));
   zz_edge(iz) = 0.5*(zz(iz-1)+zz(iz));
end
if (period==0)
   Bamp_edge(1) = 0.5*(3*Bamp(1)-Bamp(2));
   Bamp_edge(Nz+1) = 0.5*(3*Bamp(Nz)-Bamp(Nz-1));
   zz_edge(1) = 0.5*(3*zz(1)-zz(2));
   zz_edge(Nz+1) = 0.5*(3*zz(Nz)-zz(Nz-1));
else
   Bamp_edge(1) = 0.5*(Bamp(1)+Bamp(Nz));
   Bamp_edge(Nz+1) = Bamp_edge(1);
   zz_edge(1) = 0.5*(zz(1)+zz(Nz)-period);
   zz_edge(Nz+1) = zz_edge(1)+period;
end
for iz = 1:Nz
   if Bamp_edge(iz+1)~=Bamp_edge(iz)
      cellvol(1:Nt,iz) = (zz_edge(iz+1)-zz_edge(iz))/Bamp(iz)*log(Bamp_edge(iz+1)/Bamp_edge(iz))/(Bamp_edge(iz+1)/Bamp_edge(iz)-1);
   else
      cellvol(1:Nt,iz) = (zz_edge(iz+1)-zz_edge(iz))/Bamp(iz);
   end
end

matter = sum(N.*cellvol,2);
figure(13)
plot(tt,matter,'-+')
print, 'Matter balance',(max(matter)-min(matter))/matter(1)

momentum = sum(N.*V.*cellvol,2);
figure(14)
plot(tt,momentum,'-+')
print, 'Momentum balance',(max(momentum)-min(momentum))/max(abs(momentum))

e_energy = sum(N.*Te.*cellvol,2);
figure(15)
plot(tt,e_energy,'-+')
print, 'Electrons energy balance', (max(e_energy)-min(e_energy))/max(abs(e_energy))

i_energy = sum(N.*Ti.*cellvol,2);
figure(16)
plot(tt,i_energy,'-+')
print, 'Ions energy balance', (max(i_energy)-min(i_energy))/max(abs(i_energy))


%  % a little movie
%  figure(9)
%  trac1D = N;
%  mintrac = min(min(trac1D));
%  maxtrac = max(max(trac1D));
%  for it = 1:Nt
%     plot(zz,trac1D(it,:))
%     axis([min(zz) max(zz) mintrac maxtrac])
%     film1D(it) = getframe;
%  end
