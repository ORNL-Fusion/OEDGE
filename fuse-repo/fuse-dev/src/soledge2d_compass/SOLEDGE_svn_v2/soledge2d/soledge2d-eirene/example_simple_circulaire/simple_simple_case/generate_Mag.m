rmin=0.2;
rmax=0.6;
R0=2;
Btor0=4;

cd Mesh

load meshx_001.txt
load meshz_001.txt
theta=2*pi*meshz_001+pi/2;
r1=rmin+meshx_001*(rmax-rmin);

R1=R0+r1.*cos(theta);
Z1=r1.*sin(theta);

load meshx_002.txt
load meshz_002.txt
theta=2*pi*meshz_002+pi/2;
r2=rmin+meshx_002*(rmax-rmin);

R2=R0+r2.*cos(theta);
Z2=r2.*sin(theta);

[m,p]=size(meshx_002);

rs1=rmin*ones(1,p);
ths1=theta(1,:);
Rs1=R0+rs1.*cos(ths1);
Zs1=rs1.*sin(ths1);

rn1=(rmin+rmax)/2*ones(1,p);
thn1=theta(1,:);
Rn1=R0+rn1.*cos(thn1);
Zn1=rn1.*sin(thn1);

rs2=rn1;
ths2=thn1;
Rs2=Rn1;
Zs2=Zn1;

rn2=rmax*ones(1,p);
thn2=theta(1,:);
Rn2=R0+rn2.*cos(thn2);
Zn2=rn2.*sin(thn2);

Bphi1=Btor0*R0./R1;
Br1=-Btor0/4*rmin^2./(R1.*r1).*sin(theta);
Bz1=Btor0/4*rmin^2./(R1.*r1).*cos(theta);

Bphi2=Btor0*R0./R2;
Br2=-Btor0/4*rmin^2./(R2.*r2).*sin(theta);
Bz2=Btor0/4*rmin^2./(R2.*r2).*cos(theta);

Bphis1=Btor0*R0./Rs1;
Brs1=-Btor0/4*rmin^2./(Rs1.*rs1).*sin(ths1);
Bzs1=Btor0/4*rmin^2./(Rs1.*rs1).*cos(ths1);

Bphin1=Btor0*R0./Rn1;
Brn1=-Btor0/4*rmin^2./(Rn1.*rn1).*sin(thn1);
Bzn1=Btor0/4*rmin^2./(Rn1.*rn1).*cos(thn1);

Bphis2=Btor0*R0./Rs2;
Brs2=-Btor0/4*rmin^2./(Rs2.*rs2).*sin(ths2);
Bzs2=Btor0/4*rmin^2./(Rs2.*rs2).*cos(ths2);

Bphin2=Btor0*R0./Rn2;
Brn2=-Btor0/4*rmin^2./(Rn2.*rn2).*sin(thn2);
Bzn2=Btor0/4*rmin^2./(Rn2.*rn2).*cos(thn2);


cd ..
!mkdir Magnetic_input
cd Magnetic_input
save 'R_001.txt' 'R1' -ascii -double;
save 'Z_001.txt' 'Z1' -ascii -double;
save 'R_002.txt' 'R2' -ascii -double;
save 'Z_002.txt' 'Z2' -ascii -double;
save 'Bp_001.txt' 'Bphi1' -ascii -double;
save 'Bp_002.txt' 'Bphi2' -ascii -double;
save 'Br_001.txt' 'Br1' -ascii -double;
save 'Br_002.txt' 'Br2' -ascii -double;
save 'Bz_001.txt' 'Bz1' -ascii -double;
save 'Bz_002.txt' 'Bz2' -ascii -double;

!mkdir Neighs
cd Neighs
fid=fopen('North_001.txt','w');
for i=1:p
   fprintf(fid,'%18.15e\t',Brn1(i)); 
   fprintf(fid,'%18.15e\t',Bzn1(i)); 
   fprintf(fid,'%18.15e\t',Rn1(i)); 
   fprintf(fid,'%18.15e\t',Zn1(i)); 
   fprintf(fid,'%18.15e\t',Bphin1(i)); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'\n');
end
fclose(fid)
fid=fopen('South_001.txt','w');
for i=1:p
   fprintf(fid,'%18.15e\t',Brs1(i)); 
   fprintf(fid,'%18.15e\t',Bzs1(i)); 
   fprintf(fid,'%18.15e\t',Rs1(i)); 
   fprintf(fid,'%18.15e\t',Zs1(i)); 
   fprintf(fid,'%18.15e\t',Bphis1(i)); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'\n');
end
fclose(fid)
fid=fopen('North_002.txt','w');
for i=1:p
   fprintf(fid,'%18.15e\t',Brn2(i)); 
   fprintf(fid,'%18.15e\t',Bzn2(i)); 
   fprintf(fid,'%18.15e\t',Rn2(i)); 
   fprintf(fid,'%18.15e\t',Zn2(i)); 
   fprintf(fid,'%18.15e\t',Bphin2(i)); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'\n');
end
fclose(fid)
fid=fopen('South_002.txt','w');
for i=1:p
   fprintf(fid,'%18.15e\t',Brs2(i)); 
   fprintf(fid,'%18.15e\t',Bzs2(i)); 
   fprintf(fid,'%18.15e\t',Rs2(i)); 
   fprintf(fid,'%18.15e\t',Zs2(i)); 
   fprintf(fid,'%18.15e\t',Bphis2(i)); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'\n');
end
fclose(fid)

fid=fopen('East_001.txt','w');
for i=1:16
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'\n');
end
fclose(fid)
fid=fopen('East_002.txt','w');
for i=1:16
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'\n');
end
fclose(fid)
fid=fopen('West_001.txt','w');
for i=1:16
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'\n');
end
fclose(fid)
fid=fopen('West_002.txt','w');
for i=1:16
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'%18.15e\t',0); 
   fprintf(fid,'\n');
end
fclose(fid)

cd ..

cd ..
