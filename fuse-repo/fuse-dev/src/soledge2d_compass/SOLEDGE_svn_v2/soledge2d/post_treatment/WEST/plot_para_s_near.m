n0=  1.6666667E+20;
 T0=   114.015765857941     ;
 c0=   73911.6753502835     ;
 tau0=  1.700187494709652E-004;

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

cd(case_dir)

cd Results

load n002.txt
load n005.txt
load n012.txt
load n018.txt
load n015.txt
n=[n002(nl+1,2:end-1),n005(nl+1,2:end-1),...
    n012(nl+1,2:end-1),n018(nl+1,2:end-1),n015(nl+1,2:end-1)];
figure(1)
plot(s,n0*n+B,'b.-')

load M002.txt
load M005.txt
load M012.txt
load M018.txt
load M015.txt
M=[M002(nl+1,2:end-1),M005(nl+1,2:end-1),...
    M012(nl+1,2:end-1),M018(nl+1,2:end-1),M015(nl+1,2:end-1)];
figure(2)
plot(s,M+B,'ro-')
hold on
plot(s,M,'r.-')
ylim([-1.1 1.1])

load Te002.txt
load Te005.txt
load Te012.txt
load Te018.txt
load Te015.txt
Te=[Te002(nl+1,2:end-1),Te005(nl+1,2:end-1),...
    Te012(nl+1,2:end-1),Te018(nl+1,2:end-1),Te015(nl+1,2:end-1)];
load Ti002.txt
load Ti005.txt
load Ti012.txt
load Ti018.txt
load Ti015.txt
Ti=[Ti002(nl+1,2:end-1),Ti005(nl+1,2:end-1),...
    Ti012(nl+1,2:end-1),Ti018(nl+1,2:end-1),Ti015(nl+1,2:end-1)];
figure(3)
plot(s,T0*Te+B,'go-',s,T0*Ti+B,'bo-')
hold on
plot(s,T0*Te,'g.-',s,T0*Ti,'b.-')

load SN002.txt
load SN005.txt
load SN012.txt
load SN018.txt
load SN015.txt
SN=[SN002(nl+1,2:end-1),SN005(nl+1,2:end-1),...
    SN012(nl+1,2:end-1),SN018(nl+1,2:end-1),SN015(nl+1,2:end-1)];
figure(4)
plot(s,n0/tau0*SN+B,'m.-')

load SUe002.txt
load SUe005.txt
load SUe012.txt
load SUe018.txt
load SUe015.txt
SUe=[SUe002(nl+1,2:end-1),SUe005(nl+1,2:end-1),...
    SUe012(nl+1,2:end-1),SUe018(nl+1,2:end-1),SUe015(nl+1,2:end-1)];
figure(5)
plot(s,SUe+B,'c.-')

cd ../..