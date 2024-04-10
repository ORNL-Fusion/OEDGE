n0=  1.6666667E+20
 T0=   135.373690640864     
 c0=   80537.4225381667     
 tau0=  1.560314474728064E-004

nl=1;

load ps_files/Metrics/G_003.txt
load ps_files/Metrics/G_006.txt
load ps_files/Metrics/G_008.txt


G=[G_003(nl+1,2:end-1),G_006(nl+1,2:end-1),G_008(nl+1,2:end-1)];


[m,p]=size(G);

load ps_files/Mesh/chi_003.txt
load ps_files/Mesh/chi_006.txt
load ps_files/Mesh/chi_008.txt
load ps_files/Mesh/meshzb_003.txt
load ps_files/Mesh/meshzb_006.txt
load ps_files/Mesh/meshzb_008.txt
chi=[chi_003,chi_006,chi_008];
B=zeros(1,p);
for i=1:p
    if(chi(nl,i)==1)
        B(i)=NaN;
    end
end
dz=[meshzb_003(1,2:end)-meshzb_003(1,1:end-1),...
    meshzb_006(1,2:end)-meshzb_006(1,1:end-1),...
    meshzb_008(1,2:end)-meshzb_008(1,1:end-1)]

s=0;
for i=1:p-1
    s=[s,s(i)+2/(G(i)+G(i+1))*(dz(i))];
end

cd(case_dir)
cd Results

load n003.txt
load n006.txt
load n008.txt
n=[n003(nl+1,2:end-1),n006(nl+1,2:end-1),n008(nl+1,2:end-1)];
figure(1)
plot(s,n0*n+B,'b.-')

load M003.txt
load M006.txt
load M008.txt
M=[M003(nl+1,2:end-1),M006(nl+1,2:end-1),M008(nl+1,2:end-1)];
figure(2)
plot(s,M+B,'ro-')
hold on
plot(s,M,'r.-')
ylim([-1.1 1.1])

load Te003.txt
load Te006.txt
load Te008.txt
Te=[Te003(nl+1,2:end-1),Te006(nl+1,2:end-1),Te008(nl+1,2:end-1)];
load Ti003.txt
load Ti006.txt
load Ti008.txt
Ti=[Ti003(nl+1,2:end-1),Ti006(nl+1,2:end-1),Ti008(nl+1,2:end-1)];
figure(3)
plot(s,T0*Te+B,'go-',s,T0*Ti+B,'bo-')
hold on
plot(s,T0*Te,'g.-',s,T0*Ti,'b.-')

load SN003.txt
load SN006.txt
load SN008.txt
SN=[SN003(nl+1,2:end-1),SN006(nl+1,2:end-1),SN008(nl+1,2:end-1)];
figure(4)
plot(s,n0/tau0*SN+B,'m.-')

load SUe003.txt
load SUe006.txt
load SUe008.txt
SUe=[SUe003(nl+1,2:end-1),SUe006(nl+1,2:end-1),SUe008(nl+1,2:end-1)];
figure(5)
plot(s,SUe+B,'c.-')

cd ../..