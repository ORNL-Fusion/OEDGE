import h5py
import numpy as np

err_cpp=np.zeros(8, dtype=float)
err_cppr=np.zeros(8, dtype=float)
err_ctt=np.zeros(8, dtype=float)
err_cttr=np.zeros(8, dtype=float)
err_cpt=np.zeros(8, dtype=float)
err_cptr=np.zeros(8, dtype=float)
err_Jac=np.zeros(8, dtype=float)
err_Jacr=np.zeros(8, dtype=float)
err_G=np.zeros(8, dtype=float)
err_Gr=np.zeros(8, dtype=float)
res=np.zeros(8, dtype=float)

dt=1E-5

for k in range(2,10):
    res[k-2]=2**(k+1)
    Nx=2**(k+1)
    f=h5py.File('Results/'+str(k)+'_metric','r')
# cpp
    dataset = f['/zone1/cpp']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/cpp2']
    n2 = np.zeros((Nx/2,Nx,), dtype=float)
    n2=dataset[...]
    print(n)
    error=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    print(error)
    error_rel=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    err_cpp[k-2]=error
    err_cppr[k-2]=error_rel
# cpt
    dataset = f['/zone1/cpt']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/cpt2']
    n2 = np.zeros((Nx/2,Nx,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    err_cpt[k-2]=error
    err_cptr[k-2]=error_rel
# ctt
    dataset = f['/zone1/ctt']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/ctt2']
    n2 = np.zeros((Nx/2,Nx,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    err_ctt[k-2]=error
    err_cttr[k-2]=error_rel
# Jac
    dataset = f['/zone1/Jac']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/Jac2']
    n2 = np.zeros((Nx/2,Nx,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    err_Jac[k-2]=error
    err_Jacr[k-2]=error_rel
# G
    dataset = f['/zone1/G']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/G2']
    n2 = np.zeros((Nx/2,Nx,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[0:Nx-1,0:Nx/2-1]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    err_G[k-2]=error
    err_Gr[k-2]=error_rel
    f.close()


import matplotlib as mpl
mpl.use('Agg')
from pylab import *
import matplotlib.pyplot  as plt
figure(0)
line_n,= plt.loglog(res,err_cpp, 'ks-', lw=2)
line_n2,= plt.loglog(res,err_cpt, 'rv-', lw=2)
line_n3,= plt.loglog(res,err_ctt, 'b^-', lw=2)
line_n4,= plt.loglog(res,err_Jac, 'go-', lw=2)
line_n5,= plt.loglog(res,err_G, 'mp-', lw=2)
plt.grid(True,which="both",ls="-")
plt.legend([line_n, line_n2, line_n3, line_n4, line_n5],['cpp','cpt','ctt','Jac','G'])
plt.title('Metric - Absolute error / Grid3')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')
#plt.show()
plt.savefig('errors/err_abs_metric_grid3.png')

figure(1)
line_n,= plt.loglog(res,err_cppr, 'ks-', lw=2)
line_n2,= plt.loglog(res,err_cptr, 'rv-', lw=2)
line_n3,= plt.loglog(res,err_cttr, 'b^-', lw=2)
line_n4,= plt.loglog(res,err_Jacr, 'go-', lw=2)
line_n5,= plt.loglog(res,err_Gr, 'mp-', lw=2)
plt.grid(True,which="both",ls="-")
plt.legend([line_n, line_n2, line_n3, line_n4, line_n5],['cpp','cpt','ctt','Jac','G'])
plt.title('Metric - Relative error / Grid3')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')
#plt.show()
plt.savefig('errors/err_rel_metric_grid3.png')
