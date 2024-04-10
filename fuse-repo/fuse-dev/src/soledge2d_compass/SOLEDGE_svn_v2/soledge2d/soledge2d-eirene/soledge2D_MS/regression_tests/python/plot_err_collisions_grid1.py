import h5py
import numpy as np

errT=np.zeros(8, dtype=float)
errTe=np.zeros(8, dtype=float)
errTr=np.zeros(8, dtype=float)
errTer=np.zeros(8, dtype=float)
errT2=np.zeros(8, dtype=float)
errT2r=np.zeros(8, dtype=float)
res=np.zeros(8, dtype=float)



dt=1E-5

for k in range(2,10):
    res[k-2]=2**(k+1)
    Nx=2**(k+1)
    f=h5py.File('Results1/'+str(k)+'_plasma_1','r')
# temperature ions
    dataset = f['/zone1/temperature']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/temperature2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errT[k-2]=error/dt
    errTr[k-2]=error_rel/dt
    f.close()

    f=h5py.File('Results1/'+str(k)+'_plasma_2','r')
# temperature ions2
    dataset = f['/zone1/temperature']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/temperature2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errT2[k-2]=error/dt
    errT2r[k-2]=error_rel/dt
    f.close()

#electron temperature
    f=h5py.File('Results1/'+str(k)+'_plasma_0','r')
    dataset = f['/zone1/temperature']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/temperature2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errTe[k-2]=error/dt
    errTer[k-2]=error_rel/dt
    f.close()


import matplotlib as mpl
mpl.use('Agg')
from pylab import *
import matplotlib.pyplot  as plt
figure(0)
line_Te,= plt.loglog(res,errTe, 'ro-', lw=2)
line_T,= plt.loglog(res,errT, 'rv-.', lw=2)
line_T2,= plt.loglog(res,errT2, 'rs--', lw=2)
plt.grid(True,which="both",ls="-")
plt.legend([line_Te, line_T, line_T2],['Te p1','Ti $He^+$ p1','Ti $He^{2+}$ p1'])
plt.title('Collisions - Absolute error / Grid1')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')

figure(1)
line_Te,= plt.loglog(res,errTer, 'ro-', lw=2)
line_T,= plt.loglog(res,errTr, 'rv-.', lw=2)
line_T2,= plt.loglog(res,errT2r, 'rs--', lw=2)
plt.grid(True,which="both",ls="-")
plt.legend([line_Te, line_T, line_T2],['Te p1','Ti $He^+$ p1','Ti $He^{2+}$ p1'])
plt.title('Collisions - Relative error / Grid1')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')


for k in range(2,10):
    res[k-2]=2**(k+1)
    Nx=2**(k+1)
    f=h5py.File('Results2/'+str(k)+'_plasma_1','r')
# temperature ions
    dataset = f['/zone1/temperature']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/temperature2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errT[k-2]=error/dt
    errTr[k-2]=error_rel/dt
    f.close()

    f=h5py.File('Results2/'+str(k)+'_plasma_2','r')
# temperature ions2
    dataset = f['/zone1/temperature']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/temperature2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errT2[k-2]=error/dt
    errT2r[k-2]=error_rel/dt
    f.close()

#electron temperature
    f=h5py.File('Results2/'+str(k)+'_plasma_0','r')
    dataset = f['/zone1/temperature']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/temperature2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errTe[k-2]=error/dt
    errTer[k-2]=error_rel/dt
    f.close()

figure(0)
line_Te,= plt.loglog(res,errTe, 'go-', lw=2)
line_T,= plt.loglog(res,errT, 'gv-.', lw=2)
line_T2,= plt.loglog(res,errT2, 'gs--', lw=2)
plt.grid(True,which="both",ls="-")
plt.legend([line_Te, line_T, line_T2],['Te p2','Ti $He^+$ p2','Ti $He^{2+}$ p2'])
plt.title('Collisions - Absolute error / Grid1')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')

figure(1)
line_Te,= plt.loglog(res,errTer, 'go-', lw=2)
line_T,= plt.loglog(res,errTr, 'gv-.', lw=2)
line_T2,= plt.loglog(res,errT2r, 'gs--', lw=2)
plt.grid(True,which="both",ls="-")
plt.legend([line_Te, line_T, line_T2],['Te p2','Ti $He^+$ p2','Ti $He^{2+}$ p2'])
plt.title('Collisions - Relative error / Grid1')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')

for k in range(2,10):
    res[k-2]=2**(k+1)
    Nx=2**(k+1)
    f=h5py.File('Results3/'+str(k)+'_plasma_1','r')
# temperature ions
    dataset = f['/zone1/temperature']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/temperature2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errT[k-2]=error/dt
    errTr[k-2]=error_rel/dt
    f.close()

    f=h5py.File('Results3/'+str(k)+'_plasma_2','r')
# temperature ions2
    dataset = f['/zone1/temperature']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/temperature2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errT2[k-2]=error/dt
    errT2r[k-2]=error_rel/dt
    f.close()

#electron temperature
    f=h5py.File('Results3/'+str(k)+'_plasma_0','r')
    dataset = f['/zone1/temperature']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/temperature2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errTe[k-2]=error/dt
    errTer[k-2]=error_rel/dt
    f.close()


figure(0)
line_Te,= plt.loglog(res,errTe, 'bo-', lw=2)
line_T,= plt.loglog(res,errT, 'bv-.', lw=2)
line_T2,= plt.loglog(res,errT2, 'bs--', lw=2)
plt.grid(True,which="both",ls="-")
plt.legend([line_Te, line_T, line_T2],['Te p3','Ti $He^+$ p3','Ti $He^{2+}$ p3'])
plt.title('Collisions - Absolute error / Grid1')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')

figure(1)
line_Te,= plt.loglog(res,errTer, 'bo-', lw=2)
line_T,= plt.loglog(res,errTr, 'bv-.', lw=2)
line_T2,= plt.loglog(res,errT2r, 'bs--', lw=2)
plt.grid(True,which="both",ls="-")
plt.legend([line_Te, line_T, line_T2],['Te p3','Ti $He^+$ p3','Ti $He^{2+}$ p3'])
plt.title('Collisions - Relative error / Grid1')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')


figure(0)
plt.savefig('errors/err_abs_col_grid1T.png')

figure(1)
plt.savefig('errors/err_rel_col_grid1T.png')
