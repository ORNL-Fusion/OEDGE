import h5py
import numpy as np

errn=np.zeros(8, dtype=float)
errnr=np.zeros(8, dtype=float)
errG=np.zeros(8, dtype=float)
errGr=np.zeros(8, dtype=float)
errT=np.zeros(8, dtype=float)
errTe=np.zeros(8, dtype=float)
errTr=np.zeros(8, dtype=float)
errTer=np.zeros(8, dtype=float)
res=np.zeros(8, dtype=float)
errn2=np.zeros(8, dtype=float)
errn2r=np.zeros(8, dtype=float)
errG2=np.zeros(8, dtype=float)
errG2r=np.zeros(8, dtype=float)
errT2=np.zeros(8, dtype=float)
errT2r=np.zeros(8, dtype=float)



dt=1E-5

for k in range(2,10):
    res[k-2]=2**(k+1)
    Nx=2**(k+1)
    f=h5py.File('Results/'+str(k)+'_plasma_1','r')
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
#density
    dataset = f['/zone1/density']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/density2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errn[k-2]=error/dt
    errnr[k-2]=error_rel/dt
#velocity
    dataset = f['/zone1/Gamma']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/Gamma2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,2:Nx/2]-n[1:Nx,2:Nx/2]))/((Nx/2-1)*Nx)
    error_rel=np.sum(abs(n2[1:Nx,2:Nx/2]-n[1:Nx,2:Nx/2])/abs(n[1:Nx,2:Nx/2]))/((Nx/2-1)*Nx)
    errG[k-2]=error/dt
    errGr[k-2]=error_rel/dt
    f.close()

    f=h5py.File('Results/'+str(k)+'_plasma_2','r')
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
#density2
    dataset = f['/zone1/density']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/density2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2]))/(Nx**2/2)
    error_rel=np.sum(abs(n2[1:Nx,1:Nx/2]-n[1:Nx,1:Nx/2])/abs(n[1:Nx,1:Nx/2]))/(Nx**2/2)
    errn2[k-2]=error/dt
    errn2r[k-2]=error_rel/dt
#velocity2
    dataset = f['/zone1/Gamma']
    n = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n=dataset[...]
    dataset = f['/zone1/Gamma2']
    n2 = np.zeros((Nx/2+2,Nx+2,), dtype=float)
    n2=dataset[...]
    error=np.sum(abs(n2[1:Nx,2:Nx/2]-n[1:Nx,2:Nx/2]))/((Nx/2-1)*Nx)
    error_rel=np.sum(abs(n2[1:Nx,2:Nx/2]-n[1:Nx,2:Nx/2])/abs(n[1:Nx,2:Nx/2]))/((Nx/2-1)*Nx)
    errG2[k-2]=error/dt
    errG2r[k-2]=error_rel/dt
    f.close()


#electron temperature
    f=h5py.File('Results/'+str(k)+'_plasma_0','r')
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
line_n,= plt.loglog(res,errn, 'ks-', lw=2)
line_G,= plt.loglog(res,errG, 'gv-', lw=2)
line_T,= plt.loglog(res,errT, 'b^-', lw=2)
line_Te,= plt.loglog(res,errTe, 'ro-', lw=2)
line_n2,= plt.loglog(res,errn2, 'ks--', lw=3)
line_G2,= plt.loglog(res,errG2, 'gv--', lw=3)
line_T2,= plt.loglog(res,errT2, 'b^--', lw=3)
plt.grid(True,which="both",ls="-")
plt.legend([line_n, line_G, line_T, line_Te, line_n2, line_G2, line_T2],['Density','Velocity','Ti','Te','Density2','Velocity2','Ti2'])
plt.title('Collisions - Absolute error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')
plt.savefig('errors/err_abs_col_grid2.png')

figure(1)
line_n,= plt.loglog(res,errnr, 'ks-', lw=2)
line_G,= plt.loglog(res,errGr, 'gv-', lw=2)
line_T,= plt.loglog(res,errTr, 'b^-', lw=2)
line_Te,= plt.loglog(res,errTer, 'ro-', lw=2)
line_n2,= plt.loglog(res,errn2r, 'ks--', lw=3)
line_G2,= plt.loglog(res,errG2r, 'gv--', lw=3)
line_T2,= plt.loglog(res,errT2r, 'b^--', lw=3)
plt.grid(True,which="both",ls="-")
plt.legend([line_n, line_G, line_T, line_Te, line_n2, line_G2, line_T2],['Density','Velocity','Ti','Te','Density2','Velocity2','Ti2'])
plt.title('Collisions - Relative error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')
plt.savefig('errors/err_rel_col_grid2.png')
