import h5py
import numpy as np

errn=np.zeros(8, dtype=float)
errnr=np.zeros(8, dtype=float)
errG=np.zeros(8, dtype=float)
errGr=np.zeros(8, dtype=float)
res=np.zeros(8, dtype=float)
errn2=np.zeros(8, dtype=float)
errn2r=np.zeros(8, dtype=float)
errG2=np.zeros(8, dtype=float)
errG2r=np.zeros(8, dtype=float)


dt=1E-5

for k in range(2,9):
    res[k-2]=2**(k+1)
    Nx=2**(k+1)
    f=h5py.File('Results1/'+str(k)+'_plasma_1','r')
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

    f=h5py.File('Results1/'+str(k)+'_plasma_2','r')
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


import matplotlib as mpl
mpl.use('Agg')
from pylab import *
import matplotlib.pyplot  as plt
figure(0)
line_n,= plt.loglog(res,errn, 'ks-', lw=2, label='Density $He^+$')
line_G,= plt.loglog(res,errG, 'gv-', lw=2, label='Velocity $He^+$')
line_n2,= plt.loglog(res,errn2, 'ks--', lw=3, label='Density $He^{2+}$')
line_G2,= plt.loglog(res,errG2, 'gv--', lw=3, label='Velocity $He^{2+}$')
plt.grid(True,which="both",ls="-")
#plt.legend([line_n, line_G, line_n2, line_G2],['Density $He^+$','Velocity $He^+$', 'Density $He^{2+}$','Velocity $He^{2+}$'])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.gca().legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Collisions - Absolute error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')
plt.savefig('errors/err_abs_col_grid2v.png')

figure(1)
line_n,= plt.loglog(res,errnr, 'ks-', lw=2, label='Density $He^+$')
line_G,= plt.loglog(res,errGr, 'gv-', lw=2, label='Velocity $He^+$')
line_n2,= plt.loglog(res,errn2r, 'ks--', lw=3, label='Density $He^{2+}$')
line_G2,= plt.loglog(res,errG2r, 'gv--', lw=3, label='Velocity $He^{2+}$')
plt.grid(True,which="both",ls="-")
#plt.legend([line_n, line_G, line_n2, line_G2],['Density $He^+$','Velocity $He^+$', 'Density $He^{2+}$','Velocity $He^{2+}$'])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.gca().legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Collisions - Relative error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')
plt.savefig('errors/err_rel_col_grid2v.png')
