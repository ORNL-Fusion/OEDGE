import h5py
import numpy as np

errT1=np.zeros(8, dtype=float)
errTe1=np.zeros(8, dtype=float)
errTr1=np.zeros(8, dtype=float)
errTer1=np.zeros(8, dtype=float)
errT21=np.zeros(8, dtype=float)
errT2r1=np.zeros(8, dtype=float)

errT2=np.zeros(8, dtype=float)
errTe2=np.zeros(8, dtype=float)
errTr2=np.zeros(8, dtype=float)
errTer2=np.zeros(8, dtype=float)
errT22=np.zeros(8, dtype=float)
errT2r2=np.zeros(8, dtype=float)

errT3=np.zeros(8, dtype=float)
errTe3=np.zeros(8, dtype=float)
errTr3=np.zeros(8, dtype=float)
errTer3=np.zeros(8, dtype=float)
errT23=np.zeros(8, dtype=float)
errT2r3=np.zeros(8, dtype=float)
res=np.zeros(8, dtype=float)



dt=1E-5

for k in range(2,9):
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
    errT1[k-2]=error/dt
    errTr1[k-2]=error_rel/dt
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
    errT21[k-2]=error/dt
    errT2r1[k-2]=error_rel/dt
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
    errTe1[k-2]=error/dt
    errTer1[k-2]=error_rel/dt
    f.close()


import matplotlib as mpl
mpl.use('Agg')
from pylab import *
import matplotlib.pyplot  as plt
figure(0)
line_Te1,= plt.loglog(res,errTe1, 'ro-', lw=2, label='$T_e$ part1')
line_T1,= plt.loglog(res,errT1, 'rv-.', lw=2, label='$T_{He^+}$ part1')
line_T21,= plt.loglog(res,errT21, 'rs--', lw=2, label='$T_{He^{2+}}$ part1')
plt.grid(True,which="both",ls="-")
plt.title('Collisions - Absolute error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')

figure(1)
line_Ter1,= plt.loglog(res,errTer1, 'ro-', lw=2, label='$T_e$ part1')
line_Tr1,= plt.loglog(res,errTr1, 'rv-.', lw=2, label='$T_{He^+}$ part1')
line_T2r1,= plt.loglog(res,errT2r1, 'rs--', lw=2, label='$T_{He^{2+}}$ part1')
plt.grid(True,which="both",ls="-")
plt.title('Collisions - Relative error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')


for k in range(2,9):
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
    errT2[k-2]=error/dt
    errTr2[k-2]=error_rel/dt
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
    errT22[k-2]=error/dt
    errT2r2[k-2]=error_rel/dt
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
    errTe2[k-2]=error/dt
    errTer2[k-2]=error_rel/dt
    f.close()

figure(0)
line_Te2,= plt.loglog(res,errTe2, 'ko-', lw=2, label='$T_e$ part2')
line_T2,= plt.loglog(res,errT2, 'kv-.', lw=2, label='$T_{He^+}$ part2')
line_T22,= plt.loglog(res,errT22, 'ks--', lw=2, label='$T_{He^{2+}}$ part2')
plt.grid(True,which="both",ls="-")
plt.title('Collisions - Absolute error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')

figure(1)
line_Ter2,= plt.loglog(res,errTer2, 'ko-', lw=2, label='$T_e$ part2')
line_Tr2,= plt.loglog(res,errTr2, 'kv-.', lw=2, label='$T_{He^+}$ part2')
line_T2r2,= plt.loglog(res,errT2r2, 'ks--', lw=2, label='$T_{He^{2+}}$ part2')
plt.grid(True,which="both",ls="-")
plt.title('Collisions - Relative error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')

for k in range(2,9):
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
    errT3[k-2]=error/dt
    errTr3[k-2]=error_rel/dt
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
    errT23[k-2]=error/dt
    errT2r3[k-2]=error_rel/dt
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
    errTe3[k-2]=error/dt
    errTer3[k-2]=error_rel/dt
    f.close()


figure(0)
line_Te3,= plt.loglog(res,errTe3, 'bo-', lw=2, label='$T_e$ part3')
line_T3,= plt.loglog(res,errT3, 'bv-.', lw=2, label='$T_{He^+}$ part3')
line_T23,= plt.loglog(res,errT23, 'bs--', lw=2, label='$T_{He^{2+}}$ part3')
plt.grid(True,which="both",ls="-")
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.gca().legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Collisions - Absolute error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')

figure(1)
line_Ter3,= plt.loglog(res,errTer3, 'bo-', lw=2, label='$T_e$ part3')
line_Tr3,= plt.loglog(res,errTr3, 'bv-.', lw=2, label='$T_{He^+}$ part3')
line_T2r3,= plt.loglog(res,errT2r3, 'bs--', lw=2, label='$T_{He^{2+}}$ part3')
plt.grid(True,which="both",ls="-")
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.gca().legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Collisions - Relative error / Grid2')
plt.xlabel('Grid resolution')
plt.ylabel('Truncation error')


figure(0)
plt.savefig('errors/err_abs_col_grid2T.png')

figure(1)
plt.savefig('errors/err_rel_col_grid2T.png')
