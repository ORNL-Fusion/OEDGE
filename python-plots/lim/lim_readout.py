from Readout import Readout
from tkinter import filedialog, Tk
import sys


# Variables for plots.
probe_width    = 0.015
rad_cutoff     = 0.02
sep_vel_plot   = False
sep_force_plot = False
mult_runs      = False
log_center     = True
fit_exp        = True


# Supply the word 'test' to just run test file.
if len(sys.argv) > 1:
    netcdf_path = '/mnt/c/Users/Shawn/Documents/GitHub/Collector-Probes/3DLIM Scripts/colprobe-z1-001e.nc'
    dat_path    = '/mnt/c/Users/Shawn/Documents/GitHub/Collector-Probes/3DLIM Scripts/colprobe-z1-001e.dat'
    lim_path    = '/mnt/c/Users/Shawn/Documents/GitHub/Collector-Probes/3DLIM Scripts/colprobe-z1-001e.lim'

# Ask user for location to files.
else:
    root = Tk(); root.withdraw()
    netcdf_path = filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),))
    #dat_path    = filedialog.askopenfilename(filetypes=(('.dat files',   '*.dat'),))
    #lim_path    = filedialog.askopenfilename(filetypes=(('.lim files',   '*.lim'),))
    dat_path = netcdf_path.split('.nc')[0] + '.dat'
    lim_path = netcdf_path.split('.nc')[0] + '.lim'

    # Option if dat file not included.
    if dat_path == ():
        dat_path = None
    if lim_path == ():
        lim_path = None

# Controlling routine for Readout.
grid = Readout(netcdf_file=netcdf_path, dat_file=dat_path, lim_file=lim_path)
grid.print_readout()

try:
    grid.centerline(0, mult_runs, log=log_center, fit_exp=fit_exp)
except:
    print("Error: Centerline plot.")

#try:
#    grid.avg_imp_vely(1)
#except:
#    print("Error: Average impurity velocity plot.")

try:
    grid.te_contour(1)
except:
    print("Error: Te contour plot.")

try:
    grid.ne_contour(2)
except:
    print("Error: ne contour plot.")

try:
    grid.deposition_contour(3, probe_width=probe_width, rad_cutoff=rad_cutoff, side='ITF')
    grid.deposition_contour(4, probe_width=probe_width, rad_cutoff=rad_cutoff, side='OTF')
except:
    print("Error: Deposition plots.")

try:
    grid.avg_pol_profiles(5, probe_width=probe_width, rad_cutoff=rad_cutoff)
except:
    print("Error: Poloidal plot.")

try:
    #grid.imp_contour_plot(6)
    grid.imp_contour_plot_radial(6, pmin=0.10, pmax=0.13)
    #grid.imp_contour_plot_radial(6, pmin=-0.20, pmax=0.20)
except:
    print('Error: Impurity contour plot.')

try:
    #grid.force_plots(7, separate_plot=sep_force_plot, rad_loc=-0.05)
    grid.force_plots(7, separate_plot=sep_force_plot, rad_loc=0.0)
except:
    print('Error: Force plots.')

try:
    grid.vel_plots(8, vp='vp2', separate_plot=sep_vel_plot)
    #grid.vel_plots(8, vp='vp1', separate_plot=sep_vel_plot)
except:
    print('Error: VP2 plot.')

grid.show_fig()
