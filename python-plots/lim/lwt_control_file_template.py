# Name: lwt_control_file.py (lwt = LimWallToolkit)
# Author: Shawn Zamperini
# Email: zamp@utk.edu, zamperinis@fusion.gat.com
# Date: 12/20/21
#
# This is a control file for using the LimWallToolkit. Before using, somewhere
#   on your computer create a working directory and copy/paste this file into
#   it. Put the path to this directory as lwt_control_root under Header Info
#   in this file.
#
#   Example usage:
#   - The script is executed via:
#           python lwt_control_file.py 1
#     where the number following the filename indicates the usage. Options are:
#     1: Generate a 3D wall file for MAFOT by slicing a provided STL mesh at
#          each toroidal angle.
#     2: Generate a .bound for 3DLIM at a specified toroidal angle by using
#          a corresponding MAFOT output file.
#     3: Plot 3DLIM results on the flux surfaces.
#     4: Using a DIVIMP run, generate input for 3DLIM of the injection
#          probability distribution using the actual DIVIMP results.
#
# See each section below for the corresponding input options.
import sys
import numpy as np
try:
    from LimWallToolkit import LimWallToolkit
except:
    print("Error! Cannot find LimWallToolkit. Add utk-fusion/lim to your path.")
    print("  - Example: In .bashrc, add 'export PYTHONPATH=$PYTHONPATH:/Users/zamperini/github/utk-fusion/lim'")
    sys.exit()


#-------------#
# Header Info #
#-------------#
lwt_control_root = "/path/to/working/directory/shot/"

#----------------------------------------------#
# Option 1: Generate a 3D wall file for MAFOT. #
#----------------------------------------------#
stl_path          = "/Users/zamperini/My Drive/Research/Data/vessel_structure/20210826_vessel_mesh.stl"
output_mafot_file = "mafot_3d_wall.dat"
make_gif          = True
gif_path          = "wall_gif/wall.gif"  # Must create a folder called wall_gif.


#------------------------------------------------------------------#
# Option 2: Generate .bound file for 3DLIM input from MAFOT files. #
#------------------------------------------------------------------#

# These files are from running MAFOT in the +1 (top-down CCW) and -1 (top-down
# CW) directions. Which corresponds to file1 and to file2 depends on your setup.
# One way to make it work is file1 can contains all the connection lengths in
# the positive 3DLIM direction (i.e. that which replaces surface 1, L19 in the
# 3DLIM input file), and file2 be the negative values (surface 2, L21).
mafot_file1    = "lam_hires_tor240_conn+1.dat"  # CCW, Reverse = [ITF/OTF]
mafot_file2    = "lam_hires_tor240_conn-1.dat"  # CW,  Reverse = [ITF/OTF]
toroidal_angle = 240                            # MiMES = 240.

# The R bins for 3DLIM. Can have them already set or set here and copy/paste
# into 3DLIM input file.
lim_rbins = np.arange(-0.1, 0.02, 0.0025)

# The number of Pbins MUST match 2*MAXNPS+1 in 3DLIM, otherwise it will not
# work. Probably equal to 41. Also I think it needs to be symmetric. Good to
# manually set to have higher resolution near CP if that is the 3DLIM usage.
# Below example has you set the positive values and just mirrors them.
#lim_pbins = np.linspace(-0.06, 0.06, 41)
tmp = np.array([0.060, 0.050, 0.040, 0.030, 0.020, 0.015, 0.014, 0.013, 0.012, 0.011,
                0.010, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001])
lim_pbins = np.append(-tmp, [0.0])
lim_pbins = np.append(lim_pbins, tmp[::-1])

# In machine coordinates, where do you want the origin of 3DLIM to be?
# MiMES-like origin.
r_origin2 = 2.282    # Often just the R of the CP tip (with any possible shift).
z_origin2 = -0.188
# DiMES-like origin.
#r_origin2 = 1.485
#z_origin2 = -1.15

# Which machine coordinate does the 3DLIM "R" coordinate go along?
# MiMES = "R", DiMES = "Z" and None will set it along the actual R-Rsep
# direction (i.e. the plasma radial direction, but not implemented yet).
along_coord2 = "R"

# Must have .bound extension, i.e. 184527_conns.bound or ramp.bound are okay.
output_bounds_file = "167196_tor240_a2_v3.bound"
gfile_pickle_path2 = "167196_3500.pickle"
wall_path2         = "mafot_3D_wall.dat"


#-----------------------------------------------------------------#
# Option 3: Plot 3DLIM results on a 2D cross section of the wall. #
#-----------------------------------------------------------------#
lim_path     = "/path/to/lim/run.nc"
tor_angle    = 240
r_origin3    = 2.295
z_origin3    = -0.188
wall_path3   = "mafot_3D_wall.dat"
along_coord3 = "R"

# The gfile is a pickled python dictionary obtained from OMFIT. Follow the
# instructions at the bottom of this script to see how to get it.
plot_surfaces = True
gfile_pickle_path3 = "184527_3500.pickle"


#------------------------------------------------------------------------------#
# Option 4: Generate 3DLIM probability distribution for input from a DIVIMP run.
#------------------------------------------------------------------------------#
divimp_nc_path = "/path/to/divimp/run.nc"

# The R of the flux tube you want the probability distribution for.
z_fluxtube     = -0.188
r_fluxtube     = 2.282 - 0.02
show_plot      = True
reverse        = False
smooth         = True
smooth_window  = 15


#-------------#
# Run routine #
#-------------#
if __name__ == "__main__":

    print()
    print("******************")
    print("* LimWallToolkit *")
    print("******************")

    if len(sys.argv) != 2:
        print("Error: See comments at top of file for correct usage.")
        sys.exit()

    lwt = LimWallToolkit()
    if int(sys.argv[1]) == 1:
        print("Option 1: Generate a 3D wall file for MAFOT from an STL file")
        output_mafot_file = lwt_control_root + output_mafot_file
        gif_path = lwt_control_root + gif_path
        lwt.generate_3d_wall_for_mafot(stl_path=stl_path, make_gif=make_gif,
          gif_path=gif_path, output_file=output_mafot_file)

    elif int(sys.argv[1]) == 2:
        print("Option 2: Generate a .bound file for 3DLIM from provided MAFOT files with the connection lengths")
        output_bounds_file = lwt_control_root + output_bounds_file
        gfile_pickle_path2 = lwt_control_root + gfile_pickle_path2
        wall_path2 = lwt_control_root + wall_path2
        dbg = lwt.bounds_file_from_mafot(tor_angle=toroidal_angle,
          mafot_file1=mafot_file1, mafot_file2=mafot_file2, lim_rbins=lim_rbins,
          lim_pbins=lim_pbins, r_origin=r_origin2,
          z_origin=z_origin2, output_file=output_bounds_file,
          gfile_pickle_path=gfile_pickle_path2, along_coord=along_coord2,
          wall_path=wall_path2)

    elif int(sys.argv[1]) == 3:
        print("Option 3: Plot 3DLIM results on a 2D cross section of the wall.")
        wall_path3 = lwt_control_root + wall_path3
        gfile_pickle_path3 = lwt_control_root + gfile_pickle_path3
        lwt.plot_3dlim_on_wall(lim_path=lim_path, tor_angle=tor_angle,
          wall_path=wall_path3,
          plot_surfaces=plot_surfaces, gfile_pickle_path=gfile_pickle_path3,
          r_origin=r_origin3, z_origin=z_origin3, along_coord=along_coord3)

    elif int(sys.argv[1]) == 4:
        print("Option 4: Generate 3DLIM probability distribution for input from a DIVIMP run.")
        lwt.divimp_prob_dist(divimp_nc_path=divimp_nc_path,
          r_fluxtube=r_fluxtube, z_fluxtube=z_fluxtube, show_plot=show_plot,
          reverse=reverse, smooth=smooth, smooth_window=smooth_window)


"""
Additional info
Create gfile pickle object
1. Open OMFIT, load the EFIT module, and open up the GUI.
2. Select shot, time and EFIT grid size (maybe 513x513, why not?)
3. Load SNAP file from MDS+ database.
4. Choose EFIT tree, typically EFIT02 is the goal.
5. Click "Generate g-file...", defaults settings are fine.
6. Copy/paste the following code into the Command Box in OMFIT:

import pickle
import numpy as np
from os.path import expanduser
import omfit_classes
gfile = OMFIT["EFIT"]["FILES"]["gEQDSK"]
shot = int(gfile["CASE"][3].split("#")[1])
time = int(gfile["CASE"][4].split("ms")[0])
home = expanduser("~")
fname = "{}/{}_{}.pickle".format(home, shot, time)
with open(fname, "wb") as f:
    gfile_out = {}
    for k, v in gfile.items():
        if type(v) == omfit_classes.sortedDict.SortedDict:
            for k2, v2 in v.items():
                gfile_out[k2] = v2
        elif type(v) in [int, float, np.ndarray]:
            gfile_out[k] = v
    pickle.dump(gfile_out, f)
print("Saved to {}".format(fname))

7. A pickled gfile will be saved in your home directory. Move this to your
     computer. The path is what gets set to gfile_path.
"""
