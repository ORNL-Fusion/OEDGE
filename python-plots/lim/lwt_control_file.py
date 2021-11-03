# Name: lwt_control_file.py (lwt = LimWallToolkit)
# Author: Shawn Zamperini
# Email: zamp@utk.edu, zamperinis@fusion.gat.com
# Date: 10/11/21
#
# This is a control file for using the LimWallToolkit. Example usage:
#   - The script is executed via:
#           python lwt_control_file.py 1
#     where the number following the filename indicates the usage. Options are:
#     1: Generate a 3D wall file for MAFOT by slicing a provided STL mesh at
#          each toroidal angle.
#     2: Generate a .bound for 3DLIM at a specified toroidal angle by using
#          a corresponding MAFOT output file.
#
# See each section below for the corresponding input options.
import sys
import numpy as np
from LimWallToolkit import LimWallToolkit

#----------------------------------------------#
# Option 1: Generate a 3D wall file for MAFOT. #
#----------------------------------------------#
stl_path          = "/path/to/file.stl"
output_mafot_file = "mafot_3d_wall.dat"
make_gif          = False
gif_path          = "wall_gif/wall.gif"  # Must create a folder called wall_gif.


#-----------------------------------------------------------------#
# Option 2: Generate .bound file for 3DLIM input from MAFOT file. #
#-----------------------------------------------------------------#

# These files are from running MAFOT in the +1 and -1 directions. Which
# corresponds to file1 and to file2 depends on your setup. One way to make it
# work is file1 can contains all the connection lengths in the positive
# 3DLIM direction (i.e. that which replaces surface 1, L19 in the 3DLIM input
# file), and file2 be the negative values (surface 2, L21).
mafot_file1 = "/Users/zamperini/Documents/d3d_work/184527/mafot/lam_hires_tor240_conn+1.dat"
mafot_file2 = "/Users/zamperini/Documents/d3d_work/184527/mafot/lam_hires_tor240_conn-1.dat"

toroidal_angle = 240  # MiMES = 240.

# The R bins for 3DLIM. Can have them already set or set here and copy/paste
# into 3DLIM input file.
lim_rbins = np.arange(-0.1, 0.02, 0.0025)

# The number of Pbins MUST match 2*MAXNPS+1 in 3DLIM, otherwise it will not
# work. Probably equal to 41. Also I think it needs to be symmetric...
lim_pbins = np.linspace(-0.7, 0.7, 41)

# In machine coordinates, where do you want the origin of 3DLIM to be?
# MiMES-like origin.
r_origin2 = 2.295
z_origin2 = -0.188
# DiMES-like origin.
#r_origin2 = 1.485
#z_origin2 = -1.15

# Which machine coordinate does the 3DLIM "R" coordinate go along?
# MiMES = "R", DiMES = "Z" and None will set it along the actual R-Rsep
# direction (i.e. the plasma radial direction).
along_coord = "R"

# Must have .bound extension, i.e. 184527_conns.bound or ramp.bound are okay.
output_bounds_file = "/Users/zamperini/Documents/lim_runs/184527_tor240.bound"
gfile_pickle_path2 = "/Users/zamperini/Documents/d3d_work/184527/184527_3500.pickle"
wall_path2 = "/Users/zamperini/Documents/d3d_work/184527/mafot/mafot_3D_wall.dat"


#-----------------------------------------------------------------#
# Option 3: Plot 3DLIM results on a 2D cross section of the wall. #
#-----------------------------------------------------------------#
lim_path  = "/Users/zamperini/Documents/lim_runs/actual-184527-tor30.nc"
tor_angle = 30
rmrsep_origin3 = 0.037  # Printed out in option 2.
p_to_z    = [0.9766, -0.188]  # Z = p_to_z[0] * P + p_to_z[1]. Printed out in option 2.
wall_path3 = "/Users/zamperini/Documents/d3d_work/184527/mafot/mafot_3D_wall.dat"

# The gfile is a pickled python dictionary obtained from OMFIT. Follow the
# instructions at the bottom of this script to see how to get it.
plot_surfaces = True
gfile_pickle_path3 = "/Users/zamperini/Documents/d3d_work/184527/184527_3500.pickle"

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
        lwt.generate_3d_wall_for_mafot(stl_path=stl_path, make_gif=make_gif,
          gif_path=gif_path, output_file=output_mafot_file)

    elif int(sys.argv[1]) == 2:
        print("Option 2: Generate a .bound file for 3DLIM from provided MAFOT files with the connection lengths")
        dbg = lwt.bounds_file_from_mafot(tor_angle=toroidal_angle,
          mafot_file1=mafot_file1, mafot_file2=mafot_file2, lim_rbins=lim_rbins,
          lim_pbins=lim_pbins, r_origin=r_origin2,
          z_origin=z_origin2, output_file=output_bounds_file,
          gfile_pickle_path=gfile_pickle_path2, along_coord=along_coord,
          wall_path=wall_path2)

    elif int(sys.argv[1]) == 3:
        print("Option 3: Plot 3DLIM results on a 2D cross section of the wall.")
        lwt.plot_3dlim_on_wall(lim_path=lim_path, tor_angle=tor_angle,
          rmrsep_origin=rmrsep_origin3, p_to_z=p_to_z, wall_path=wall_path3,
          plot_surfaces=plot_surfaces, gfile_pickle_path=gfile_pickle_path3)


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
