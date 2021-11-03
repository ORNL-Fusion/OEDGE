import sys
import math
import pickle
import trimesh
import netCDF4
import pandas                 as pd
import numpy                  as np
import matplotlib.pyplot      as plt
import matplotlib.path        as mplpath
from tqdm                     import tqdm
from collections              import Counter
from shapely.geometry         import Point
from shapely.geometry.polygon import Polygon
from scipy.interpolate        import RBFInterpolator
from matplotlib.colors        import LogNorm, Normalize


class LimWallToolkit:

    def __init__(self):
        pass

    def read_3d_wall(self, wall_path):
        """
        Reads in the 3D wall file that gets passed into MAFOT and returns it
        as a python dictionary, where each key is the toroidal angle. Useful for
        plotting.

        Input
        wall_path (str): Path to the 3D wall file.

        Output
        wall (dict): Dictionary of (R, Z) wall coordinates at each toroidal
          degree. Keys are ints of the degree.
        """

        # Read in the 3D wall coordinates.
        wall = {}
        with open(wall_path) as f:
            count = 0
            tor_angle = 0
            rs = []; zs = []
            for line in f:

                # Skip first three lines of header info and the first # of
                # coordinates.
                if count < 4:
                    count += 1
                    continue

                # If we've hit a single number, this indicates the start of the
                # next set of coordinates, so save what we got and restart
                # arrays.
                ls = line.split()
                if len(ls) == 1:

                    # Add on the first point again so it completes the loop.
                    rs.append(rs[0])
                    zs.append(zs[0])

                    wall[tor_angle] = [rs, zs]
                    rs = []; zs = []
                    tor_angle += 1

                # Save the pair of values in our coordinate lists.
                elif len(ls) == 2:
                    rs.append(float(ls[0]))
                    zs.append(float(ls[1]))
                else:
                    print("Error! Line not recognized: {}".format(ls))
                    break
                count += 1

        # Save the last entry since it won't have a number after it to trigger
        # the saving code.
        wall[tor_angle] = [rs, zs]

        # Add a 360 entry, which is just identical to 0 so things come full
        # circle.
        wall[360] = wall[0]

        return wall

    def generate_3d_wall_for_mafot(self, stl_path, angle_step=1.0,
      make_gif=False, gif_path="wall_gif/wall.gif",
      output_file="mafot_3d_wall.dat"):
        """
        Routine to generate coordinates of each toroidal slice of a supplied
        3D STL wall mesh, formatted for input to MAFOT via the -W flag.

        Input
        stl_path (str): Path to the STL file.
        angle_step (float): Currently only 1.0 is supported due to it being
          hard-coded into MAFOT, but this may change in the future.
        make_gif (bool): Whether or not to generate a gif to see if the point
          sorting worked out okay. Fun.
        output_file (str): Name of the file to be saved to, ready for MAFOT.

        Output
        A file is saved to output_file of the wall coordinates, ready for input
          to MAFOT via the -W flag.
        """

        # Sort points poloidally, clockwise.
        def sort_polar(center, rs, zs, poly):

            # Pop the entries that fall within the polygon. These are the ones
            # we want to sort.
            points = [Point(rs[i], zs[i]) for i in range(0, len(rs))]
            inside = np.array([poly.contains(p) for p in points])

            # The Points object consists of an __array__ attribute, and numpy
            # does not like creating an array of arrays that may be different
            # shapes (i.e. np.array(points)), so the suggested fix is to
            # created an empty array of type object and then fill it.
            points_np = np.empty(len(points), dtype=object)
            for i in range(0, len(points)):
                points_np[i] = points[i]

            # Create a new array that we will put the sorted points into,
            # starting at first_idx.
            new_points = points_np[~inside]
            sort_points = points_np[inside]
            first_idx = np.argwhere(inside==True)[0]

            # Sort the points according to polar angle. We could in theory use
            # the centroid but for greater flexibility best to just allow the
            # user defined center coordinate. 2*pi - angle to go clockwise since
            # that just feels natural.
            pol_ang = [2 * np.pi - math.atan2(p.y - center[1], p.x - center[0]) for p in sort_points]
            sort_idx = np.argsort(pol_ang)

            # Insert into the new_points array at the first index where the
            # sorted points were pulled from.
            new_points = np.insert(new_points, first_idx, sort_points[sort_idx])

            # Break back down into r and z for output.
            rs = [p.x for p in new_points]
            zs = [p.y for p in new_points]

            return rs, zs

        # First load in the STL file with a soft exit if failed.
        try:
            mesh = trimesh.load_mesh(stl_path)
            print("Loaded mesh file: {}".format(stl_path.split("/")[-1]))
        except ValueError:
            print("Error! STL file not found: {}".format(stl_path))
            return None

        # Require watertightness.
        if not mesh.is_watertight:
            print("Error! Mesh must be watertight.")
            return None

        # Polygons defining regions to sort and re-sort points.
        polys = []
        polys.append(Polygon([(-99, -99),    (99,   -99),   (99,   99),    (-99,  99)]))
        polys.append(Polygon([(2.00, 0.35),  (2.60, 0.35),  (2.60, 1.10),  (2.00, 1.10)]))
        polys.append(Polygon([(1.80, -1.10), (2.60, -1.10), (2.60, -0.35), (1.80, -0.35)]))
        polys.append(Polygon([(1.45, 1.08),  (1.45, 1.30),  (1.55, 1.30),  (1.55, 1.08)]))
        polys.append(Polygon([(1.20, 1.12),  (1.38, 1.14),  (1.47, 1.38),  (1.20, 1.38)]))
        polys.append(Polygon([(1.32, 1.28),  (1.32, 1.38),  (1.44, 1.38),  (1.44, 1.28)]))
        polys.append(Polygon([(1.69, 0.95),  (1.69, 1.33),  (2.21, 1.33),  (2.21, 0.95)]))

        # We want to take slices of the 3D volume to get a successive sequence of
        # 2D cross-sections, one at each degree. We do this by giving trimesh the
        # surface normal of the intersecting plane.
        sections = []
        degrees = np.arange(0.0, 360.0, angle_step)
        print("Generating cross sections...")
        for i in tqdm(range(0, len(degrees))):
            deg = degrees[i]

            # We add 90 since we want the surface normal at this angle (which is
            # 90 degrees away).
            x_norm = np.cos(np.radians(deg + 90))
            y_norm = np.sin(np.radians(deg + 90))

            # Can use this to take slices at each degree, returning the X, Y
            # coordinates.
            slice = trimesh.intersections.mesh_plane(mesh, plane_origin=[0,0,0],
              plane_normal=[x_norm, y_norm, 0])

            # Slice is a sequence of 3D lines, where each entry is the start and
            # and coordinates of the line [(X0, Y0, Z0), (X1, Y1, Z1)]. Each line
            # segment is connected, so each coordinate is repeated for the two
            # line segements it shares, but it's not consistent as to if the 1's
            # connect to the next segment's 0's or 1's, so we will go through and
            # pull out one coordinate at a time, avoiding duplicates.
            points = []
            for pointset in slice:
                for point in pointset:
                    point = list(point)
                    if point not in points:
                        points.append(point)
            points = np.array(points)

            # At this point can define R and Z, in m.
            r = np.sqrt(np.square(points[:,0]) + np.square(points[:,1])) / 1000
            z = points[:,2] / 1000

            # Keep the correct half of the cross-section.
            if deg > 0 and deg <= 180:
                keep = points[:,0] < 0
            else:
                keep = points[:,0] > 0
            x_cs = points[:,0][keep]
            y_cs = points[:,1][keep]
            r = r[keep]
            z = z[keep]

            # Sort using the middle of the vessel as a center. This doesn't get all
            # of them in the correct order, but it's a start.
            # Open thought: Could a sorting algorithm be made such that you divide
            # a region into subregions, and then sort points in a region clockwise,
            # where you change the number of regions as a knob? Hmmmm...
            r, z = sort_polar((1.1, 0.0), r, z, polys[0])

            # Re-sort the port locations. The middle port seems to get sorted just
            # fine via the first sort. This sorts the ports above and below the
            # midplane.
            r, z = sort_polar((2.00, 0.70), r, z, polys[1])
            r, z = sort_polar((2.00, -0.60), r, z, polys[2])

            # Re-sort the SAS region.
            r, z = sort_polar((1.44, 1.10), r, z, polys[3])

            # Re-sort the upper divertor area.
            r, z = sort_polar((1.30, 1.14), r, z, polys[4])

            # Re-sort the UOB area.
            r, z = sort_polar((1.32, 1.32), r, z, polys[5])

            # Re-sort this top port area.
            r, z = sort_polar((1.95, 1.00), r, z, polys[6])

            # Machine coordinates are clockwise so 360 - phi.
            phi_mach = 360 - deg

            sections.append({"degree":deg, "x_norm":x_norm, "y_norm":y_norm,
              "x_cs":x_cs, "y_cs":y_cs, "r":r, "z":z, "phi_mach":phi_mach})

        # Create the MAFOT input wall file. First see what the maximum number of
        # points are for a wall section.
        max_pts = 0
        for i in range(0, len(sections)):
            pts = len(sections[i]["r"])
            if pts > max_pts:
                max_pts = pts

        # Open file and print the header info.
        with open(output_file, "w") as f:
            f.write("# DIII-D 3D wall: (R,Z) for every degree: phi = 0,...,359;" + \
              " leading integer gives number of points per plane\n")
            f.write("# The first integer gives the maximum number of points in" + \
              " any plane\n")
            f.write(str(max_pts) + "\n")

            # Then write the R, Z coordinates with the number of points first.
            for i in range(0, len(sections)):
                num_pts = len(sections[i]["r"])
                f.write(str(num_pts) + "\n")
                for j in range(0, num_pts):
                    f.write("{:7.5f} {:7.5f}\n".format(sections[i]["r"][j], sections[i]["z"][j]))

        # Optional flag to make a gif of the wall as you progress around the vessel.
        if make_gif:
            import imageio
            import os

            fnames = []
            print("Creating gif...")
            for i in tqdm(range(0, len(sections))):
                fig, ax = plt.subplots()
                ax.plot(sections[i]["r"], sections[i]["z"], color="r", zorder=9)
                ax.scatter(sections[i]["r"], sections[i]["z"], color="k", zorder=10, s=5)
                for poly in polys[1:]:
                    ax.plot(*poly.exterior.xy, color="b")

                ax.set_aspect("equal")
                ax.text(2.0, 1.3, r"$\phi$={:}$^\circ$".format(sections[i]["phi_mach"]))
                ax.set_xlim([0.95, 2.75])
                ax.set_ylim(-1.5, 1.5)

                fname = "{}/{}.png".format(gif_path.split("/")[0], i)
                #fname = "wall_plots/{}.png".format(i)
                fnames.append(fname)
                fig.savefig(fname)
                plt.close(fig)

            # Now build the gif and clean up by removing the plots.
            print("Saving as gif...")
            with imageio.get_writer(gif_path, mode="I") as writer:
                for fname in fnames:
                    image = imageio.imread(fname)
                    writer.append_data(image)
            for fname in fnames:
                os.remove(fname)


    def bounds_file_from_mafot_old(self, tor_angle, mafot_file1, mafot_file2,
      lim_rbins, lim_pbins, r_origin, z_origin, output_file,
      gfile_pickle_path, show_plot=True):
        """
        DEFUNCT. DON'T USE ME.
        This function generates a .bound file, which is used as input to 3DLIM.
        This is so that a MAFOT run calculating the connection lengths
        with a realistic 3D wall geometry can be formatted for use in 3DLIM.
        Some approximations are made in this function, as the rectilinear
        approximation in 3DLIM is getting pushed to the limits at this point.

        Input
        tor_angle (float): The toroidal angle at which this 3DLIM simulation
          is to be ran for. I.e. if we are doing a collector probe on MiMES,
          then we would enter 240.
        mafot_file1 (str): The MAFOT file run in the direction of what we are
          calling the positive direction in 3DLIM (either +1 or -1, depends
          on how we set the problem up).
        mafot_file2 (str): The other direction (+1 or -1).
        lim_rbins (list/array): The R bins for the 3DLIM input file. You could
          either have the bins predetermined, or just set them here and copy/
          paste them into the input file afterwards.
        lim_pbins (list/array): Similar, but the bins for the poloidal (or more
          accurately, perpendicular) direction. WARNING: The number of bins
          here matters due to 3DLIM code limitations. It must equal 2*MAXNPS+1,
          which is commonly equal to 41. If you do not do this then 3DLIM will
          not work correctly with the resulting .bound file.
        rmrsep_origin (float): The R-Rsep value, of the origin of the 3DLIM
          volume. I.e. if we are simulating a collector probe, then rmrsep_origin
          would be the R-Rsep value of the tip of the probe.
        z_origin (float): Likewise, but the Z value.
        output_file (str):
        gfile_pickle_path (str):
        show_plot (bool):
        """

        # Useful to check some of the data.
        debug_dict = {}

        # Load data into DataFrame.
        print("Loading MAFOT runs...")
        columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
          "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
        try:
            df1 = pd.read_csv(mafot_file1, skiprows=52, names=columns,
              delimiter="\t")
        except FileNotFoundError:
            print("Error: Unable to find file: {}".format(mafot_file1))
            print("Exiting")
            sys.exit()
        try:
            df2 = pd.read_csv(mafot_file2, skiprows=52, names=columns,
              delimiter="\t")
        except FileNotFoundError:
            print("Error: Unable to find file: {}".format(mafot_file2))
            print("Exiting")
            sys.exit()

        # Also read the file to pull out the number of R and Z coords.
        with open(mafot_file1) as f:
            for line in f:
                if line[:10] == "# R-grid: ":
                    numrs1 = int(line.split(":")[1])
                if line[:10] == "# Z-grid: ":
                    numzs1 = int(line.split(":")[1])
                    break
        with open(mafot_file2) as f:
            for line in f:
                if line[:10] == "# R-grid: ":
                    numrs2 = int(line.split(":")[1])
                if line[:10] == "# Z-grid: ":
                    numzs2 = int(line.split(":")[1])
                    break

        # Reshape into 2D arrays. Reasonable assumption that both use the
        # same R, Z.
        if df1.shape != df2.shape:
            print("Error: Please use the same number of R and Z coordinates" + \
              " for both MAFOT runs.")
            #numrs1 = len(df1["R (m)"].unique())
            #numrs2 = len(df2["R (m)"].unique())
            #numzs1 = len(df1["Z (m)"].unique())
            #numzs2 = len(df2["Z (m)"].unique())
            print("            | # R's | # Z's |")
            print("mafot_file1 | {:5} | {:5} |".format(numrs1, numzs1))
            print("mafot_file2 | {:5} | {:5} |".format(numrs2, numzs2))
            print("Exiting")
            sys.exit()
        r = df1["R (m)"].unique()
        z = df1["Z (m)"].unique()
        #r = df1["R (m)"][:numrs1]
        #z = df1["Z (m)"].unique()
        R, Z = np.meshgrid(r, z)

        # Pull out each's connection length and pitch angles.
        l1 = df1["Lconn (km)"].values * 1000  # km to m
        l2 = df2["Lconn (km)"].values * 1000
        p1 = df1["pitch angle"].values
        p2 = df2["pitch angle"].values
        L1 = l1.reshape(len(r), len(z))
        L2 = l2.reshape(len(r), len(z))
        pitch1 = p1.reshape(len(r), len(z))
        pitch2 = p2.reshape(len(r), len(z))

        # G-file needed for the equilibrium related info. This is not the actual
        # gfile, but rather a pickled version from the OMFIT EFIT module. See
        # comment at bottom of lwt_control_file.py file for instructions.
        with open(gfile_pickle_path, "rb") as f:
            gfile = pickle.load(f)
        gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
        psin = gfile["PSIRZ_NORM"]
        R_sep = gfile["RBBBS"]
        Z_sep = gfile["ZBBBS"]
        debug_dict["gfile"] = gfile

        # Go through one location at a time and find out how far it is from the
        # nearest point on the separatrix (R-Rsep).
        print("Calculating distances...")
        dists = np.zeros(gR.shape)
        for i in range(0, gR.shape[0]):
            for j in range(0, gR.shape[1]):
                d = np.sqrt(np.square(gR[i, j] - R_sep) + np.square(gZ[i, j] - Z_sep))
                if psin[i, j] < 1:
                    dists[i, j] = -d.min()
                else:
                    dists[i, j] = d.min()
        debug_dict["gR"] = gR
        debug_dict["gZ"] = gZ
        debug_dict["dists"] = dists

        # Given the (R, Z) in machine coordinates of the 3DLIM origin, find out
        # what this R-Rsep value is with an interpolation function.
        print("Finding 3DLIM origin R-Rsep...")
        int_R = gR.flatten()[::10]
        int_Z = gZ.flatten()[::10]
        int_dists = dists.flatten()[::10]
        coords = list(zip(int_R, int_Z))
        f_dist = RBFInterpolator(coords, int_dists, smoothing=10)
        rmrsep_origin = f_dist(np.array([r_origin, z_origin]).reshape(-1, 2))[0]
        print("R-Rsep 3DLIM origin is: {:.3f}".format(rmrsep_origin))

        # 3DLIM R-coordinate is actually ~R-Rsep, just off by some constant
        # that depends on how far away our simulation volume is from the
        # separatrix. This means
        rmrseps = [rmrsep_origin - rb for rb in lim_rbins]

        # The corresponding Z values from the 3DLIM P (for perpendicular) bins
        # are a bit more work since the P direction is not exactly in the Z
        # direction and is angled off by the pitch angle of the field line. This
        # means the P direction does not stay at a single toroidal angle, and
        # instead slightly jumps across different angles.
        # Here we use the pitch angle at the origin. We are approximating that
        # the pitch angle does not change along the field line, otherwise a
        # rectangular grid like 3DLIM does not work. Square into a circle kinda
        # thing.
        # LIMITATION: This assumes that lim_rbins are essentially horizontal,
        # i.e. pretty much at the outboard midplane.
        dist = np.sqrt(np.square(R - r_origin) + np.square(Z - z_origin))
        nearest_idx = np.where(dist == np.nanmin(dist))
        origin_pitch = pitch1[nearest_idx]
        zs = [float(z_origin + pb * np.cos(origin_pitch)) for pb in lim_pbins]
        print("pitch1 = {:.3f}     pitch2 = {:.3f}".format(pitch1[nearest_idx][0],
          pitch2[nearest_idx][0]))
        debug_dict["rmrseps"] = rmrseps
        debug_dict["zs"] = zs

        # Do a linear fit so in later scripts we can convert between lim_pbins
        # and Z. Returned coefficients here are highest power first.
        pfit = np.polyfit(lim_pbins, zs, 1)

        # Make an interpolation function R(RMRS, Z) so that we can see what the
        # corresponding machine R value for each 3DLIM bin is. Cut off
        # everything to the left of the plasma center, otherwise it's an ill-
        # defined problem using distance as an input (fails the 2D version of
        # the straight line test). Use every 10 points too to prevent memory
        # from getting out of hand.
        print("Interpolating 3DLIM coordinates onto grid...")
        print("  Note: Using every 10th EFIT data point")
        right = gR > gfile["RMAXIS"]
        int_R = gR[right].flatten()[::10]
        int_Z = gZ[right].flatten()[::10]
        int_dists = dists[right].flatten()[::10]
        coords = list(zip(int_dists, int_Z))
        f_R = RBFInterpolator(coords, int_R, smoothing=10)
        rmrseps_2d, zs_2d = np.meshgrid(rmrseps, zs)
        R_lim_2d = f_R(list(zip(rmrseps_2d.flatten(), zs_2d.flatten()))).reshape(zs_2d.shape)
        debug_dict["f_R"] = f_R
        debug_dict["zs_2d"] = zs_2d
        debug_dict["R_lim_2d"] = R_lim_2d

        # Create empty bounds array and fill it with the closest connection
        # length. The approximation above, where it is assumed the pitch angle
        # does not change across toroidal angles is still made here. Got a
        # little confused with the array dimensions here, but I think it's
        # right.
        print("Determining bounds...")
        #bounds1 = np.zeros((len(rmrseps), len(zs)))
        #bounds2 = np.zeros((len(rmrseps), len(zs)))
        bounds1 = np.zeros(R_lim_2d.shape).T
        bounds2 = np.zeros(R_lim_2d.shape).T
        for ir in range(0, R_lim_2d.shape[0]):
            for iz in range(0, R_lim_2d.shape[1]):
                dist = np.sqrt(np.square(R - R_lim_2d[ir][iz]) + np.square(Z - zs_2d[ir][iz]))
                nearest_idx = np.where(dist == np.nanmin(dist))
                #bounds1[ir][iz] = L1[nearest_idx]
                #bounds2[ir][iz] = -L2[nearest_idx]
                bounds1[iz][ir] = L1[nearest_idx]
                bounds2[iz][ir] = -L2[nearest_idx]
        debug_dict["bounds1"] = bounds1
        debug_dict["bounds2"] = bounds2

        # Print to a text file ready for input to 3DLIM.
        with open(output_file, "w") as f:
            print("Writing to {}".format(output_file))
            f.write("Absorbing boundary data. Rows are for each radial bin, columns are for each poloidal bin. Created via LimPlotToolkit.bounds_file_from_mafot.py.\n")
            f.write("Dimensions: {} {}\n".format(len(rmrseps), len(zs)))
            f.write("yabsorb1a boundary (all values should be positive):\n")
            for ir in range(0, len(rmrseps)):
                f.write("{:5.2f}".format(bounds1[ir][0]) + ' ' + ' '.join(["{:5.2f}".format(b) for b in bounds1[ir][1:]]) + '\n')
            f.write("yabsorb2a boundary (all values should be negative):\n")
            for ir in range(0, len(rmrseps)):
                f.write("{:5.2f}".format(bounds2[ir][0]) + ' ' + ' '.join(["{:5.2f}".format(b) for b in bounds2[ir][1:]]) + '\n')

        # Print out the bin values to be copy/pasted to input file.
        print("R bins for copy/paste to 3DLIM input file:")
        print("Number of bins: {}".format(len(lim_rbins)))
        for r in lim_rbins:
            print("{:.3f}".format(r))
        print()
        print("P bins for copy/paste to 3DLIM input file:")
        print("Number of bins: {}".format(len(lim_pbins)))
        for p in lim_pbins:
            print("{:.3f}".format(p))
        print()
        nybins = 150
        print("Option for Y Bins")

        # Don't include zero as 3DLIM requires you don't include it.
        print("Number of bins: {}".format(nybins-1))
        lim_ybins = np.linspace(0, max(bounds1.max(), bounds2.max()), nybins)
        for y in lim_ybins[1:]:
            print("{:.2f}".format(y))
        print()
        print("Some additional comments")
        print("  - Change AW to something less than {:.3f}".format(lim_rbins[0]))
        print("  - Change A to something greater than {:.3f}".format(lim_rbins[-1]))
        print("  - Maximum bounds value is {:.3f}".format(max(bounds1.max(), bounds2.max())))
        print("These values are needed for option 3. Write them down!")
        print("  - Z = {:.4f} * lim_pbins + {:.4f}".format(pfit[0], pfit[1]))
        print("  - R-Rsep 3DLIM origin = {:.3f}".format(rmrsep_origin))

        if show_plot:

            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

            # Negative bounds2 just so we can plot the relative magnitudes
            # of each direction. The negative is just needed for the .bound file
            # which is already done.
            vmin = min(bounds1.min(), -bounds2.min())
            vmax = max(bounds1.max(), -bounds2.max())
            cont1 = ax1.pcolormesh(zs, rmrseps, bounds1, shading="auto", vmin=vmin,
              vmax=vmax, cmap="inferno")
            cont2 = ax2.pcolormesh(zs, rmrseps, -bounds2, shading="auto", vmin=vmin,
              vmax=vmax, cmap="inferno")

            cbar = fig.colorbar(cont2, ax=[ax1, ax2], location="right")
            #cbar.set_label("Distance from 3DLIM origin", fontsize=12)

            ax1.set_title("Positive Bounds")
            ax2.set_title("Negative Bounds")
            ax1.set_xlabel("Z (m)", fontsize=12)
            ax1.set_ylabel("R (m)", fontsize=12)
            ax2.set_xlabel("Z (m)", fontsize=12)
            #fig.tight_layout()

            total_bounds = bounds1 - bounds2  # Minuc bc bounds2 is negative.
            cont3 = ax3.pcolormesh(zs, rmrseps, total_bounds, shading="auto",
              vmin=total_bounds.min(), vmax=total_bounds.max(), cmap="inferno")
            ax3.set_title("Total Field Line Length")
            ax3.set_xlabel("Z (m)")
            cbar = fig.colorbar(cont3, ax=ax3)
            plt.show()

            # Another plot with the 3DLIM points mapped to the equlibrium.
            fig, ax = plt.subplots()
            cont = ax.pcolormesh(gR, gZ, dists, shading="auto", cmap="coolwarm",
              vmin=-1, vmax=1)
            cbar = fig.colorbar(cont)
            ax.contour(gR, gZ, dists, levels=[0], colors="k")
            ax.scatter(R_lim_2d.flatten(), zs_2d.flatten(), s=3, color="k")
            ax.set_aspect("equal")
            ax.spines["top"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            fig.tight_layout()
            fig.show()

        return debug_dict

    def get_mach_coords(self, gR, gZ, psin, mid_r, mid_z, lim_pbins,
      along_coord=None):
        """
        Using a single machine coordinate identifying the middle of the poloidal
        bins (mid_r, mid_z), return arrays of len(lim_pbins) containing the
        rest of the machine coordinates (R, Z) for each poloidal bin.

        Input
        gR (array): Machine R coordinates from gfile.
        gZ (array): Machine Z coordinates from gfile.
        psin (array): Normalized psi values at each gR, gZ.
        mid_r (float): The machine R coordinate of where Pbin = 0.
        mid_z (float): The machine Z coordinate of where Pbin = 0.
        lim_pbins (array): The lim Pbins under consideration.
        along_coord (str): One of "R", "Z" or None. This is choosing along which
          coordinate direction the 3DLIM R bins are along. More practically,
          at MiMES this would be "R", and at DiMES this would be "Z". If None,
          then the 3DLIM R coordinate is what it naturally would be, i.e. the
          coordinate would be along the true R-Rsep direction, which is always
          perpendicular to the separatrix.

        Output
        {lim_machR, lim_machZ} (dict): The machine coordinates at each Pbin
          location along the field line that intersects (mid_r, mid_z).
        """


        # Create 2D array of distance from input coordinate and find where
        # the distance is a minimum to find our psin we care about.
        dist = np.sqrt(np.square(mid_r-gR) + np.square(mid_z-gZ))
        psin_mid = psin[np.where(dist == dist.min())]

        fig, ax = plt.subplots()

        # Bit of hardcoding to ensure we pick the right region of flux tubes.
        # Have only tested for USN MiMES and DiMES so far.
        if type(along_coord) == type(None):
            print("Error: along_coord None not implemented yet.")

        # Essentially hard-coding the MiMES area here.
        elif along_coord == "R":
            keep = gR > 1.50
            keep = keep[0,:]
            gR_keep = gR[:,keep]
            gZ_keep = gZ[:,keep]
            psin_keep = psin[:,keep]

        # Essentially hard-coding the DiMES area here.
        elif along_coord == "Z":
            keep = gZ < 0.00
            keep = keep[:,0]
            gR_keep = gR[keep]
            gZ_keep = gZ[keep]
            psin_keep = psin[keep]
        else:
            print("Error: along_coord \"{}\" not recognized".format(along_coord))
            sys.exit()

        # Able to hijack the contour function to grab our field line coordinates.
        cont2 = ax.contour(gR_keep, gZ_keep, psin_keep, levels=[psin_mid],
          colors="r")
        mid_line_rz = cont2.allsegs[0][0]

        # Now calculate distance along our field line (which is in the poloidal
        # direction).
        dist = np.sqrt(np.square(mid_line_rz[:,0] - mid_r) +
          np.square(mid_line_rz[:,1] - mid_z))
        close_mid_idx = np.where(dist == dist.min())[0][0]
        close_mid = mid_line_rz[close_mid_idx]

        # To hold the corresponding R, Z coordinate for each P bin.
        lim_machR = np.zeros(lim_pbins.shape)
        lim_machZ = np.zeros(lim_pbins.shape)

        pbin_center_idx = np.where(lim_pbins == 0)[0][0]

        # Go in one direction from the origin along the field line.
        pbin_idx = pbin_center_idx
        tot_dist = 0
        prev_coord = close_mid
        for d in mid_line_rz[close_mid_idx:]:

            # Tally distance travelled along field line, assigning an R, Z
            # coordinate each time we pass a P bin.
            tmp_dist = np.sqrt(np.square(d[0] - prev_coord[0]) + np.square(d[1]
              - prev_coord[1]))
            tot_dist += tmp_dist
            if tot_dist >= np.abs(lim_pbins[pbin_idx]):
                lim_machR[pbin_idx] = d[0]
                lim_machZ[pbin_idx] = d[1]
                pbin_idx += 1
                if pbin_idx > len(lim_pbins) - 1:
                    break
            prev_coord = d

        # Likewise in the other direction.
        pbin_idx = pbin_center_idx
        tot_dist = 0
        prev_coord = close_mid
        for d in mid_line_rz[:close_mid_idx][::-1]:
            tmp_dist = np.sqrt(np.square(d[0] - prev_coord[0]) + np.square(d[1]
              - prev_coord[1]))
            tot_dist += tmp_dist
            if tot_dist >= np.abs(lim_pbins[pbin_idx]):
                lim_machR[pbin_idx] = d[0]
                lim_machZ[pbin_idx] = d[1]
                pbin_idx -= 1
                if pbin_idx < 0:
                    break
            prev_coord = d

        plt.close()
        return {"lim_machR":lim_machR, "lim_machZ":lim_machZ}

    def bounds_file_from_mafot(self, tor_angle, mafot_file1, mafot_file2,
      lim_rbins, lim_pbins, r_origin, z_origin, output_file,
      gfile_pickle_path, show_plot=True, along_coord=None, wall_path=None):
        """
        This function generates a .bound file, which is used as input to 3DLIM.
        This is so that a MAFOT run calculating the connection lengths
        with a realistic 3D wall geometry can be formatted for use in 3DLIM.
        Some approximations are made in this function, as the rectilinear
        approximation in 3DLIM is getting pushed to the limits at this point.
        This assumes the the perpendicular direction in 3DLIM is actually the
        poloidal (is isn't actaully, it's off by the pitch angle of the field
        lines), and that there is no flux expansion along the field lines. This
        is not true either, and it is unclear right now what it means.

        Input
        tor_angle (float): The toroidal angle at which this 3DLIM simulation
          is to be ran for. I.e. if we are doing a collector probe on MiMES,
          then we would enter 240.
        mafot_file1 (str): The MAFOT file run in the direction of what we are
          calling the positive direction in 3DLIM (either +1 or -1, depends
          on how we set the problem up).
        mafot_file2 (str): The other direction (+1 or -1).
        lim_rbins (list/array): The R bins for the 3DLIM input file. You could
          either have the bins predetermined, or just set them here and copy/
          paste them into the input file afterwards.
        lim_pbins (list/array): Similar, but the bins for the poloidal (or more
          accurately, perpendicular) direction. WARNING: The number of bins
          here matters due to 3DLIM code limitations. It must equal 2*MAXNPS+1,
          which is commonly equal to 41. If you do not do this then 3DLIM will
          not work correctly with the resulting .bound file.
        rmrsep_origin (float): The R-Rsep value, of the origin of the 3DLIM
          volume. I.e. if we are simulating a collector probe, then rmrsep_origin
          would be the R-Rsep value of the tip of the probe.
        z_origin (float): Likewise, but the Z value.
        output_file (str): Path to where to save the .bound file.
        gfile_pickle_path (str): Path to a pickled gfile created from the script
          at the bottom of lwt_control_file.py.
        show_plot (bool): Whether or not to show a plot vizualizing what we've
          done. Generally leave as True.
        along_coord (str): One of "R", "Z" or None. This is choosing along which
          coordinate direction the 3DLIM R bins are along. More practically,
          at MiMES this would be "R", and at DiMES this would be "Z". If None,
          then the 3DLIM R coordinate is what it naturally would be, i.e. the
          coordinate would be along the true R-Rsep direction, which is always
          perpendicular to the separatrix.
        wall_path (str): Optional. Path to the 3D wall file (probably called
          mafot_3D_wall.dat)

        Output
        debug_dict (dict): A dictionary with some of the data used throughout
          the function for debugging purposes.
        """

        # Useful to check some of the data.
        debug_dict = {}

        # Load data into DataFrame.
        print("Loading MAFOT runs...")
        columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
          "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
        try:
            df1 = pd.read_csv(mafot_file1, skiprows=52, names=columns,
              delimiter="\t")
        except FileNotFoundError:
            print("Error: Unable to find file: {}".format(mafot_file1))
            print("Exiting")
            sys.exit()
        try:
            df2 = pd.read_csv(mafot_file2, skiprows=52, names=columns,
              delimiter="\t")
        except FileNotFoundError:
            print("Error: Unable to find file: {}".format(mafot_file2))
            print("Exiting")
            sys.exit()

        # Also read the file to pull out the number of R and Z coords.
        with open(mafot_file1) as f:
            for line in f:
                if line[:10] == "# R-grid: ":
                    numrs1 = int(line.split(":")[1])
                if line[:10] == "# Z-grid: ":
                    numzs1 = int(line.split(":")[1])
                    break
        with open(mafot_file2) as f:
            for line in f:
                if line[:10] == "# R-grid: ":
                    numrs2 = int(line.split(":")[1])
                if line[:10] == "# Z-grid: ":
                    numzs2 = int(line.split(":")[1])
                    break

        # Reshape into 2D arrays. Reasonable assumption that both use the
        # same R, Z.
        if df1.shape != df2.shape:
            print("Error: Please use the same number of R and Z coordinates" + \
              " for both MAFOT runs.")
            #numrs1 = len(df1["R (m)"].unique())
            #numrs2 = len(df2["R (m)"].unique())
            #numzs1 = len(df1["Z (m)"].unique())
            #numzs2 = len(df2["Z (m)"].unique())
            print("            | # R's | # Z's |")
            print("mafot_file1 | {:5} | {:5} |".format(numrs1, numzs1))
            print("mafot_file2 | {:5} | {:5} |".format(numrs2, numzs2))
            print("Exiting")
            sys.exit()
        r = df1["R (m)"].unique()
        z = df1["Z (m)"].unique()
        #r = df1["R (m)"][:numrs1]
        #z = df1["Z (m)"].unique()
        R, Z = np.meshgrid(r, z)

        # Pull out each's connection length and pitch angles.
        l1 = df1["Lconn (km)"].values * 1000  # km to m
        l2 = df2["Lconn (km)"].values * 1000
        p1 = df1["pitch angle"].values
        p2 = df2["pitch angle"].values
        L1 = l1.reshape(len(r), len(z))
        L2 = l2.reshape(len(r), len(z))
        pitch1 = p1.reshape(len(r), len(z))
        pitch2 = p2.reshape(len(r), len(z))

        # G-file needed for the equilibrium related info. This is not the actual
        # gfile, but rather a pickled version from the OMFIT EFIT module. See
        # comment at bottom of lwt_control_file.py file for instructions.
        with open(gfile_pickle_path, "rb") as f:
            gfile = pickle.load(f)
        gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
        psin = gfile["PSIRZ_NORM"]
        R_sep = gfile["RBBBS"]
        Z_sep = gfile["ZBBBS"]
        debug_dict["gfile"] = gfile

        lim_machRs = np.zeros((len(lim_rbins), len(lim_pbins)))
        lim_machZs = np.zeros((len(lim_rbins), len(lim_pbins)))
        for i in range(0, len(lim_rbins)):
            lim_r = lim_rbins[i]
            if type(along_coord) == type(None):
                print("Error: along_coord None not implemented yet.")
                sys.exit()
            elif along_coord == "R":
                # Get the corresponding machine R coordinate.
                machR = r_origin - lim_r
                machZ = z_origin
            elif along_coord == "Z":
                # If DiMES, then the 3DLIM R coord is actually parallel to Z.
                machZ = z_origin - lim_r
                machR = r_origin
            else:
                print("Error: along_coord \"{}\" not recognized".format(along_coord))
                sys.exit()

            mach_coords = self.get_mach_coords(gR, gZ, psin, machR, machZ,
              lim_pbins, along_coord)
            lim_machRs[i] = mach_coords["lim_machR"]
            lim_machZs[i] = mach_coords["lim_machZ"]
        debug_dict["lim_machRs"] = lim_machRs
        debug_dict["lim_machZs"] = lim_machZs

        # Issue a warning if the 3DLIM coordinates once mapped to the machine
        # have coordinates outside the range in which MAFOT was run.
        warn = False
        if lim_machRs.min() < r.min():
            warn = True
        elif lim_machRs.max() > r.max():
            warn = True
        elif lim_machZs.min() < z.min():
            warn = True
        elif lim_machZs.max() > z.max():
            warn = True
        if warn:
            print("Warning! 3DLIM coordinates are outside of the MAFOT grid!")
            print("       |        R       |        Z       |")
            print(" MAFOT | [{:>5.2f}, {:>5.2f}] | [{:>5.2f}, {:>5.2f}] |".format(
              r.min(), r.max(), z.min(), z.max()))
            print(" 3DLIM | [{:>5.2f}, {:>5.2f}] | [{:>5.2f}, {:>5.2f}] |".format(
              lim_machRs.min(), lim_machRs.max(), lim_machZs.min(),
              lim_machZs.max()))


        # Create empty bounds array and fill it with the closest connection
        # length. The approximation above, where it is assumed the pitch angle
        # does not change across toroidal angles is still made here. Got a
        # little confused with the array dimensions here, but I think it's
        # right.
        print("Determining bounds...")
        #bounds1 = np.zeros((len(rmrseps), len(zs)))
        #bounds2 = np.zeros((len(rmrseps), len(zs)))
        bounds1 = np.zeros(lim_machRs.shape)
        bounds2 = np.zeros(lim_machRs.shape)
        for ir in range(0, lim_machRs.shape[0]):
            for iz in range(0, lim_machRs.shape[1]):
                dist = np.sqrt(np.square(R - lim_machRs[ir][iz]) + np.square(Z - lim_machZs[ir][iz]))
                nearest_idx = np.where(dist == np.nanmin(dist))
                #print(len(nearest_idx))
                #if len(nearest_idx) > 1:
                #    i = nearest_idx[0][0]
                #    k = nearest_idx[0][1]
                #    nearest_idx = [[i], [k]]
                #bounds1[ir][iz] = L1[nearest_idx]
                #bounds2[ir][iz] = -L2[nearest_idx]
                bounds1[ir][iz] = L1[nearest_idx]
                bounds2[ir][iz] = -L2[nearest_idx]
        debug_dict["bounds1"] = bounds1
        debug_dict["bounds2"] = bounds2

        # Print to a text file ready for input to 3DLIM.
        with open(output_file, "w") as f:
            print("Writing to {}".format(output_file))
            f.write("Absorbing boundary data. Rows are for each radial bin, columns are for each poloidal bin. Created via LimPlotToolkit.bounds_file_from_mafot.py.\n")
            f.write("Dimensions: {} {}\n".format(len(lim_rbins), len(lim_pbins)))
            f.write("yabsorb1a boundary (all values should be positive):\n")
            for ir in range(0, len(lim_rbins)):
                f.write("{:5.2f}".format(bounds1[ir][0]) + ' ' + ' '.join(["{:5.2f}".format(b) for b in bounds1[ir][1:]]) + '\n')
            f.write("yabsorb2a boundary (all values should be negative):\n")
            for ir in range(0, len(lim_rbins)):
                f.write("{:5.2f}".format(bounds2[ir][0]) + ' ' + ' '.join(["{:5.2f}".format(b) for b in bounds2[ir][1:]]) + '\n')

        # Print out the bin values to be copy/pasted to input file.
        print("R bins for copy/paste to 3DLIM input file:")
        print("Number of bins: {}".format(len(lim_rbins)))
        for r in lim_rbins:
            print("{:.3f}".format(r))
        print()
        print("P bins for copy/paste to 3DLIM input file:")
        print("Number of bins: {}".format(len(lim_pbins)))
        for p in lim_pbins:
            print("{:.3f}".format(p))
        print()
        nybins = 150
        print("Option for Y Bins")

        # Don't include zero as 3DLIM requires you don't include it.
        print("Number of bins: {}".format(nybins-1))
        lim_ybins = np.linspace(0, max(bounds1.max(), bounds2.max()), nybins)
        for y in lim_ybins[1:]:
            print("{:.2f}".format(y))
        print()
        print("Some additional comments")
        print("  - Change AW to something less than {:.3f}".format(lim_rbins[0]))
        print("  - Change A to something greater than {:.3f}".format(lim_rbins[-1]))
        print("  - Maximum bounds value is {:.3f}".format(max(bounds1.max(), bounds2.max())))
        print("These values are needed for option 3. Write them down!")
        print("  - R-Rsep 3DLIM origin = {:.3f}".format(r_origin))
        if warn:
            print()
            print("***********************************************************")
            print("* ERROR: MAFOT bounds exceeded! Scroll up to check error. *")
            print("***********************************************************")

        if show_plot:

            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

            # Negative bounds2 just so we can plot the relative magnitudes
            # of each direction. The negative is just needed for the .bound file
            # which is already done.
            vmin = min(bounds1.min(), -bounds2.min())
            vmax = max(bounds1.max(), -bounds2.max())
            cont1 = ax1.pcolormesh(lim_pbins, lim_rbins, bounds1, shading="auto", vmin=vmin,
              vmax=vmax, cmap="inferno")
            cont2 = ax2.pcolormesh(lim_pbins, lim_rbins, -bounds2, shading="auto", vmin=vmin,
              vmax=vmax, cmap="inferno")

            cbar = fig.colorbar(cont2, ax=[ax1, ax2], location="right")
            #cbar.set_label("Distance from 3DLIM origin", fontsize=12)

            ax1.set_title("Positive Bounds")
            ax2.set_title("Negative Bounds")
            ax1.set_xlabel("Z (m)", fontsize=12)
            ax1.set_ylabel("R (m)", fontsize=12)
            ax2.set_xlabel("Z (m)", fontsize=12)
            #fig.tight_layout()

            total_bounds = bounds1 - bounds2  # Minuc bc bounds2 is negative.
            cont3 = ax3.pcolormesh(lim_pbins, lim_rbins, total_bounds, shading="auto",
              vmin=total_bounds.min(), vmax=total_bounds.max(), cmap="inferno")
            ax3.set_title("Total Field Line Length")
            ax3.set_xlabel("Z (m)")
            cbar = fig.colorbar(cont3, ax=ax3)
            plt.show()

            # Load wall for our toroidal angle.
            wall = self.read_3d_wall(wall_path)
            wall_coords = wall[tor_angle]

            # Another plot with the 3DLIM points mapped to the equlibrium.
            fig, ax = plt.subplots(figsize=(5,8))
            ax.plot(wall_coords[0], wall_coords[1], color="k", zorder=3)
            #cont = ax.pcolormesh(gR, gZ, dists, shading="auto", cmap="coolwarm",
            #  vmin=-1, vmax=1)
            #cbar = fig.colorbar(cont)
            #ax.contour(gR, gZ, dists, levels=[0], colors="k")
            ax.contour(gR, gZ, psin, levels=np.linspace(0.95, 1.10, 8),
              colors="k", linewidths=1, zorder=30)
            ax.contour(gR, gZ, psin, levels=[1], colors="k", linewidths=3,
              zorder=40)
            #ax.scatter(lim_machRs.flatten(), lim_machZs.flatten(), s=3, color="k")
            masked_total = np.ma.masked_where(total_bounds<=0, total_bounds)
            cont = ax.pcolormesh(lim_machRs, lim_machZs, masked_total,
              cmap="inferno", shading="auto", vmin=total_bounds.min(),
              vmax=total_bounds.max(), zorder=50)
            ax.set_aspect("equal")
            ax.spines["top"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            fig.tight_layout()
            fig.show()
            plt.show()

        return debug_dict

    def plot_3dlim_on_wall(self, lim_path, tor_angle, rmrsep_origin, p_to_z,
      wall_path, plot_surfaces=False, gfile_pickle_path=None, show_plot=True):
        """
        NOT WORKING YET
        This takes a 3DLIM runs, with some input about where the simulation
        takes place, and then overlays it on a 2D cross section of the wall
        at the toroidal angle.

        Inputs
        lim_path (str): Path to the 3DLIM run that was setup for this particular
          toroidal angle.
        tor_angle (float): The toroidal angle that this 3DLIM run was setup for.
        r_origin (float): The R value, in machine coordinates, of the 3DLIM
          volume. I.e. if we are simulating a collector probe, then r_origin
          would be the R value of the tip of the probe.
        p_to_z ([float, float]): This is returned in bounds_file_from_mafot. It
          is the linear fit from Z = m * P + b, i.e. convert from P to Z. The
          first entry is m, and the second is b.
        wall_path (str): Path to the 3D wall file (probably called
          mafot_3D_wall.dat)
        plot_surfaces (bool): Whether or not to include flux surfaces on the
          plot.
        gfile_path (str): Path to a pickled object that has been created ahead
          of time from one of Shawn's scripts. NEEDS WORK on explaining how this
          is obtained.

        Outputs
        Returns a dictionary of the relevant data, as well as the created figure.
        """

        # Load in the 3DLIM results.
        try:
            nc = netCDF4.Dataset(lim_path)
        except FileNotFoundError:
            print("Error: Unable to find file: {}".format(lim_path))
            print("Exiting")

        ps = nc.variables["PS"][:].data
        xs = nc.variables["XS"][:].data
        ys = nc.variables["YS"][:].data
        ywids = nc.variables["YWIDS"][:].data
        xwids = nc.variables["XWIDS"][:].data
        pwids = nc.variables["PWIDS"][:].data
        ddlim3 = nc.variables["DDLIM3"][:].data
        vp = nc.variables["velplasma_4d_1"][:].data

        # Sum across all charge states.
        print("Summing across charge states...")
        ddlim3 = ddlim3.sum(axis=1)

        # Calculate centers of bins. Special treatment for the Y coordinate
        # since it needs to be mirrored and a zero added.
        rad_locs = xs - xwids / 2
        pol_locs = ps - pwids / 2
        tmp      = ys - ywids / 2
        par_locs = np.append(np.append(-tmp[::-1], 0), tmp)

        # Trim trailing/leading zeros off the arrays.
        y_keep_start = np.nonzero(par_locs)[0].min()
        y_keep_end   = np.nonzero(par_locs)[0].max() + 1
        x_keep_start = np.nonzero(rad_locs)[0].min()
        x_keep_end   = np.nonzero(rad_locs)[0].max() + 1
        p_keep_start = np.nonzero(pol_locs)[0].min()
        p_keep_end   = np.nonzero(pol_locs)[0].max() + 1
        rad_locs = rad_locs[x_keep_start:x_keep_end]
        pol_locs = pol_locs[p_keep_start:p_keep_end]
        par_locs = par_locs[y_keep_start:y_keep_end]
        ddlim3 = ddlim3[p_keep_start:p_keep_end, y_keep_start:y_keep_end,
          x_keep_start:x_keep_end]

        # Convert from 3DLIM R (which is along R-Rsep but not R-Rsep) to normal
        # R-Rsep.
        rmrseps = [rmrsep_origin - rb for rb in rad_locs]

        # Convert R to actual R coordinate, and P to Z coordinate.
        #rad_locs = r_origin - rad_locs
        zs = p_to_z[0] * pol_locs + p_to_z[1]  # From P to Z.

        # As done in the above function, we need to go from R-Rsep to machine
        # R. This is done with the gfile and interpolation function.
        with open(gfile_pickle_path, "rb") as f:
            gfile = pickle.load(f)
        gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
        psin = gfile["PSIRZ_NORM"]
        R_sep = gfile["RBBBS"]
        Z_sep = gfile["ZBBBS"]
        #debug_dict["gfile"] = gfile

        print("Calculating distances...")
        dists = np.zeros(gR.shape)
        for i in range(0, gR.shape[0]):
            for j in range(0, gR.shape[1]):
                d = np.sqrt(np.square(gR[i, j] - R_sep) + np.square(gZ[i, j] - Z_sep))
                if psin[i, j] < 1:
                    dists[i, j] = -d.min()
                else:
                    dists[i, j] = d.min()
        #debug_dict["gR"] = gR
        #debug_dict["gZ"] = gZ
        #debug_dict["dists"] = dists

        print("Interpolating 3DLIM coordinates onto grid...")
        print("  Note: Using every 10th EFIT data point")
        right = gR > gfile["RMAXIS"]
        int_R = gR[right].flatten()[::10]
        int_Z = gZ[right].flatten()[::10]
        int_dists = dists[right].flatten()[::10]
        coords = list(zip(int_dists, int_Z))
        f_R = RBFInterpolator(coords, int_R, smoothing=10)
        rmrseps_2d, zs_2d = np.meshgrid(rmrseps, zs)
        R_lim_2d = f_R(list(zip(rmrseps_2d.flatten(), zs_2d.flatten()))).reshape(zs_2d.shape)
        #debug_dict["f_R"] = f_R
        #debug_dict["zs_2d"] = zs_2d
        #debug_dict["R_lim_2d"] = R_lim_2d

        # Need to use the middle, i.e. Y = 0, since this is where the simulation
        # was setup for. Can't jump across toroidal angles in this case since
        # it would require different connection lengths.
        mid = int(len(par_locs) / 2)
        ddlim3_mid = ddlim3[:, mid-1:mid+2, :].sum(axis=1)  # Shape is (z_loc, rad_locs).

        # Now load in the 2D cross section of the wall at our angle.
        wall = self.read_3d_wall(wall_path)
        wall_coords = wall[tor_angle]

        masked_data = np.ma.masked_where(ddlim3_mid <= 0, ddlim3_mid)

        fig, ax = plt.subplots()
        ax.plot(wall_coords[0], wall_coords[1], color="k", zorder=3)
        ax.contourf(R_lim_2d, zs_2d, masked_data, cmap="inferno",
          norm=LogNorm(), zorder=4)
        ax.set_aspect("equal")

        if plot_surfaces:

            # Easiest way here would be to plot a normal contour plot at
            # specified psin steps, and to mask anything outside of the vessel.
            bbpath = mplpath.Path(list(zip(wall_coords[0], wall_coords[1])))
            bbpath_mask = ~bbpath.contains_points(np.array(list(zip(gR.flatten(),
              gZ.flatten()))))
            psin_masked = np.ma.masked_array(psin.flatten(),
              mask=bbpath_mask).reshape(psin.shape)
            ax.contour(gR, gZ, psin_masked, colors="k",
              levels=np.linspace(1.0, 1.5, 7), zorder=2)
            ax.contour(gR, gZ, psin_masked, colors="k",
              levels=[1.0], linewidths=2, zorder=1)

        fig.tight_layout()
        if show_plot:
            plt.show()

        # Return everything needed to make the plot elsewhere.
        return {"R_lim_2d":R_lim_2d, "zs_2d":zs_2d, "imp_density":masked_data,
          "wall":wall_coords, "rad_locs":rad_locs, "pol_locs":pol_locs,
          "par_locs":par_locs, "ddlim3":ddlim3, "gR":gR, "gZ":gZ, "psin_masked":psin_masked}
