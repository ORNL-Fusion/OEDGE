import os
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
      output_file="mafot_3d_wall.dat", axis_offset=180):
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
        axis_offset (float): The CAD file may not have the X axis at 0 degrees,
          and may instead have it at, say 180. The this will make sure that
          zero degrees is actually the machine zero degrees if it is not in the
          mesh file. This is normal angle units, where it increases in the
          CCW direction (machine angle is CW).

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

            # If there are no points within the polygon...
            if inside.sum() == 0:
                return rs, zs

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
        polys.append(Polygon([(0.90, 0.08),  (0.90, 0.18),  (1.15, 0.18),  (1.15, 0.08)]))

        # We want to take slices of the 3D volume to get a successive sequence of
        # 2D cross-sections, one at each degree. We do this by giving trimesh the
        # surface normal of the intersecting plane.
        sections = []
        degrees = np.arange(0.0, 360.0, angle_step)
        print("Generating cross sections...")
        for i in tqdm(range(0, len(degrees))):

            # Apply the offset to the axis here, as well as making sure the
            # angle (in normal, not machine) is between 0 and 360.
            deg = (degrees[i] - axis_offset) % 360
            #deg = degrees[i] % 360

            # Weird behavior that at zero degrees barely any points are
            # registered, but 360 is fine. I dunno man.
            if deg == 0:
                deg = 360

            # We add 90 since we want the surface normal at this angle (which is
            # 90 degrees away).
            # Something wrong here...
            #x_norm = np.cos(np.radians(deg + 180))
            #y_norm = np.sin(np.radians(deg + 180))
            x_norm = np.cos(np.radians(deg))
            y_norm = np.sin(np.radians(deg))

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
                keep = points[:,0] > 0
            else:
                keep = points[:,0] < 0
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

            # Re-sort the the small step in the inner wall.
            r, z = sort_polar((1.15, 0.00), r, z, polys[7])

            # Machine coordinates are clockwise so 360 - phi. axis_offset just
            # to get it to line up with actual machine coordinates (trial and
            # error).
            #phi_mach = (360 - deg + axis_offset) % 360
            phi_mach = 360 - degrees[i]
            print("{}-{} --> {} --> {} ({:.2f}, {:.2f})".format(degrees[i], axis_offset, deg, phi_mach, x_norm, y_norm))

            sections.append({"degree":degrees[i], "x_norm":x_norm, "y_norm":y_norm,
              "x_cs":x_cs, "y_cs":y_cs, "r":r, "z":z, "phi_mach":phi_mach})

        # Sort sections according to the angle (not machine angle) before
        # writing to MAFOT, which I have ASSUMED wants it in order of normal
        # angle, not machine.
        sections = sorted(sections, key=lambda i:i["degree"])

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

        # Optional flag to make a gif of the wall as you progress around the
        # vessel. Warning! This process uses a lot of memory because of all the
        # plots it makes.
        if make_gif:
            import imageio
            import os

            # Define a standard backend, especially if on Mac since Mac does a
            # horrible job clearing out memory after plots have been closed
            # (and we're making a lot of plots here!).
            # This just causes a crash on MacOS though, so for now will have to
            # deal with the memory leak.
            #matplotlib.use('TkAgg')

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

                #fname = "{}/{}.png".format(gif_path.split("/")[0], i)
                fname = "wall_gif/{}.png".format(i)
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

    def get_mach_coords(self, gR, gZ, psin, mid_r, mid_z, lim_pbins,
      along_coord=None, contour_idx=0):
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
        contour_idx (int):

        Output
        {lim_machR, lim_machZ} (dict): The machine coordinates at each Pbin
          location along the field line that intersects (mid_r, mid_z).
        """


        # Create 2D array of distance from input coordinate and find where
        # the distance is a minimum to find our psin we care about.
        dist = np.sqrt(np.square(mid_r-gR) + np.square(mid_z-gZ))
        psin_mid = psin[np.where(dist == dist.min())]
        if len(psin_mid) > 1:
            psin_mid = psin_mid[1]

        fig, ax = plt.subplots()

        """
        # Bit of hardcoding to ensure we pick the right region of flux tubes.
        # Have only barely tested for USN MiMES and DiMES so far.
        if type(along_coord) == type(None):
            print("Error: along_coord None not implemented yet.")

        # Essentially hard-coding the MiMES area here.
        elif along_coord == "R":
            keep = gR > r_limit
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
        """
        gR_keep = gR
        gZ_keep = gZ
        psin_keep = psin

        #print("gR")
        #print(gR_keep)
        #print("gZ")
        #print(gZ_keep)
        #print("psin")
        #print(psin_keep)
        #print("psin_mid")
        #print(psin_mid)

        # Able to hijack the contour function to grab our field line coordinates.
        cont2 = ax.contour(gR_keep, gZ_keep, psin_keep, levels=[psin_mid],
          colors="r")

        """
        # Difficult thing to programatically detect, but on the further out psin
        # values you can have multiple contours, which means grabbing
        # mid_line_rz might select not the one you want...
        if len(cont2.allsegs[0]) > 1:

            # This is not a fix, just a one-time fix that worked so I could avoid
            # this issue.
            mid_line_rz = cont2.allsegs[0][1]
        else:
            mid_line_rz = cont2.allsegs[0][0]
        """
        mid_line_rz = cont2.allsegs[0][contour_idx]

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
      gfile_pickle_path, show_plot=True, along_coord=None, wall_path=None,
      max_conn=None, r_limit=1.50):
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
        max_conn (float): Optional. Set a maximum on the connection length
          values.
        r_limit (float): Generally shouldn't need to change, but may need to
          if you're using a strange shaped equlibrium (maybe a negative
          triangularity shot for example).

        divimp_nc_path

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
                elif line[:10] == "# Z-grid: ":
                    numzs1 = int(line.split(":")[1])
                elif line[:8] == "# Rmin: ":
                    rmin1 = float(line.split(":")[1])
                elif line[:8] == "# Rmax: ":
                    rmax1 = float(line.split(":")[1])
                elif line[:8] == "# Zmin: ":
                    zmin1 = float(line.split(":")[1])
                elif line[:8] == "# Zmax: ":
                    zmax1 = float(line.split(":")[1])
                    break
        with open(mafot_file2) as f:
            for line in f:
                if line[:10] == "# R-grid: ":
                    numrs2 = int(line.split(":")[1])
                elif line[:10] == "# Z-grid: ":
                    numzs2 = int(line.split(":")[1])
                elif line[:8] == "# Rmin: ":
                    rmin2 = float(line.split(":")[1])
                elif line[:8] == "# Rmax: ":
                    rmax2 = float(line.split(":")[1])
                elif line[:8] == "# Zmin: ":
                    zmin2 = float(line.split(":")[1])
                elif line[:8] == "# Zmax: ":
                    zmax2 = float(line.split(":")[1])
                    break

        # If the origin is outside the MAFOT grid, it's not gonna work obvi.
        if r_origin < rmin1 or r_origin > rmax1:
            print("Error! Origin is outside MAFOT grid. Re-run MAFOT so that" + \
            " it covers the origin.")
            print("r_origin = {:.4f}   MAFOT range R=[{:.2f}, {:.2f}]".format(r_origin, rmin1, rmax1))
            sys.exit()
        if z_origin < zmin1 or z_origin > zmax1:
            print("Error! Origin is outside MAFOT grid. Re-run MAFOT so that" + \
            " it covers the origin.")
            print("z_origin = {:.4f}   MAFOT range Z=[{:.2f}, {:.2f}]".format(z_origin, zmin1, zmax1))
            sys.exit()

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

        # Due to rounding errors, just reassign the R's and Z's manually instead
        # of calling unique() on the dataframe.
        r = np.linspace(rmin1, rmax1, numrs1)
        z = np.linspace(zmin1, zmax1, numzs1)
        #r = df1["R (m)"].unique()
        #z = df1["Z (m)"].unique()
        #r = df1["R (m)"][:numrs1]
        #z = df1["Z (m)"].unique()
        R, Z = np.meshgrid(r, z)

        # Pull out each's connection length and pitch angles.
        l1 = df1["Lconn (km)"].values * 1000  # km to m
        l2 = df2["Lconn (km)"].values * 1000
        p1 = df1["pitch angle"].values
        p2 = df2["pitch angle"].values
        #L1 = l1.reshape(len(r), len(z))
        #L2 = l2.reshape(len(r), len(z))
        #pitch1 = p1.reshape(len(r), len(z))
        #pitch2 = p2.reshape(len(r), len(z))
        L1 = l1.reshape(len(z), len(r))
        L2 = l2.reshape(len(z), len(r))
        pitch1 = p1.reshape(len(z), len(r))
        pitch2 = p2.reshape(len(z), len(r))

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

        #if along_coord == "R":
        lim_machRs = np.zeros((len(lim_rbins), len(lim_pbins)))
        lim_machZs = np.zeros((len(lim_rbins), len(lim_pbins)))
        #elif along_coord == "Z":
        #    lim_machRs = np.zeros((len(lim_pbins), len(lim_rbins)))
        #    lim_machZs = np.zeros((len(lim_pbins), len(lim_rbins)))
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

            # Make sure that mach_coords is within the MAFOT grid. If the 3DLIM
            # bins just simply extend outside the grid, then a smaller bin
            # width or larger MAFOT grid needs to be made. This is on the user.
            # But sometimes a psin contour in get_mach_coords can be found at
            # multiple locations, so try again with the next contour and hope it
            # is within the MAFOT bounds.
            contour_idx = -1
            while True:
                contour_idx += 1
                mach_coords = self.get_mach_coords(gR, gZ, psin, machR, machZ,
                  lim_pbins, along_coord, contour_idx=contour_idx)

                if contour_idx > 999:
                    print("Error! Max contour_idx reached (>999). Unsure what the fix is here.")
                    sys.exit()

                # Failed. Point is not in MAFOT grid.
                if (mach_coords["lim_machR"] < r.min()).sum():
                    continue
                elif (mach_coords["lim_machR"] > r.max()).sum():
                    continue
                elif (mach_coords["lim_machZ"] < z.min()).sum():
                    continue
                elif (mach_coords["lim_machZ"] > z.max()).sum():
                    continue

                # Passed. Point is in MAFOT grid.
                else:
                    break

            #if along_coord == "R":
            lim_machRs[i] = mach_coords["lim_machR"]
            lim_machZs[i] = mach_coords["lim_machZ"]
            #elif along_coord == "Z":
            #    lim_machRs[:,i] = mach_coords["lim_machR"]
            #    lim_machZs[:,i] = mach_coords["lim_machZ"]
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
            print("Exiting. Check returned dictionary for info.")

            # Make a plot of the equilibrium with the selected points overlaid.
            fig, ax = plt.subplots()
            ax.contour(gR, gZ, psin, levels=np.linspace(0.95, 1.10, 8),
              colors="k", linewidths=1, zorder=30)
            ax.contour(gR, gZ, psin, levels=[1], colors="k", linewidths=3,
              zorder=40)
            ax.scatter(lim_machRs, lim_machZs, color="r", s=10)
            ax.set_aspect("equal")
            fig.tight_layout()
            plt.show()

            return debug_dict

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
        for ir in range(0, bounds1.shape[0]):
            for iz in range(0, bounds1.shape[1]):
                dist = np.sqrt(np.square(R - lim_machRs[ir][iz]) + np.square(Z - lim_machZs[ir][iz]))
                nearest_idx = np.where(dist == np.nanmin(dist))
                #print(len(nearest_idx))
                #if len(nearest_idx) > 1:
                #    i = nearest_idx[0][0]
                #    k = nearest_idx[0][1]
                #    nearest_idx = [[i], [k]]
                #bounds1[ir][iz] = L1[nearest_idx]
                #bounds2[ir][iz] = -L2[nearest_idx]
                try:
                    bounds1[ir][iz] = L1[nearest_idx]
                    bounds2[ir][iz] = -L2[nearest_idx]
                except:
                    print("{} {} {}".format(ir, iz, nearest_idx))
                    sys.exit()

        # Until I can figure out this bug in 3DLIM where ions that hit the end
        # seemingly restart at the top, set the last bound row to 0's.
        #print("Bug workaround: First bound row set to 0's")
        #bounds1[0] = np.zeros(len(bounds1[0]))
        #bounds2[0] = np.zeros(len(bounds2[0]))

        # If a maximum connection length value has been given, cap the bounds.
        if type(max_conn) != type(None):
            bounds1 = np.clip(bounds1, 0, max_conn)
            bounds2 = np.clip(bounds2, -max_conn, 0)

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

        # Print out the bin values to be copy/pasted to input file. Don't include
        # last value as it is assumed to be the CAP value, which is input
        # elsewhere in the input file.
        print("R bins for copy/paste to 3DLIM input file:")
        print("Number of bins: {}".format(len(lim_rbins)))
        for r in lim_rbins[:-1]:
            print("{:.4f}".format(r))
        print()
        print("P bins for copy/paste to 3DLIM input file:")
        print("Number of bins: {}".format(len(lim_pbins)))
        for p in lim_pbins:
            print("{:.4f}".format(p))
        print()
        nybins = 150
        print("Option for Y Bins")

        # Don't include zero as 3DLIM requires you don't include it.
        print("Number of bins: {}".format(nybins-1))
        lim_ybins = np.linspace(0, max(np.abs(bounds1).max(), np.abs(bounds2).max()), nybins)
        for y in lim_ybins[1:]:
            print("{:.2f}".format(y))
        print()
        print("Some additional comments")
        print("  - Change AW to something less than {:.3f}".format(lim_rbins[0]))
        print("  - Change A to something greater than {:.3f}".format(lim_rbins[-1]))
        #print("  - Maximum bounds value is {:.3f}".format(max(np.abs(bounds1).max(), np.abs(bounds2).max())))
        print("  - Bounds range for L21 and L19 is ({:.2f}, {:.2f})".format(bounds2.min(), bounds1.max()))
        #print("  - A bug workaround is applied that sets the first row to 0's")
        #print("These values are needed for option 3. Write them down!")
        #print("  - R-Rsep 3DLIM origin = {:.3f}".format(r_origin))
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
            vmin = min(bounds1.min(), (-bounds2).min())
            vmax = max(bounds1.max(), (-bounds2).max())
            #if along_coord == "R":
            cont1 = ax1.pcolormesh(lim_pbins, lim_rbins, bounds1, shading="auto", vmin=vmin,
              vmax=vmax, cmap="inferno")
            cont2 = ax2.pcolormesh(lim_pbins, lim_rbins, -bounds2, shading="auto", vmin=vmin,
              vmax=vmax, cmap="inferno")

            # Need to swap directions so thatt he plot makes sense from a DiMES
            # perspective.
            #elif along_coord == "Z":
            #    cont1 = ax1.pcolormesh(lim_rbins, lim_pbins, bounds1, shading="auto", vmin=vmin,
            #      vmax=vmax, cmap="inferno")
            #    cont2 = ax2.pcolormesh(lim_rbins, lim_pbins, -bounds2[::-1, :], shading="auto", vmin=vmin,
            #      vmax=vmax, cmap="inferno")

            cbar = fig.colorbar(cont2, ax=[ax1, ax2], location="right")
            #cbar.set_label("Distance from 3DLIM origin", fontsize=12)

            ax1.set_title("Positive Bounds")
            ax2.set_title("Negative Bounds")
            #if along_coord == "R":
            ax1.set_xlabel("3DLIM Binormal (m)", fontsize=12)
            ax1.set_ylabel("3DLIM Radial (m)", fontsize=12)
            ax2.set_xlabel("3DLIM Binormal (m)", fontsize=12)
            ax2.set_yticklabels([])
            #elif along_coord == "Z":
            #    ax1.set_xlabel("R (m)", fontsize=12)
            #    ax1.set_ylabel("Z (m)", fontsize=12)
            #    ax2.set_xlabel("R (m)", fontsize=12)
            #    ax2.set_yticklabels([])
            #fig.tight_layout()

            #if along_coord == "R":
            total_bounds = bounds1 - bounds2  # Minus bc bounds2 is negative.
            cont3 = ax3.pcolormesh(lim_pbins, lim_rbins, total_bounds, shading="auto",
              vmin=total_bounds.min(), vmax=total_bounds.max(), cmap="inferno")
            ax3.set_xlabel("3DLIM Binormal (m)")
            #elif along_coord == "Z":
            #    total_bounds = bounds1 - bounds2[::-1, :]  # Minus bc bounds2 is negative.
            #    cont3 = ax3.pcolormesh(lim_rbins, lim_pbins, total_bounds, shading="auto",
            #      vmin=total_bounds.min(), vmax=total_bounds.max(), cmap="inferno")
            #    ax3.set_xlabel("R (m)")
            ax3.set_title("Total Field Line Length")

            cbar = fig.colorbar(cont3, ax=ax3)
            plt.show()

            # Load wall for our toroidal angle.
            wall = self.read_3d_wall(wall_path)
            wall_coords = wall[tor_angle]
            debug_dict["wall_coords"] = wall_coords

            # Another plot with the 3DLIM points mapped to the equlibrium.
            fig, ax = plt.subplots(figsize=(5,8))
            ax.plot(wall_coords[0], wall_coords[1], color="k", zorder=3)
            #cont = ax.pcolormesh(gR, gZ, dists, shading="auto", cmap="coolwarm",
            #  vmin=-1, vmax=1)
            #cbar = fig.colorbar(cont)
            #ax.contour(gR, gZ, dists, levels=[0], colors="k")
            ax.contour(gR, gZ, psin, levels=np.linspace(0.95, 1.15, 8),
              colors="k", linewidths=1, zorder=30)
            ax.contour(gR, gZ, psin, levels=[1], colors="k", linewidths=3,
              zorder=40)
            #ax.scatter(lim_machRs.flatten(), lim_machZs.flatten(), s=3, color="k")
            #if along_coord == "R":
            masked_total = np.ma.masked_where(total_bounds<=0, total_bounds)
            #elif along_coord == "Z":
            #    masked_total = np.ma.masked_where(total_bounds.T<=0, total_bounds.T)

            # This generates a warning, too lazy to figure out how to fix.
            #cont = ax.pcolormesh(lim_machRs, lim_machZs, masked_total,
            #  cmap="inferno", shading="auto", vmin=total_bounds.min(),
            #  vmax=total_bounds.max(), zorder=50)

            #if along_coord == "R":
            cont = ax.contourf(lim_machRs, lim_machZs, masked_total,
              cmap="inferno", vmin=total_bounds.min(),
              vmax=total_bounds.max(), zorder=50)
            #elif along_coord == "Z":
            #    cont = ax.contourf(lim_machRs, lim_machZs, bounds1,
            #      cmap="inferno", vmin=bounds1.min(),
            #      vmax=bounds1.max(), zorder=50)
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


    def plot_3dlim_on_wall(self, lim_path, tor_angle, wall_path, r_origin,
      z_origin, along_coord=None, plot_surfaces=False, gfile_pickle_path=None,
      show_plot=True):
        """
        This takes a 3DLIM runs, with some input about where the simulation
        takes place, and then overlays it on a 2D cross section of the wall
        at the toroidal angle.

        Inputs
        lim_path (str): Path to the 3DLIM run that was setup for this particular
          toroidal angle.
        tor_angle (float): The toroidal angle that this 3DLIM run was setup for.
        wall_path (str): Path to the 3D wall file (probably called
          mafot_3D_wall.dat)
        r_origin (float): Machine R coordinate of the 3DLIM origin.
        z_origin (float): Machine Z coordinate of the 3DLIM origin.
        along_coord (str): One of "R", "Z" or None. Chooses whether the radial
          coordinate of 3DLIM is along the R direction (MiMES), Z direction
          (DiMES) or just the actual radial plasma coordinate (None, but this
          one isn't implemented yet).
        plot_surfaces (bool): Whether or not to include flux surfaces on the
          plot.
        gfile_path (str): Path to a pickled object that has been created ahead
          of time from one of Shawn's scripts. NEEDS WORK on explaining how this
          is obtained.
        show_plot (bool): Show the plot at the end. Useful if you just want the
          data returned for your own plots.

        Outputs
        Returns a dictionary of the relevant data, as well as the created figure.
        """

        debug_dict = {}

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
        lim_rbins = xs
        lim_pbins = ps
        rad_locs = xs - xwids / 2
        pol_locs = ps - pwids / 2
        tmp      = ys - ywids / 2
        par_locs = np.append(np.append(-tmp[::-1], 0), tmp)

        # Trim trailing/leading zeros off the arrays.
        lrb_keep_start = np.nonzero(lim_rbins)[0].min()
        lrb_keep_end   = np.nonzero(lim_rbins)[0].max() + 1
        lpb_keep_start = np.nonzero(lim_pbins)[0].min()
        lpb_keep_end   = np.nonzero(lim_pbins)[0].max() + 1
        y_keep_start = np.nonzero(par_locs)[0].min()
        y_keep_end   = np.nonzero(par_locs)[0].max() + 1
        x_keep_start = np.nonzero(rad_locs)[0].min()
        x_keep_end   = np.nonzero(rad_locs)[0].max() + 1
        p_keep_start = np.nonzero(pol_locs)[0].min()
        p_keep_end   = np.nonzero(pol_locs)[0].max() + 1
        lim_rbins = lim_rbins[lrb_keep_start:lrb_keep_end]
        lim_pbins = lim_pbins[lpb_keep_start:lpb_keep_end]
        rad_locs = rad_locs[x_keep_start:x_keep_end]
        pol_locs = pol_locs[p_keep_start:p_keep_end]
        par_locs = par_locs[y_keep_start:y_keep_end]
        ddlim3 = ddlim3[p_keep_start:p_keep_end, y_keep_start:y_keep_end,
          x_keep_start:x_keep_end]
        debug_dict["lim_rbins"] = lim_rbins
        debug_dict["lim_pbins"] = lim_pbins

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

        # Need to use the middle, i.e. Y = 0, since this is where the simulation
        # was setup for. Can't jump across toroidal angles in this case since
        # it would require different connection lengths.
        mid = int(len(par_locs) / 2)
        ddlim3_mid = ddlim3[:, mid-1:mid+2, :].sum(axis=1)  # Shape is (z_loc, rad_locs).

        # Now load in the 2D cross section of the wall at our angle.
        wall = self.read_3d_wall(wall_path)
        wall_coords = wall[tor_angle]

        masked_data = np.ma.masked_where(ddlim3_mid <= 0, ddlim3_mid).T
        debug_dict["masked_ddlim3"] = masked_data
        debug_dict["wall_coords"] = wall_coords

        fig, ax = plt.subplots()
        ax.plot(wall_coords[0], wall_coords[1], color="k", zorder=3)
        #ax.contourf(R_lim_2d, zs_2d, masked_data, cmap="inferno",
        #  norm=LogNorm(), zorder=4)
        ax.contourf(lim_machRs, lim_machZs, masked_data, cmap="inferno",
          norm=LogNorm(), zorder=4)
        ax.set_aspect("equal")

        # Easiest way here would be to plot a normal contour plot at
        # specified psin steps, and to mask anything outside of the vessel.
        bbpath = mplpath.Path(list(zip(wall_coords[0], wall_coords[1])))
        bbpath_mask = ~bbpath.contains_points(np.array(list(zip(gR.flatten(),
          gZ.flatten()))))
        psin_masked = np.ma.masked_array(psin.flatten(),
          mask=bbpath_mask).reshape(psin.shape)
        debug_dict["psin_masked"] = psin_masked

        if plot_surfaces:
            ax.contour(gR, gZ, psin_masked, colors="k",
              levels=np.linspace(1.0, 1.5, 7), zorder=2)
            ax.contour(gR, gZ, psin_masked, colors="k",
              levels=[1.0], linewidths=2, zorder=1)

        fig.tight_layout()
        if show_plot:
            plt.show()

        # Return everything needed to make the plot elsewhere.
        return debug_dict

    def divimp_prob_dist(self, divimp_nc_path, r_fluxtube, z_fluxtube,
      show_plot=True, reverse=False, smooth=False, smooth_window=None):
      """
      Using a DIVIMP run, get the probability distribution along a field line
      specified by r_fluxtube, z_fluxtube. The distribution is specified by
      the impurity density * the volume of the cell (area * 2piR). This is
      the total number of particles represented in this cell, from which the
      3DLIM probability distribution is based off of. ABSFAC as well.

      --- Inputs ---
      divimp_nc_path (str): Path to the DIVIMP netcdf file.
      r_fluxtube (float): The R coordinate of the flux tube you want impurity
        results from.
      z_fluxtube (float): Likewise.
      show_plot (bool): Whether to show a plot or not.
      reverse (bool): Whether to reverse the S values or not. Depending on your
        3DLIM setup, S = 0 may be referring to the wrong target, so set this to
        True to swap things around.
      smooth (bool): Return savgol filtered results (shown on plot too).
      smooth_window (int): The smoothing window for the savgol filter. Must be
        odd.
      """

      # Add to path the path to oedge_plots and import.
      sys.path.append(os.getcwd().split("lim")[0] + "oedge")
      import oedge_plots

      # First load the DIVIMP run.
      op = oedge_plots.OedgePlots(divimp_nc_path)
      cell = op.find_ring_knot(r_fluxtube, z_fluxtube, return_cell=True)
      ring, knot = op.find_ring_knot(r_fluxtube, z_fluxtube)
      print("Selecting ring {}".format(int(ring)))

      if show_plot:
          fig = op.plot_contour_polygon("Ring", vmin=ring-1, vmax=ring+1)
          ax = fig.axes[0]
          ax.scatter([r_fluxtube], [z_fluxtube], marker="o", s=50, color="k")
          ax.plot(cell[:,0], cell[:,1])
          fig.show()

      # Grab the impurity density and cell volume along our selected ring.
      s, imp = op.along_ring(int(ring), "DDLIMS", charge="all", plot_it=False)
      s, vol = op.along_ring(int(ring), "KVOLS", plot_it=False)
      tot_imp = imp * vol

      # If smoothing then do that with a savgol filter. Default window size
      # is to do 15 windows.
      if smooth:
          from scipy.signal import savgol_filter
          if type(smooth_window) == type(None):
              smooth_window = int(len(tot_imp)/15)
              if smooth_window % 2 == 0:
                  smooth_window += 1
          tot_imp_smooth = savgol_filter(tot_imp, smooth_window, 2)

      # Print out the density for copy/paste.
      print("Input for option Z02: {}".format(len(tot_imp)))
      if reverse:

          # Need to offset the S values by the min otherwise the S values will
          # not be the same.
          s   = s.min() + s.max() - s[::-1]
          tot_imp = tot_imp[::-1]
          if smooth:
              tot_imp_smooth = tot_imp_smooth[::-1]

      for i in range(0, len(tot_imp)):
          if smooth:
              print("{:<6.2f} {:.3e}".format(s[i], tot_imp_smooth[i]))
          else:
              print("{:<6.2f} {:.3e}".format(s[i], tot_imp[i]))

      # It's really on the user to make sure what they're inputting makes sense
      # for their simulation.
      print("\nImportant details!")
      print("  - Length of flux tube is {:.2f} m. ".format(s.max()) + \
        "The distance between\n      the injection region (Y0L - Y0S) " + \
        "should be\n      close to this number.")
      print("  - You need to make sure S = 0 corresponds to right side\n" + \
        "      in 3DLIM! Use reverse = True if the DIVIMP results \n" + \
        "      need to be flipped.")
      print("  - smooth_window = {}".format(int(smooth_window)))
      print("  - ABSFAC = {:.3e}".format(tot_imp.sum()))
      #print("  - sum(nz) = {:.3e}".format(imp.sum()))

      if show_plot:
          fig, ax = plt.subplots()
          ax.plot(s, tot_imp, color="k")
          if smooth:
              ax.plot(s, tot_imp_smooth, color="tab:red")
          ax.set_xlabel("Distance from target")
          ax.set_ylabel("Number of impurity ions")
          fig.tight_layout()
          plt.show()

      return {"s":s, "tot_imp":tot_imp, "imp":imp, "vol":vol, "ring":ring, "knot":knot}
