"""
Author  : Shawn Zamperini
Email   : zamp@utk.edu, zamperinis@fusion.gat.com
Updated : 5/25/22

This script started with some code written by Jake Nichols, but has moved on
to be it's own standalone script. It provides a framework to plot OEDGE output
data, mostly found in the netCDF file. The script oedge_plots_gui acts
as an interface to this class full of functions, though there is nothing
stopping someone from using this script on it's own.

List of functions:

__init__
  - Initiate an OedgePlots object by passing in the path to the NetCDF file.

add_dat_file
  - Pass in the path to the .dat file. The GUI should automatically find it if
      it's in the same folder as the NetCDF file.

read_data_2d
  - Read data from the NetCDF file into a form easily manipulated for plotting
      in plot_contour_polygon.

read_data_2d_kvhs
  - Same as above, but uses the .dat file to add on an additional, ad-hoc
      flow velocity (option T13 and associated parameters).

get_sep
  - Creates lines for plotting the separatrix in plot_contour_polygon.

calculate_forces
  - DIVIMP does not seem to output the forces in a meaningful way yet, so instead
      we just recalculate them here. Returns in a 2D form that can be used in
      plot_contour_polygon.

plot_contour_polygon
  - The primary function in this script. This is the one that creates the 2D plots
      of whatever you want from the NetCDF file, as well as the above mentioned
      plots that require a bit more data. There are many input parameters to change,
      so it's a very flexible function. Users of the GUI just get a plug-n-play
      usage of it.

plot_lp_input
  - Not implemented yet.

create_ts_from_omfit
  - This is a function that is a part of a larger workflow that uses the OMFIT
      module OMFITprofiles to pull Thomson scattering data for comparison to an
      OSM background solution. See the GitHub for details on this workflow.

compare_ts
  - This one uses data created from the above function and creates a PDF of
      comparisons.

check_ts
  - Helper function to make sure the TS data after being mapped to a common
      flux surface looks reasonable.

find_ring_knot
  - Helper function to find which ring, knot a given R, Z location is on.

fake_probe
  - See what a hypothetical probe (Langmuir, Mach) inserted into the plasma
      would look like.

along_ring
  - Choose a ring number and plot some data along it from target to target.
      Requires the variable name form the NetCDF file currently.

input_for_midpoint_shift
  - This enables the user to mess around with different offsets from the
      midpoint of the grid related to option G55. This is useful if you want
      to shift the midpoint from the default (halfway between the targets) to
      a more reaosnable location such as the OMP.

collector_probe
  - An experimental function that has yet to show it actually makes sense. The
      goal is to make a number of assumptions from a DIVIMP run that would
      apply to reproducing a double-sided collector probe deposition profile.
      This may not be appropriate, and regardless is subpar compared to a
      full-fledged 3DLIM simulation.
"""

import netCDF4
import warnings
import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt
import pandas            as pd
from collections         import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages


# Some issues can be fixed by enabling this backend, though I can't say what
# exactly is going on since this is some low-level stuff.
#mpl.use("TkAgg")

# A nice looking font.
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.family'] = 'DejaVu Sans'

# These are the "Tableau 20" colors as RGB. Note: I have learned this is
# unnecessary since matplotlib now includes this color palette.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

# A very small number.
epsilon = 1e-30

class OedgePlots:

    def __init__(self, netcdf_path):
        """
        Initialization for class object. Loads in the relevant netCDF variables
        needed for the plots.

        netcdf_path: Path location of .nc file.
        """

        # Load in the netCDF file.
        self.nc = netCDF4.Dataset(netcdf_path)
        self.ncpath = netcdf_path

        # Initialize variable to hold dat file.
        self.dat_file = None

        # Load in some netCDF data that is used a lot.
        self.rs     = self.nc['RS'][:]
        self.zs     = self.nc['ZS'][:]
        self.nrs    = self.nc['NRS'][:]
        self.nks    = self.nc['NKS'][:]
        self.area   = self.nc['KAREAS'][:]
        self.korpg  = self.nc['KORPG'][:]
        self.rvertp = self.nc['RVERTP'][:]
        self.zvertp = self.nc['ZVERTP'][:]
        self.rvesm  = self.nc['RVESM'][:]
        self.zvesm  = self.nc['ZVESM'][:]
        self.irsep  = self.nc['IRSEP'][:]
        self.qtim   = self.nc['QTIM'][:]
        self.kss    = self.nc['KSS'][:]
        self.kfizs  = self.nc['KFIZS'][:]
        self.ksmaxs = self.nc['KSMAXS'][:]
        self.irsep  = self.nc['IRSEP'][:]
        self.irwall = self.nc['IRWALL'][:]
        self.crmb   = self.nc['CRMB'][:]
        self.crmi   = self.nc['CRMI'][:]
        self.emi    = 1.602E-19
        self.cion   = self.nc['CION'][:]

        try:
            self.absfac = self.nc['ABSFAC'][:]
        except:
            print("Warning: Can't load ABSFAC (okay if DIVIMP was not run).")
            self.absfac = 1.0

        # Create a mesh of of the corners of the each cell/polygon in the grid.
        mesh = []
        num_cells = 0

        # Scan through the rings.
        for ir in range(self.nrs):

            # Scan through the knots.
            for ik in range(self.nks[ir]):

                # Get the cell index of this knot on this ring.
                index = self.korpg[ir,ik] - 1

                # Only if the area of this cell is not zero append the corners.
                if self.area[ir,ik] != 0.0:
                    vertices = list(zip(self.rvertp[index][0:4], self.zvertp[index][0:4]))
                    mesh.append(vertices)
                    num_cells = num_cells + 1

                    # Print out a warning is the cell center is not within the vertices.
                    cell = mpl.path.Path(list(vertices))
                    r = self.rs[ir, ik]
                    z = self.zs[ir, ik]
                    if not cell.contains_point([r, z]):
                        print("Error: Cell center not within vertices.")
                        print("  (ir, ik)    = ({}, {})".format(ir, ik))
                        print("  Vertices    = {}".format(vertices))
                        print("  Cell center = ({}, {})".format(r, z))

        # Save the results in the class.
        self.num_cells = num_cells
        self.mesh = mesh

        # Add the .dat file with initialization as well.
        #try:
        dat_path = netcdf_path.split(".nc")[0] + ".dat"
        self.add_dat_file(dat_path)
        #except:
        #    print("Error: Can't locate .dat file.")

    def __repr__(self):
        message = 'OedgePlots Object\n' + \
                  '  Title:       ' + self.nc['TITLE'][:].data.tostring().decode('utf-8') + '\n' + \
                  '  Date Run:    ' + self.nc['JOB'][:].data.tostring().strip()[:8].decode('utf-8') + '\n' + \
                  '  Grid:        ' + self.nc['EQUIL'][:].data.tostring().decode('utf-8') + '\n' + \
                  '  A Cool Dude: You\n' + \
                  '  Description: ' + self.nc['DESC'][:].data.tostring().decode('utf-8') + '\n'

        return message

    def add_dat_file(self, dat_path):
        """
        Quick function to read in the dat_file into the class. The file is not
        suitible for pandas or anything like that, so just read it in as a text
        with newlines at the end of each line.

        Input
        dat_path: Path location of .dat file.
        """
        self.datpath = dat_path
        with open(dat_path) as f:
            self.dat_file = f.read()

        try:
            self.cpu_time = int(self.dat_file.split('TOTAL CPU TIME USED' + \
                             '     (S)')[1].split('\n')[0])
        except:
            pass

    def read_data_2d(self, dataname, charge=None, scaling=1.0, fix_fill=False,
        no_core=False):
        """
        Reads in 2D data into a 1D array, in a form that is then passed easily
        to PolyCollection for plotting.

        Input
        dataname : The 2D data as named in the netCDF file.
        charge   : The charge state to be plotted, if applicable.
        scaling  : Scaling factor to apply to the data, if applicable. Secret
                     option is 'Ring' to just return the ring number at each cell.
        fix_fill : Swap out the default huge number that it seems is put in when
                     there is no data with just a zero. Better for plotting.
        no_core:   Exclude any core data.

        Output
        data : The data in a 2D format compatible with plot_contour_polygon.
        """

        # Get the 2D data from the netCDF file. Some special options included.
        # Normalized s coordinate of each cell.
        if dataname == 'Snorm':
            raw_data = self.nc['KSB'][:] / self.ksmaxs[:, None]
            raw_data = np.abs(raw_data + raw_data / 1.0 - 1.0)

        # Special flag for this array which has a time dimension.
        elif dataname == "blob_counts_time":
            raw_data = self.nc[dataname][:].sum(axis=2)

        # For ZEFFS, I'm pretty sure the first two entries of the first dimention
        # are just there to help calculate the third, which is Zeff.
        elif dataname == "ZEFFS":
            raw_data = self.nc["ZEFFS"][2]

        else:
            raw_data = self.nc[dataname][:]

        data = np.zeros(self.num_cells)

        # 'all' just sums up all charge states. So here instead of loop for speed.
        if charge == 'all':
            if dataname.lower() == "ddlims":

                # The first entry are "primary" neutrals, and the second are
                # "total" neutrals. The difference between these may just be
                # total includes things like self-sputtering perhaps. But
                # regardless, it sounds like we only want to include the second
                # entry onwards. E.g.
                # 0 - C0primary, 1 - C0total, 2 - C1+, ... 7 - C6+.
                #raw_data = raw_data[:-1].sum(axis=0)
                raw_data = raw_data[1:].sum(axis=0)
            else:
                raw_data = raw_data[1:].sum(axis=0)

        count = 0
        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):
                if self.area[ir, ik] != 0.0:

                    # Just make the data the ring number (+1 because they aren't
                    # named starting at zero, but 1).
                    if scaling == 'Ring':
                        data[count] = ir + 1
                    elif scaling == 'Knot':
                        data[count] = ik + 1
                    else:

                        # If charge is specifed, this will be the first
                        # dimension. See above, but for raw_data, charge=0 will
                        # be neutrals and charge=1 = Z1+ and onwards.
                        if charge in [None, 'all']:
                            data[count] = raw_data[ir][ik] * scaling
                        else:
                            #data[count] = raw_data[charge + 1][ir][ik] * scaling
                            data[count] = raw_data[charge][ir][ik] * scaling

                        # If we don't want the core data in the plot, then let's
                        # replace it with a small number.
                        if no_core:
                            if ir < self.irsep - 1:
                                data[count] = epsilon
                    count = count + 1

        # This will fix an issue where the netCDF will fill in values (that are
        # maybe outside the normal grid?) with 9.99E36. Swap these out with zeros
        # if desired.
        if fix_fill:
            fill_val = data.max()
            data[np.where(data == fill_val)] = None

        return data

    def read_data_2d_kvhs_t13(self, no_core=False):
        """
        Special function for plotting the flow velocity. This
        is because some DIVIMP options (T13, T31, T37?, T38?...) add
        additional flow values not reflected in KVHS. These additional values
        are in the .dat file, and thus it is required to run this function.

        Input
        no_core : Exclude the core data.

        Output
        data : The data in a 2D format compatible with plot_contour_polygon.
        """

        # Make sure .dat file has been loaded.
        if self.dat_file == None:
            print("Error: .dat file not loaded in. Run 'add_dat_file' first.")

        try:
            pol_opt = float(self.dat_file.split('POL DRIFT OPT')[1].split(':')[0])
            if pol_opt == 0.0:
                print("Error: Poloidal drift option T13 was not on for this run.")
                return -1

            # Get the relevant table for the extra drifts out of the .dat file.
            add_data = self.dat_file.split('TABLE OF DRIFT REGION BY RING - RINGS ' + \
                                            'WITHOUT FLOW ARE NOT LISTED\n')[1]. \
                                            split('DRIFT')[0].split('\n')

            # Split the data between the spaces, put into DataFrame.
            add_data = [line.split() for line in add_data]
            add_df = pd.DataFrame(add_data[1:-1], columns=['IR', 'Vdrift (m/s)',
                                  'S_START (m)', 'S_END (m)'], dtype=float). \
                                  set_index('IR')

            # Get the 2D data from the netCDF file.
            dataname = 'KVHS'
            scaling = 1.0 / self.qtim
            raw_data = self.nc[dataname][:]
            data = np.zeros(self.num_cells)

            # Convert the 2D data (ir, ik) into 1D for plotting in the PolyCollection
            # matplotlib function.
            count = 0
            for ir in range(self.nrs):
                for ik in range(self.nks[ir]):
                    if self.area[ir, ik] != 0.0:

                        # Put the data from this [ring, knot] into a 1D array.
                        data[count] = raw_data[ir][ik] * scaling

                        # If this ring has additional drifts to be added.
                        if ir in add_df.index:

                            # Then add the drift along the appropriate s (or knot) range.
                            if self.kss[ir][ik] > add_df['S_START (m)'].loc[ir] and \
                               self.kss[ir][ik] < add_df['S_END (m)'].loc[ir]:

                               data[count] = data[count] + add_df['Vdrift (m/s)'].loc[ir]

                        # If we don't want the core data in the plot, then let's
                        # replace it with just None.
                        if no_core:
                            if ir < self.irsep - 1:
                                #data[count] = None
                                data[count] = epsilon

                        count = count + 1

            return data

        except IndexError:
            print('Error: Was T13 on for this run?')
            return None

    def get_sep(self):
        """
        Return collection of lines to be plotted with LineCollection method
        of matplotlib for the separatrix.

        Output
        lines : List of coordinates to draw the separatrix in a format friendly
                 for LineCollection.
        """

        # Get (R, Z) coordinates of separatrix.
        rsep = self.rvertp[self.korpg[self.irsep-1,:self.nks[self.irsep-1]]][:,0]
        zsep = self.zvertp[self.korpg[self.irsep-1,:self.nks[self.irsep-1]]][:,0]
        nsep = len(rsep)
        lines=[]

        # Construct separatrix as a series of pairs of coordinates (i.e. a line
        # between the coordinates), to be plotted. Don't connect final point to first.
        for i in range(nsep-2):
            lines.append([(rsep[i], zsep[i]), (rsep[i+1], zsep[i+1])])

        return lines

    def calculate_forces(self, force, vz_mult = 0.0, charge=None, no_core=False):
        """
        Return 2D representations of the parallel forces (FF, FiG, FeG, FPG, FE)
        in the same format as read_data_2d for easy plotting. This data is not
        returned in the NetCDF file, so we calculate it here. Ideally we would
        use exactly what DIVIMP calculates, but that doesn't seem to be here yet.

        Input
        force   : One of 'FF', 'FiG', 'FeG', 'FPG', 'FE' or 'Fnet' (to sum
                  them all).
        charge  : Charge of impurity ion needed for some of the forces.
        vz_mult : Fraction of the background velocity used for vz (needed only
                  in FF).
        no_core : Exclude core data.

        Output
        force : The force data in a 2D format compatible with plot_contour_polygon.
        """

        # Temperature and density values for calculations.
        te = self.read_data_2d('KTEBS', no_core=no_core)
        ti = self.read_data_2d('KTIBS', no_core=no_core)
        ne = self.read_data_2d('KNBS',  no_core=no_core)
        col_log = 15
        qe      = 1.602E-19
        amu_kg  = 1.66E-27
        fact = np.power(self.qtim, 2) * qe / (self.crmi * amu_kg)

        # See if T13 was on. Important for FF calculations.
        try:
            pol_opt = float(self.dat_file.split('POL DRIFT OPT')[1].split(':')[0])
            if pol_opt == 0.0:
                t13 = False
            elif pol_opt == 1.0:
                print("Additional poloidal drift option T13 was ON.")
                t13 = True

        except:
            print("Warning: Unable to determine if T13 was on/off. Possible .dat \
                   file was not loaded?")
            t13 = False

        # Friction force calculations.
        if force.lower() in ['ff', 'fnet']:

            # Slowing down time.
            try:
                tau_s = 1.47E13 * self.crmi * ti * np.sqrt(ti / self.crmb) / \
                        ((1 + self.crmb / self.crmi) * ne * np.power(charge, 2) * col_log)

                #print("FF:  {:.2f} s".format(tau_s))
            except TypeError:
                print("Error: Charge of ion needed for FF calculation.")
                return None

            # TODO: Currently will assume impurity velocity is some fraction of the
            # background velocity (default zero), though this is obviously not true
            # and a better implementation would use the real values wherever they are.
            # This would require having the velocity distribution of the impurities.
            if t13:
                vi = self.read_data_2d_kvhs_t13()
            else:

                # Don't forget KVHS is scaled by 1/QTIM to get m/s.
                scaling = 1.0 / self.qtim
                vi = self.read_data_2d('KVHS', scaling=scaling, no_core=no_core)

            vz = vz_mult * vi

            # Calculate the force.
            ff = self.crmi * amu_kg * (vi - vz) / tau_s

            # If not 'fnet', then go ahead and return.
            if force.lower() == 'ff':
                return ff

        # TODO: Pressure gradient force calculations.
        if force.lower() in ['fpg', 'fnet']:

            # Parallel collisional diffusion time.
            #tau_par = 1.47E13 * self.crmi * ti * np.sqrt(ti / self.crmb) / \
            #          (ne * np.power(charge, 2) * col_log)

            # Need to probably figure this out.
            fpg = 0

            # If not 'fnet', then go ahead and return.
            if force.lower() == 'fpg':
                return fpg

        # Electron temperature gradient force calculations.
        if force.lower() in ['feg', 'fnet']:
            alpha = 0.71 * np.power(charge, 2)

            # The electron temperature gradient from the code. I don't understand
            # this scaling factor, but Jake uses it to get into presumably eV/m.
            kfegs = self.read_data_2d('KFEGS', scaling = qe / fact, no_core=no_core)

            # Calculate the force.
            feg = alpha * kfegs

            # If not 'fnet', then go ahead and return.
            if force.lower() == 'feg':
                return feg

        # Ion temperature gradient force calculations.
        if force.lower() in ['fig', 'fnet']:
            mu = self.crmi / (self.crmi + self.crmb)
            beta = 3 * (mu + 5 * np.sqrt(2) * np.power(charge, 2) * \
                   (1.1 * np.power(mu, 5/2) - 0.35 * np.power(mu, 3/2)) - 1) / \
                   (2.6 - 2 * mu + 5.4 * np.power(mu, 2))
            print("FiG: Beta = {:.2f}".format(beta))

            # Again, must scale the gradient with this scaling factor.
            kfigs = self.read_data_2d('KFIGS', scaling = qe / fact, no_core=no_core)

            # Calculate the force.
            fig = beta * kfigs

            # If not 'fnet', then go ahead and return.
            if force.lower() == 'fig':
                return fig

        # Electric field force calculations.
        if force.lower() in ['fe', 'fnet']:

            # This also gets scaled by a factor as well in Jake's code.
            e_pol = self.read_data_2d('E_POL', scaling = qe / fact, no_core=no_core)
            fe = charge * qe * e_pol

            # If not 'fnet', then go ahead and return.
            if force.lower() == 'fe':
                return fe

        # Net force calculation.
        if force.lower() == 'fnet':
            return ff + fpg + feg + fig + fe

    def plot_contour_polygon(self, dataname, charge=None, scaling=1.0,
                             normtype='linear', cmap='plasma', xlim=[0.9, 2.4],
                             ylim = [-1.4, 1.4], plot_sep=True, levels=None,
                             cbar_label=None, fontsize=16, lut=21,
                             smooth_cmap=False, vmin=None, vmax=None,
                             show_cp=None, ptip=None, show_mr=False,
                             fix_fill=False, own_data=None, no_core=False,
                             vz_mult=0.0, wall_data=None, show_grid=True,
                             figsize=(5, 6)):

        """
        Create a standalone figure using the PolyCollection object of matplotlib.
        This is the main function here for plotting 2D data, and you should be
        able to control most whatever you want through the options since this
        was made to be a versatile function.

        Input
        dataname:    The netCDF variable name to plot. Some special datanames
                       will perform extra data handling: KVHSimp, ...
        charge:      The charge state to be plotted, if applicable.
        scaling:     Scaling factor to apply to the data, if applicable.
        normtype:    One of 'linear', 'log', ... of how to normalize the data on
                       the plot.
        cmap:        The colormap to apply to the plot. Uses standard matplotlib
                       names.
        xlim:        X range of axes.
        ylim:        Y range of axes.
        plot_sep:    Include separatrix on plot or not.
        levels:      Number of levels for colorbar (needs work).
        cbar_label:  Label for the colorbar.
        fontsize:    Size of font for labels.
        lut:         Number of chunks to break the colormap into.
        smooth_cmap: Choose whether to break colormap up into chunks or not.
        vmin/vmax:   Option to choose own vmin, vmax for the colorbars. Useful
                       when plots need tweaking.
        show_cp:     Choose if collector probes are to be shown. Can have
                       multiple probes as a list. 1 = MiMES, 2 = Top, 3 = DiMES.
        ptip:        The tip (either R or Z, depends on which probe you picked)
                       of the probe. Can pass as list where each entry corresponds
                       to the one in show_cp.
        show_mr:     Option to show the metal rings or not from MRC-I.
        fix_fill:    The defaults for some masked_arrays is to fill in blank
                       values with 1e36. I think zero is more appropriate for
                       plotting correctly.
        own_data:    Bring your own data to work day. This could be used if you
                      want to plot the ratio of two data point or something. So if
                      you wanted to ratio of two datasets returned by read_data_2d,
                      then you could say own_data = data1 / data2.
        no_core:     Don't include data in the core region.
        vz_mult:     An optional input for if you're plotting the friction force.
                       See calculate_forces for more info.
        wall_data:   Can supply your own R, Z wall coordinates, say if you need
                       modify it to remove part of a limiter or something. Format
                       is a tuple (Rcoords, Zcoords) in meters.
        show_grid:   Whether to show the grid or not.

        Output
        fig : The plotted Figure object.
        """

        # Make sure show_cp and ptip is in list form if not.
        if type(show_cp) is not list:
            show_cp = [show_cp]
        if type(ptip) is not list:
            ptip = [ptip]

        # Option to provide own dataset from read_data_2d, say like a ratio
        # between two of them or something that you calculated ahead of itme.
        if own_data is not None:
            data = own_data

        # Get the data through normal means.
        else:

            # Read in the data into a form for PolyCollection. Account for
            # special options.
            # Flow velocity with additional velocity specified by T13.
            if dataname == 'KVHSimp':
                data = self.read_data_2d_kvhs_t13(no_core=no_core)

            # Special option to plot the ring numbers.
            elif dataname == 'Ring':
                data = self.read_data_2d('KTEBS', scaling='Ring', no_core=no_core)
            elif dataname == 'Knot':
                data = self.read_data_2d('KTEBS', scaling='Knot', no_core=no_core)

            # Divide the background velocity by the sounds speed to get the Mach number.
            elif dataname == 'KVHSimp - Mach':
                te   = self.read_data_2d('KTEBS', no_core=no_core)
                ti   = self.read_data_2d('KTIBS', no_core=no_core)
                kvhs = self.read_data_2d_kvhs_t13()
                mi   = self.crmb * 931.494*10**6 / (3e8)**2  # amu --> eV s2 / m2
                cs   = np.sqrt((te + ti) / mi)
                data = kvhs / cs  # i.e. the Mach number.

            elif dataname == 'KVHS - Mach':
                te   = self.read_data_2d('KTEBS', no_core=no_core)
                ti   = self.read_data_2d('KTIBS', no_core=no_core)
                kvhs = self.read_data_2d('KVHS',  no_core=no_core, charge=charge,
                                         scaling=scaling, fix_fill=fix_fill)
                mi   = self.crmb * 931.494*10**6 / (3e8)**2  # amu --> eV s2 / m2
                cs   = np.sqrt((te + ti) / mi)
                print("CRMB = {} amu".format(self.crmb))
                data = kvhs / cs  # i.e. the Mach number.

            # Special function for plotting the forces on impuirties.
            elif dataname.lower() in ['ff', 'fig', 'feg', 'fpg', 'fe', 'fnet']:
                data = self.calculate_forces(dataname, charge=charge,
                                             no_core=no_core, vz_mult=vz_mult)

            # Everything else in the netCDF file.
            else:
                data = self.read_data_2d(dataname, charge, scaling, fix_fill, no_core)

        # Remove any cells that have nan values.
        not_nan_idx = np.where(~np.isnan(data))[0]
        mesh = np.array(self.mesh)[not_nan_idx, :, :]
        data = data[not_nan_idx]

        # Check if the array is all zeros. If so don't plot anything because
        # it will cause errors below.
        if not data.any():
            print("Error: Data is all zeros. Plot not made.")
            return None

        # Create a good sized figure with correct proportions.
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot(111)

        if normtype == 'linear':
            if vmin == None:
                vmin = data.min()
            if vmax == None:
                vmax = data.max()
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        elif normtype == 'log':

            # Set zero data to the minimum non-zero value. Strictly speaking
            # this isn't correct, but as far as the colormap is concerned
            # it doesn't matter.
            data[data == 0.0] = data[data != 0.0].min()
            #data = np.ma.masked_where(data==0.0, data)
            if vmin == None:
                vmin = data.min()
            if vmax == None:
                vmax = data.max()
            norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

        elif normtype == 'symlin':
            if vmin == None:
                vmin = -np.abs(data).max()
            if vmax == None:
                vmax = -vmin
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            cmap = 'coolwarm'

        elif normtype == 'symlog':
            data[data == 0.0] = data[data != 0.0].min()
            if vmin == None:
                vmin = -np.abs(data[~np.isnan(data)]).max()
            if vmax == None:
                vmax = -vmin
            norm = mpl.colors.SymLogNorm(linthresh=0.01 * vmax, vmin=vmin, vmax=vmax, base=10)
            cmap = 'coolwarm'

        # Shamelessly copied from StackOverflow. The end of nipy_spectral is
        # grey, which makes it look like there's a hole in the largest data. Fix
        # this by just grabbing a subset of the colormap from 0.0-0.95, leaving
        # out the last portion that's grey.
        if cmap == 'nipy_spectral':
            cmap_obj = plt.get_cmap(cmap)
            new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
              'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_obj.name, a=0, b=0.95),
              cmap_obj(np.linspace(0, 0.95, 500)), N=lut)
            cmap = new_cmap

        # Choose whether to discretize the colormap or not first.
        if smooth_cmap:
            scalar_map = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        else:
            scalar_map = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.get_cmap(cmap, lut=lut))

        # Create PolyCollection object.
        if show_grid:
            edgecolors="none"
        else:
            edgecolors="face"
        coll = mpl.collections.PolyCollection(mesh, array=data,
                                              cmap=scalar_map.cmap,
                                              norm=scalar_map.norm,
                                              edgecolors=edgecolors)

        # Add the PolyCollection to the Axes object.
        ax.add_collection(coll)

        # Plot the wall.
        if wall_data is None:

            # Drop all the (0, 0)'s and append the first point on the end so it
            # doesn't leave a gap in the plot.
            keep_idx = np.where(np.logical_and(self.rvesm[0] != 0, self.zvesm[0] != 0))
            rvesm = np.append(self.rvesm[0][keep_idx], self.rvesm[0][keep_idx][0])
            zvesm = np.append(self.zvesm[0][keep_idx], self.zvesm[0][keep_idx][0])
            #ax.plot(self.rvesm[0][keep_idx], self.zvesm[0][keep_idx], color='k', linewidth=1)
        else:
            rvesm = wall_data[0]
            zvesm = wall_data[1]
        ax.plot(rvesm, zvesm, color='k', linewidth=1)

        # Get the separatrix coordinates as a collection of lines and plot.
        if plot_sep:
            sep = self.get_sep()
            sc = mpl.collections.LineCollection(sep, color='k')
            ax.add_collection(sc)

        # Use correct amount of levels for colorbar, if specified.
        if levels is not None:
            cbar = fig.colorbar(coll, ax=ax, boundaries=levels, ticks=levels, extend='both', shrink=0.8)
        else:
            cbar = fig.colorbar(coll, ax=ax, extend='both', shrink=0.8)

        # Add colorbar label.
        if cbar_label is None:
            cbar.ax.set_ylabel(dataname, fontsize=fontsize)
        else:
            cbar.ax.set_ylabel(cbar_label, fontsize=fontsize)

        # A bit larger colorbar tick label font size.
        cbar.ax.tick_params(labelsize=14)

        # Option to add collector probes to plots.
        if show_cp:
            if ptip == None:
                print("Error: Location of tip of probe not given.")

        # MiMES probe.
        cp_width  = 0.03
        facecolor = 'grey'
        edgecolor = 'black'
        if 1 in show_cp:
            ptip_tmp = ptip[show_cp.index(1)]
            lower_xy = (ptip_tmp, -0.185 - cp_width/2.0)
            width    = 2.38 - ptip_tmp
            height   = cp_width
            rect = mpl.patches.Rectangle(lower_xy, width=width, height=height,
                                         facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(rect)

        # Top crown probe. Point is to simulate a DiMES probe if it were USN.
        if 2 in show_cp:
            ptip_tmp = ptip[show_cp.index(2)]
            lower_xy = (1.485 - cp_width/2.0, ptip_tmp)
            width    = cp_width
            height   = ptip_tmp - 1.172
            rect = mpl.patches.Rectangle(lower_xy, width=width, height=-height,
                                         facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(rect)

        # DiMES probe.
        if 3 in show_cp:
            ptip_tmp = ptip[show_cp.index(3)]
            lower_xy = (1.485 - cp_width/2.0, -1.25)
            width    = cp_width
            height   = ptip_tmp - (-1.25)
            rect = mpl.patches.Rectangle(lower_xy, width=width, height=height,
                                         facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(rect)

        # Option to show metal rings. Floor (1.32-1.37) Shelf (1.404, 1.454).
        if show_mr:
            tile_height = 0.01
            facecolor   = 'red'
            edgecolor   = 'black'
            floor_xy    = (1.32, -1.363 - tile_height)
            floor_width = 1.37 - 1.32
            shelf_xy    = (1.404, -1.250 - tile_height)
            shelf_width = 1.454 - 1.404
            floor_rect  = mpl.patches.Rectangle(floor_xy, width = floor_width,
                                               height=tile_height,
                                               facecolor=facecolor,
                                               edgecolor=edgecolor)
            shelf_rect  = mpl.patches.Rectangle(shelf_xy, width = shelf_width,
                                               height=tile_height,
                                               facecolor=facecolor,
                                               edgecolor=edgecolor)
            ax.add_patch(floor_rect)
            ax.add_patch(shelf_rect)

        # Organize plot, add labels.
        ax.axis('equal')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('R (m)', fontsize=fontsize)
        ax.set_ylabel('Z (m)', fontsize=fontsize)
        fig.tight_layout()
        fig.show()

        return fig

    def plot_lp_input(self):
        """
        May need to find the data in the .lim or .dat file.
        """
        pass

    def create_ts_from_omfit(self, shot, t_window, omfit_path=None,
        output_path=None, filt_abv_avg=None, smooth=False, smooth_window=11,
        filter=False, plot_it=True, core_shift=0.0, div_shift=0.0,
        div_sas_shift=0.0):
        """
        This creates a similar file as create_ts, the file that is used to create
        PDF's of plots comparing the OEDGE solution to the TS data, except instead
        of directly pulling the TS data, it relies on a file generated from
        OMFITprofiles. This is because OMFITprofiles does a number of slicing
        and filtering techniques that are superior to anything I can write. See the
        file create_omfit_excel.py or the GitHub for instructions on generating
        this input file (most likely it's called omfit_excel_file.csv).

        Filtering order is median --> average --> savgol.

        Input
        shot         : Right now can just handle a single shot, so input it here.
                        Though nothing is stopping you from running this multiple
                        times and just copy/pasting the two Excel files together,
                        one for each shot.
        t_window     : Restrict the data to just data between these two times.
                        Input as list or tuple, e.g. [2000, 5000]
        omfit_path   : Path to the Excel file created from create_omfit_excel.
        output_path  : Name of the Excel file that is output, and ready to be
                        used in compare_ts. Leaving as None will generate a
                        generic filename for you.
        filt_abv_avg : Inevitably, some of the data from OMFIT will still have
                        spikes in it from ELMs. This parameter will filter out
                        anything above the average Te/ne value for that chord.
                        Ex. filt_abv_avg = 1.5 will exclude anything above 1.5 *
                        average value, or in other words anything that is 50%
                        above the average.
        smooth        : Apply a savgol_filter to each chord to help smooth the
                         data out some.
        smooth_window : The window length for the savgol filter.
        filter        : Apply a median filter to the data using the same window
                         as smooth_window. Applied after savgol filtering.
        plot_it       : (bool) Option to plot the data to see if smoothing and
                         all went alright.
        core_shift    : (float) Shift in psin to apply to the core data.
        div_shift     : (float) Shift in psin to apply to the divertor (R-1)
                         data.
        div_sas_shift : (float) Shift in psin to apply to the SAS divertor (R+1)
                         data.

        Output
        output_path : The filename the data was saved in.
        """

        # Helper function to load data into a dataframe.
        def load_omfit_file(path, shot):
            if path.split(".")[1] == "csv":
                print("Loading csv...")
                df = pd.read_csv(path)
            elif path.split(".")[1] == "xlsx":
                print("Loading Excel...")
                df = pd.read_excel(path)
            else:
                print("Error: Cannot recognize filetype of omfit_file.")

            # Add column with shot number.
            return df

        # Pass in a list of omfit files, each being paired up with a corresponding
        # list of shots. Combine into one large DataFrame.
        if type(shot) == int:
            shot = [shot]
        if len(omfit_path) != len(shot):
            print("Error: omfit_path and shot must be the same length.")
        omfit_df = pd.DataFrame()
        for i in range(0, len(omfit_path)):
            omfdf = load_omfit_file(omfit_path[i], shot[i])
            omfdf["shot"] = np.full(len(omfdf), shot[i])
            omfit_df = pd.concat([omfit_df, omfdf])

        # Need to reset index since later operations assume unique indices.
        omfit_df.reset_index(inplace=True)

        # Apply shifts to the core and divertor data if specified.
        core_r1_shift = omfit_df[omfit_df["subsystem"] == "core_r+1"]["psin"] + core_shift
        core_r0_shift = omfit_df[omfit_df["subsystem"] == "core_r+0"]["psin"] + core_shift
        div_shift = omfit_df[omfit_df["subsystem"] == "divertor_r-1"]["psin"] + div_shift
        div_sas_shift = omfit_df[omfit_df["subsystem"] == "divertor_r+1"]["psin"] + div_sas_shift
        omfit_df.update(core_r1_shift)
        omfit_df.update(core_r0_shift)
        omfit_df.update(div_shift)
        omfit_df.update(div_sas_shift)

        # Remove data outside of the time window.
        keep = np.logical_and(omfit_df["time"].astype(float)>=t_window[0],
          omfit_df["time"].astype(float)<=t_window[1])
        omfit_df = omfit_df[keep]

        # Save raw data for plotting comparison later.
        raw_omfit_df = omfit_df.copy()

        # Flag to print warning at the end of everything if window size was changed.
        window_warning = False

        # Apply a median filter to each chord. This does not remove any data
        # points, as the median filter returns an array of equal size.
        if filter:
            print("Applying median filter for each channel...")
            from scipy.signal import medfilt
            pd.options.mode.chained_assignment = None

            # For each system...
            for system in ["core_r+1", "divertor_r-1", "divertor_r+1"]:
                try:
                    sys_df = omfit_df[omfit_df["subsystem"]==system]
                except:
                    print("No data for {}".format(system))
                    continue

                # For each channel...
                for chan in sys_df["channel"].unique():

                    # Get the data points from the chord, replace with a median
                    # filtered signal.
                    chan_df = sys_df[sys_df["channel"]==chan]
                    if smooth_window > len(chan_df["te"]):
                        smaller_window = int(len(chan_df["te"]) / 1.5)
                        if smaller_window % 2 == 0:
                            smaller_window -= 1
                            if smaller_window <= 0:
                                smaller_window = 1
                        smooth_te = medfilt(chan_df["te"], smaller_window)
                        smooth_ne = medfilt(chan_df["ne"], smaller_window)
                        print("  {}: channel {:d}: Decreasing window size to {:d}".format(system, int(chan), smaller_window))
                        window_warning = True
                    else:
                        smooth_te = medfilt(chan_df["te"], smooth_window)
                        smooth_ne = medfilt(chan_df["ne"], smooth_window)

                    chan_df["te"] = smooth_te
                    chan_df["ne"] = smooth_ne
                    omfit_df.update(chan_df)

        # Perform filter above average methodology. This process will iterate
        # and continue to remove data points until everything is within the
        # specified multiplier of the average.
        if filt_abv_avg != None:
            print("Filtering points above average theshold...")
            print("  Before len(omfit_df) = {:d}".format(len(omfit_df)))

            # Will rebuild the filtered data into a new DataFrame, and then
            # overwrite omfit_df with it.
            new_omfit_df = pd.DataFrame()

            # For each system...
            for system in ["core_r+1", "divertor_r-1", "divertor_r+1"]:
                try:
                    sys_df = omfit_df[omfit_df["subsystem"]==system]
                except:
                    print("No data for {}".format(system))
                    continue

                # For each channel...
                for chan in sys_df["channel"].unique():

                    # Get the data points from the chord.
                    chan_df = sys_df[sys_df["channel"]==chan]
                    points_dropped_te = 0
                    while True:

                        # Find points above the threshold.
                        drop_te = (chan_df["te"] > filt_abv_avg * chan_df["te"].mean())
                        drop_te = drop_te[drop_te == True].index

                        # Nothing left above the threshold.
                        if len(drop_te) == 0:
                            break

                        # Remove from this channel. As this is a view into the
                        # original omfit_df, this should be enough.
                        chan_df["te"].drop(drop_te, axis=0, inplace=True)
                        points_dropped_te += len(drop_te)

                    # Same for ne.
                    points_dropped_ne = 0
                    while True:
                        drop_ne = (chan_df["ne"] > filt_abv_avg * chan_df["ne"].mean())
                        drop_ne = drop_ne[drop_ne == True].index
                        if len(drop_ne) == 0:
                            break
                        chan_df["ne"].drop(drop_ne, axis=0, inplace=True)
                        points_dropped_ne += len(drop_ne)

                    #omfit_df.update(chan_df)
                    print("  {} {:d}: ".format(system, int(chan)) + \
                      "Dropped {:} ne ".format(points_dropped_ne) + \
                      "and {:} Te points ".format(points_dropped_te) + \
                      "({:} total)".format(points_dropped_te + points_dropped_ne))

                    new_omfit_df = pd.concat([new_omfit_df, chan_df])

            # Reassemble omfit_df of the filtered channels.
            omfit_df = new_omfit_df
            print("  After len(omfit_df) = {:}".format(len(omfit_df)))


        # Smooth with just a simple third order savgol filter.
        if smooth:
            print("Smoothing data for each channel...")
            from scipy.signal import savgol_filter
            pd.options.mode.chained_assignment = None

            # For each system...
            for system in ["core_r+1", "divertor_r-1", "divertor_r+1"]:
                try:
                    sys_df = omfit_df[omfit_df["subsystem"]==system]
                except:
                    print("No data for {}".format(system))
                    continue

                # For each channel...
                for chan in sys_df["channel"].unique():

                    # Get the data points from the chord, replace with a savgol
                    # filtered signal.
                    chan_df = sys_df[sys_df["channel"]==chan]
                    try:
                        smooth_te = savgol_filter(chan_df["te"], smooth_window, 2, mode="constant", cval=chan_df["te"].median())
                        smooth_ne = savgol_filter(chan_df["ne"], smooth_window, 2, mode="constant", cval=chan_df["ne"].median())
                    except Exception as e:

                        # Use smaller window size the size of the whole channel.
                        # Make sure it's an odd number too.
                        smaller_window = int(len(chan_df["te"]) / 1.5)
                        if smaller_window % 2 == 0:
                            smaller_window -= 1
                        if smaller_window < 2:
                            print("  {}: channel {:d}: Only 1 data point! Skip smoothing.".format(system, int(chan)))
                            continue
                        print("  {}: channel {:d}: Decreasing window size to {:d}".format(system, int(chan), smaller_window))
                        smooth_te = savgol_filter(chan_df["te"], smaller_window, 2, mode="constant", cval=chan_df["te"].median())
                        smooth_ne = savgol_filter(chan_df["ne"], smaller_window, 2, mode="constant", cval=chan_df["ne"].median())
                        window_warning = True
                    chan_df["te"] = smooth_te
                    chan_df["ne"] = smooth_ne
                    omfit_df.update(chan_df)

        # First, let's find the rings, knots and S values for each measurement location.
        print("Comparing to grid...")
        rings = []; knots = []; s = []; systems = []; channels = []
        for row in range(0, len(omfit_df)):

            # Coordinates for this TS measurement.
            ts_psin = omfit_df.iloc[row]['psin']
            ts_r = omfit_df.iloc[row]['r']
            ts_z = omfit_df.iloc[row]['z']

            # Find what ring, and then knot, this psin corresponds to. Why
            # wouldn't I just use find_ring_knot? Because R, Z is constant so
            # we'd get no variation of the ring it falls on! Using psin first
            # enables spreading the wealth some.
            ring = self.find_ring_from_psin(ts_psin)
            knot = self.find_knot_on_ring(ring, ts_r, ts_z)
            #ring, knot = self.find_ring_knot(ts_r, ts_z)

            # Get the S coordinate.
            s_tmp = self.kss[ring][knot]

            # Append to lists. The +1 is because the rings are named starting at
            # 1, not 0, so when writing them down we want to have the right naming
            # convention. Annoying, I know, but proper.
            rings.append(ring + 1)
            knots.append(knot + 1)
            s.append(s_tmp)

            # Change the system name to our naming convention.
            sys = omfit_df.iloc[row]['subsystem']
            if sys in ['core_r+1', 'core_r+0']:
                systems.append('core')
            elif sys == 'divertor_r-1':
                systems.append('divertor')
            elif sys == 'tangential_r+0':
                systems.append('tangential')
            elif sys == "divertor_r+1":
                systems.append('divertor_sas')

        # Print warning for if window size changed.
        if window_warning:
            message = "Warning! Window size was changed for at least one " + \
            "chord. This probably means it missed out on at least one of " + \
            "the filtering methods and won't look as good. Options are " + \
            "to either just remove this chord completely from your output " + \
            "file, or better yet to pass in multiple omfit_path's to " + \
            "combine the data from repeat shots so the window shortening " + \
            "procedure won't kick in (because you'll have more data)."
            print(message)

        # Organize into the output dataframe with the expected format.
        te   = omfit_df['te'].values
        ne   = omfit_df['ne'].values
        psin = omfit_df['psin'].values
        r    = omfit_df['r'].values
        z    = omfit_df['z'].values
        time = omfit_df['time'].values
        channels = omfit_df['channel'].values
        #shot = np.full(len(time), shot)
        shot = omfit_df["shot"].values
        data = np.vstack((s, te, ne, rings, knots, systems, psin, channels, shot, r, z, time))

        self.output_df = pd.DataFrame(data.T, columns=('S (m)', 'Te (eV)', 'ne (m-3)',
                                 'Ring', 'Cell', 'System', 'Psin', 'Channel', 'Shot', 'R (m)',
                                 'Z (m)', 'Time'))

        # Make sure the data types of each column are correct.
        self.output_df['S (m)']    = self.output_df['S (m)'].astype(float)
        self.output_df['Te (eV)']  = self.output_df['Te (eV)'].astype(float)
        self.output_df['ne (m-3)'] = self.output_df['ne (m-3)'].astype(float)
        self.output_df['Ring']     = self.output_df['Ring'].astype(int)
        self.output_df['Cell']     = self.output_df['Cell'].astype(int)
        self.output_df['System']   = self.output_df['System'].astype(np.str)
        self.output_df['Psin']     = self.output_df['Psin'].astype(float)
        self.output_df['Channel']  = self.output_df['Channel'].astype(float).astype(int)
        self.output_df['Shot']     = self.output_df['Shot'].astype(float).astype(int)
        self.output_df['R (m)']    = self.output_df['R (m)'].astype(float)
        self.output_df['Z (m)']    = self.output_df['Z (m)'].astype(float)
        self.output_df['Time']     = self.output_df['Time'].astype(float)

        # Plot the resulting data to see how things turned out.
        if plot_it:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,6), sharex=True)
            core    = self.output_df[self.output_df["System"]=="core"]
            core_x  = core["Psin"]
            core_te = core["Te (eV)"]
            core_ne = core["ne (m-3)"]
            divertor    = self.output_df[self.output_df["System"]=="divertor"]
            divertor_x  = divertor["Psin"]
            divertor_te = divertor["Te (eV)"]
            divertor_ne = divertor["ne (m-3)"]
            divertor_sas    = self.output_df[self.output_df["System"]=="divertor_sas"]
            divertor_sas_x  = divertor_sas["Psin"]
            divertor_sas_te = divertor_sas["Te (eV)"]
            divertor_sas_ne = divertor_sas["ne (m-3)"]
            core_raw = raw_omfit_df[raw_omfit_df["subsystem"]=="core_r+1"]
            core_x_raw = core_raw["psin"]
            core_te_raw = core_raw["te"]
            core_ne_raw = core_raw["ne"]
            divertor_raw = raw_omfit_df[raw_omfit_df["subsystem"]=="divertor_r-1"]
            divertor_x_raw = divertor_raw["psin"]
            divertor_te_raw = divertor_raw["te"]
            divertor_ne_raw = divertor_raw["ne"]
            divertor_sas_raw = raw_omfit_df[raw_omfit_df["subsystem"]=="divertor_r+1"]
            divertor_sas_x_raw = divertor_sas_raw["psin"]
            divertor_sas_te_raw = divertor_sas_raw["te"]
            divertor_sas_ne_raw = divertor_sas_raw["ne"]
            marker = "."
            alpha = 0.6
            ax1.plot(core_x, core_te, marker, label="Core", color="tab:purple", alpha=alpha)
            ax1.plot(divertor_x, divertor_te, marker, label="DTS", color="tab:red", alpha=alpha)
            ax1.plot(divertor_sas_x, divertor_sas_te, marker, label="SAS", color="tab:green", alpha=alpha)
            ax3.plot(core_x, core_ne, marker, label="Core", color="tab:purple", alpha=alpha)
            ax3.plot(divertor_x, divertor_ne, marker, label="DTS", color="tab:red", alpha=alpha)
            ax3.plot(divertor_sas_x, divertor_sas_ne, marker, label="SAS", color="tab:green", alpha=alpha)
            ax2.plot(core_x_raw, core_te_raw, marker, color="tab:purple", alpha=alpha)
            ax2.plot(divertor_x_raw, divertor_te_raw, marker, color="tab:red", alpha=alpha)
            ax2.plot(divertor_sas_x_raw, divertor_sas_te_raw, marker, color="tab:green", alpha=alpha)
            ax4.plot(core_x_raw, core_ne_raw, marker, color="tab:purple", alpha=alpha)
            ax4.plot(divertor_x_raw, divertor_ne_raw, marker, color="tab:red", alpha=alpha)
            ax4.plot(divertor_sas_x_raw, divertor_sas_ne_raw, marker, color="tab:green", alpha=alpha)
            ax2.set_ylim(ax1.get_ylim())
            ax4.set_ylim(ax3.get_ylim())
            ax1.set_ylabel("Te (eV)", fontsize=14)
            ax3.set_ylabel("ne (m-3)", fontsize=14)
            ax3.set_xlabel("Psin", fontsize=14)
            ax4.set_xlabel("Psin", fontsize=14)
            ax1.set_title("Filtered/Smoothed", fontsize=14)
            ax2.set_title("Raw Data", fontsize=14)
            ax1.legend()
            fig.tight_layout()
            fig.show()
            print("plot shown")

        # Output to a new Excel file.
        print("Saving to Excel...")
        if output_path == None:
            output_path = 'ts' + str(shot[0]) + '_from_omfit.xlsx'
        self.output_df.to_excel(output_path)
        return output_path

    def add_rcp_data(self, ts_path, rcp_path, rcp_sheet, output_path=None):
        """
        This function will take the Excel file generated from
        create_ts_from_omfit and add on plunging probe data, if available. This
        will interpolate between RCP measurement points to ensure there is a
        comparison point on every cell the RCP path goes through. This function
        can be run more than once, say to add MiMES RCP data and then X-point
        RCP data.

        Input
        ts_path (str): The path to the file output from create_ts_from_omfit.
        rcp_path (str): Path to the RCP data Excel file. Columns MUST be the following,
          R (m), Z (m), ne (m-3), Te (eV), Mach, vi (m/s), qpar (W/m2). If any data is
          not available just include the column name anyways and leave the rows
          beneath it blank.

        Output
        output_path (str): Path to the save the file to. If None will generate
          a generic filename for you.

        """

        rcp_cols = ["R (m)", "Z (m)", "ne (m-3)", "Te (eV)", "Mach",
          "vi (m/s)", "qpar (W/m2)"]

        print("Loading files...")
        ts_df  = pd.read_excel(ts_path)
        rcp_df = pd.read_excel(rcp_path, sheet_name=rcp_sheet, usecols=rcp_cols)

        # Find which ring, knot each measurement location is on on the grid.
        print("Comparing to grid...")
        rings = []; knots = []; s = []; systems = []
        for row in range(0, len(rcp_df)):

            # Using the R, Z locations, find the ring, knot for each location.
            rcp_r = rcp_df.iloc[row]['R (m)']
            rcp_z = rcp_df.iloc[row]['Z (m)']
            ring, knot = self.find_ring_knot(rcp_r, rcp_z)

            # Get the S coordinate.
            s_tmp = self.kss[ring][knot]

            # Append to lists. The +1 is because the rings are named starting at
            # 1, not 0, so when writing them down we want to have the right naming
            # convention. Annoying, I know, but proper.
            rings.append(ring + 1)
            knots.append(knot + 1)
            s.append(s_tmp)

        # Hijack the systems column for indicating it's RCP.
        systems = np.full(len(rcp_df), "rcp")

    def print_core_ts(self, ncpath):
        """
        This function will use the data gathered from OMFIT to print out
        data for copy/paste for the core TS input in the DIVIMP input file.
        The file it needs is the NetCDF file you can export at the end of the
        OMFITprofiles GUI via "Export OMFITprofiles to NetCDF". pass that in
        here.

        --- Inputs ---
        omfit_path (str): Path to the csv file created from create_omfit_excel.
        nbins (int): The number of psin bins to bin the data into. Averages
          of each bin are returned.
        """

        nc = netCDF4.Dataset(ncpath)
        psin = nc.variables["psi_n"][:]
        te = nc.variables["T_e"][:].mean(axis=0)
        ne = nc.variables["n_e"][:].mean(axis=0)
        for p, t, n in zip(psin, te, ne):
            print("{:5.3f} {:7.2f} {:7.2f} {:8.2e} 0".format(p, t, t, n))



    def compare_ts(self, ts_filename, rings, show_legend='all', nrows=3,
        ncols=2, bin_width=1.0, output_file='my_ts_comparison.pdf',
        rad_bin_width=0.0025, core_sweep_bug_time=None, filter_zeros=True,
        dashed_rings=True):
        """
        Function to create plots to compare the OEDGE results to Thomson
        scattering. A PDF file is output with all the graphs for each ring. This
        function only really makes sense with L-mode shots since it has no ELM
        filtering capablities.

        Input
        ts_filename   : The Excel file created from "create_ts_from_omfit" above. Has all the
                         TS data mapped to S.
        rings         : OEDGE rings to compare the TS against. INDEXED STARTING AT 1, NOT 0.
        show_legend   : One of 'all' or 'short'. 'all' will show which shot is
                         which color, and 'short' will only show OEDGE and the
                         average TS values, binned accordinagly with std. dev.
        nrows         : Number of rows of plots for the PDF file.
        ncols         : Number of columns of plots for the PDF file.
        bin_width     : Width of each bin, in meters, to bin the TS data along S
                         into. This makes getting averages and a standard deviation
                         possible, making comparing TS to OEDGE much more meaningful.
        output_file   : What the save the PDF as.
        rad_bin_width : Bin width for the R-Rsep OMP plots to bin the TS data into.
        core_sweep_bug_time : Handles a bug in the ThomsonClass script. See
                         explanation in the comment about 15 lines down from here.
        filter_zeros  : Get rid of some zeros that aren't real data.
        dashed_rings  : Put dashed lines on the plots for every tenth ring.
        """

        # Warning that pops up but is unecessary.
        plt.rcParams.update({'figure.max_open_warning': 0})

        # Load in the DataFrame from Excel file. Could maybe save time by seeing
        # if user already has it loaded in from create_ts above, but meh.
        self.s_df_all = pd.read_excel(ts_filename)
        df = self.s_df_all

        # There is a slight bug in the ThomsonClass script that is fortunately
        # easy to work around (love u pandas). To map to a common reference time,
        # it maps each location to the X-point. This allows 2D graphs with the
        # divertor system, but can make it seem like the core is sweeping a
        # range, when it really isn't. Use this variable to indicate the start
        # of the strike point sweep so we can exclude that core data.
        if core_sweep_bug_time is not None:
            drop_idx = np.logical_and(df['System'] == 'core', df['Time'] > core_sweep_bug_time).values
            df = df.iloc[~drop_idx]

        # Get rid of some random zeros that muck things up.
        if filter_zeros:
            df = df[df['Te (eV)']  > 0.0]
            df = df[df['ne (m-3)'] > 0.0]

        # Determine if USN/LSN by looking at the Z coordinate of the X-point.
        if self.nc["ZXP"][:].data > 0.0:
            config = "usn"
        else:
            config = "lsn"
        print("Plasma in {} configuration...".format(config.upper()))

        # Variables to keep track of which Axes we are on.
        ngraphs = nrows * ncols
        graph_count = 0

        # Create pdf filename. Open up file to save Figures to.
        pdf_filename = output_file
        with PdfPages(pdf_filename) as pdf:

            # Get the R, Z of the OMP.
            z0 = self.nc['Z0'][:].data
            r0 = self.nc['R0'][:].data

            # Get Rsep at the OMP.
            sep_ring = self.nc['IRSEP'][:] - 1  # Subtract 1 for indexing means.
            omp_side = self.rs[sep_ring] > r0
            zs_omp   = self.zs[sep_ring][omp_side]
            rs_omp   = self.rs[sep_ring][omp_side]
            dist     = np.abs(zs_omp - z0)
            #dist         = np.abs(self.zs[sep_ring] - z0)
            rsepomp_cell = np.where(dist == dist.min())[0]
            rsepomp      = rs_omp[rsepomp_cell]
            zsepomp      = zs_omp[rsepomp_cell]
            #rsepomp      = self.rs[sep_ring, rsepomp_cell]

            # Get the Thomson R's and Z's. Only upstream and rings outside IRSEP.
            # Downstream comparisons will have to be with Langmuir probes and
            # the parallel to s graphs that utilize the divertor TS.
            tmp_df  = df[np.logical_and(df['Ring'] >= self.irsep, df['Ring'] < 100)] # Just to stop far rings from goofing up the plot.
            ts_r    = tmp_df[tmp_df['System'] == 'core']['R (m)'].values
            ts_z    = tmp_df[tmp_df['System'] == 'core']['Z (m)'].values
            ts_ring = tmp_df[tmp_df['System'] == 'core']['Ring'].values  # Subtract 1 for indexing means.
            ts_te   = tmp_df[tmp_df['System'] == 'core']['Te (eV)'].values
            ts_ne   = tmp_df[tmp_df['System'] == 'core']['ne (m-3)'].values
            ts_cell = tmp_df[tmp_df['System'] == 'core']['Cell'].values
            ts_romp = np.zeros(len(ts_r))

            # First time will be the ref_time since that's HOW I MADE IT WORK DAMMIT.
            ref_time = tmp_df[tmp_df['System'] == 'core']['Time'].values[0]
            self.tmp_df = tmp_df
            """
            # Find which knot corresponds to the OMP. This loop may need work...
            for i in range(0, len(ts_r)):

                if ts_ring[i] >= self.irwall:
                    continue

                # Only want the OMP side, so take cells on the right half of R0.
                omp_side = self.rs[ts_ring[i]] > r0

                #dist = np.abs(self.zs[ts_ring[i]][omp_side] - z0)
                #omp_cell = np.where(dist == dist.min())[0]
                #romp = self.rs[ts_ring[i]][omp_side][omp_cell]
                #print("Ring {}. romp = {}".format(ts_ring[i], romp))

                # Find R-Rsep OMP.
                #ts_romp[i] = romp - rsepomp

                # We are making a bit of an approximation here. If we were to
                # take the closest knot to the OMP for each ring, then it's possible
                # that each ring's corresponding knot could either be slightly above
                # or below the OMP, making R-Rsep OMP negative (which is obivously
                # not true). So instead we assume the OMP cell of each ring is the
                # same cell as the separatrix one calculated above. So in effect
                # the "radial" profile at the OMP we end up with may actually
                # be at a slight angle upwards or downwards, depending on the grid,
                # but we assume it to be negligible (probably a safe assumption).
                try:
                    romp = self.rs[ts_ring[i]][omp_side][rsepomp_cell]
                except:
                    print('R0           = {:.2f}'.format(r0))
                    print('ts_ring[{}]  = {}'.format(i, ts_ring[i]))
                    print('omp_side     = {}'.format(omp_side))
                    print('rsepomp_cell = {}'.format(rsepomp_cell))
                zomp = self.zs[ts_ring[i]][omp_side][rsepomp_cell]
                #old_settings = np.seterr(all='ignore')
                romp_approx = np.sqrt(np.power(romp-rsepomp, 2) + np.power(zomp-zsepomp, 2))
                #np.seterr(**old_settings)
                ts_romp[i] = romp_approx

            # Now need to get the OEDGE values at the OMP.
            more_rings = np.arange(sep_ring, self.irwall)
            oedge_romps = np.array([])
            oedge_teomp = np.array([])
            oedge_neomp = np.array([])

            # Get the R location where core Thomson is taken.
            tmp_df2 = tmp_df[tmp_df['System'] == 'core']
            ts_r_ref = tmp_df2[tmp_df2['Time'] == ref_time]['R (m)'].values[0]

            # Array to hold Romp values for plotting dashed lines for every
            # 10th ring (i.e. 20, 30, 40, ...).
            oedge_ring_dashed = np.array([])
            oedge_ring_romps  = np.array([])

            for i in range(0, len(more_rings)):

                # Get the ring we are interested in.
                ring = more_rings[i]
                #print(ring)

                try:


                    # Get the RminRsep OMP values for these values mapped to the OMP.
                    omp_side = self.rs[ring] > r0
                    #dist = np.abs(self.zs[ring][omp_side] - z0)
                    #omp_cell = np.where(dist == dist.min())[0]
                    #romp = self.rs[ring][omp_side][omp_cell]
                    #print("Ring {} omp_cell = {}, {}, {}".format(ring, omp_cell, romp, romp-rsepomp))
                    #if romp-rsepomp < 1e-4:
                    #    omp_cell = omp_cell - 1
                    #    romp = self.rs[ring][omp_side][omp_cell]
                    #    print("  Updated to {}, {}, {}".format(omp_cell, romp, romp-rsepomp))
                    #oedge_romps = np.append(oedge_romps, romp - rsepomp)


                    # Overwrite with the approximation already explained above.
                    romp = self.rs[ring][omp_side][rsepomp_cell]
                    zomp = self.zs[ring][omp_side][rsepomp_cell]
                    #old_settings = np.seterr(all='ignore')
                    romp_approx = np.sqrt(np.power(romp-rsepomp, 2) + np.power(zomp-zsepomp, 2))
                    #np.seterr(**old_settings)
                    oedge_romps = np.append(oedge_romps, romp_approx)

                    if ring % 10 == 0:
                        oedge_ring_dashed = np.append(oedge_ring_dashed, int(ring))
                        #oedge_ring_romps  = np.append(oedge_ring_romps, romp - rsepomp)
                        oedge_ring_romps  = np.append(oedge_ring_romps, romp_approx)
                        #print("Added to dashed rings: {} at {}".format(ring, romp-rsepomp))

                    # Get the Te, ne values along the R value where TS is taken. Need
                    # to find which cell on this ring this is at. Just to be safe only
                    # grab the portion above the OMP.
                    above_omp = self.zs[ring] > z0
                    dist = np.abs(self.rs[ring][above_omp] - ts_r_ref)
                    close_cell = np.where(dist == dist.min())[0]
                    oedge_teomp = np.append(oedge_teomp, self.nc['KTEBS'][:][ring][above_omp][close_cell])
                    oedge_neomp = np.append(oedge_neomp, self.nc['KNBS'][:][ring][above_omp][close_cell])

                except ValueError:
                    print("Maybe ring {} isn't to the right of the OMP.".format(ring))
                    pass

            # Sort the OEDGE data for plotting.
            idx = np.argsort(oedge_romps)
            oedge_romps = oedge_romps[idx]
            oedge_teomp = oedge_teomp[idx]
            oedge_neomp = oedge_neomp[idx]

            # Make figures, put them into the PDF.
            fig, axs = plt.subplots(2, 1, sharex=True)
            axs = axs.flatten()

            for data in ['Te', 'ne']:
                if data == 'Te':
                    y = ts_te
                    oedge_y = oedge_teomp
                    graph_num = 0
                else:
                    y = ts_ne
                    oedge_y = oedge_neomp
                    graph_num = 1

                # Bin the Thomson data some to get averages and error bars.
                bins = np.arange(ts_romp.min(), ts_romp.max()+rad_bin_width, rad_bin_width)
                digi = np.digitize(ts_romp, bins)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    bin_means = [y[digi == i].mean() for i in range(1, len(bins))]
                    bin_stds = [y[digi == i].std() for i in range(1, len(bins))]
                    bin_centers = [ts_romp[digi == i].mean() for i in range(1, len(bins))]
                #bin_centers = bins[:-1] + rad_bin_width / 2.0

                # Save max values for setting the plot limits.
                if data == 'Te':
                    te_bin_max = max(bin_means)
                else:
                    ne_bin_max = max(bin_means)

                line1, = axs[graph_num].plot(ts_romp, y, '.', color='lightgrey', alpha=0.35, zorder=-32)
                line2  = axs[graph_num].errorbar(bin_centers, bin_means, fmt='k.',
                                                 yerr=bin_stds,
                                                 capthick=1.5, capsize=3)
                line3, = axs[graph_num].plot(oedge_romps, oedge_y, 'k')
                axs[graph_num].legend((line2, line3), ('Thomson', 'OEDGE'),
                                        loc='upper center', ncol=3,
                                        bbox_to_anchor=(0.5, 1.1),
                                        fancybox=True, shadow=True)

                # Vertical dashed lines for every 10th ring.
                if dashed_rings:
                    axs[graph_num].vlines(oedge_ring_romps, 0, 1e30, colors='k', linestyles='dashed')

                    # Annotate with the rings number of this line.
                    ann_locs_te = list(zip(oedge_ring_romps-0.001, np.full(len(oedge_ring_romps), 150*0.77)))
                    ann_locs_ne = list(zip(oedge_ring_romps-0.001, np.full(len(oedge_ring_romps), 3e19*0.77)))
                    for i in range(0, len(oedge_ring_dashed)):
                        if graph_num == 0:
                            axs[graph_num].annotate(str(int(oedge_ring_dashed[i])), ann_locs_te[i], bbox=dict(facecolor='white', edgecolor='black'))
                        else:
                            axs[graph_num].annotate(str(int(oedge_ring_dashed[i])), ann_locs_ne[i], bbox=dict(facecolor='white', edgecolor='black'))

            axs[0].set_xlim([0, ts_romp.max() * 1.1])
            axs[0].set_ylim([0, 100])
            #axs[0].set_ylim([0, te_bin_max * 1.3])
            axs[1].set_ylim([0, 1.5e19])
            #axs[1].set_ylim([0, ne_bin_max * 1.3])
            axs[1].set_xlabel('R-Rsep OMP (m)')
            axs[0].set_ylabel('Te (eV)')
            axs[1].set_ylabel('ne (m-3)')
            fig.tight_layout()
            # papertype depreciated.
            #pdf.savefig(fig, papertype='letter')
            pdf.savefig(fig)
            """
            # Make parallel to s graphs. One loop for Te, one loop for ne.
            for oedge_data in ['Te', 'ne']:
                graph_count = 0

                # Go through one ring at a time.
                for ring in rings:

                    # Index only the rows of the DataFrame with that ring.
                    ring_df = df[df['Ring'] == ring]  # +1 since "ring" is zero-indexed here.

                    # Use modulo to determine if this ring starts a new fig/pdf page.
                    if graph_count % ngraphs == 0:
                        fig, axs = plt.subplots(nrows, ncols, sharex=False, figsize=(10,8))
                        axs = axs.flatten()

                    # Find what color this shot is from the above shots_enum. Elaborate...
                    color = 6
                    shots = df['Shot'].unique()
                    for shot in shots:

                        # Get only the portion matching this shot.
                        smaller_ring_df = ring_df[ring_df['Shot'] == shot]
                        s_core = smaller_ring_df[smaller_ring_df['System'] == 'core']['S (m)']
                        s_div  = smaller_ring_df[smaller_ring_df['System'] == 'divertor']['S (m)']

                        # Choose correct data for Y axis.
                        if oedge_data == 'Te':
                            y = smaller_ring_df['Te (eV)']
                            y_core = smaller_ring_df[smaller_ring_df['System'] == 'core']['Te (eV)']
                            y_div  = smaller_ring_df[smaller_ring_df['System'] == 'divertor']['Te (eV)']
                        if oedge_data == 'ne':
                            y = smaller_ring_df['ne (m-3)']
                            y_core = smaller_ring_df[smaller_ring_df['System'] == 'core']['ne (m-3)']
                            y_div  = smaller_ring_df[smaller_ring_df['System'] == 'divertor']['ne (m-3)']

                        # Plot it with a specific color. 2 because those colors are better.
                        # zorder is needed because of a matplotlib bug where the
                        # errorbars render under these data points, which we don't
                        # want.
                        #ls = axs[graph_count].plot(smaller_ring_df['S (m)'], y, '.',
                        #                           color=tableau20[color], label=shot,
                        #                           alpha=1.0, zorder=-32)
                        ls = axs[graph_count].plot(s_core, y_core, '.',
                                                   color=tableau20[color], label=shot,
                                                   alpha=0.5, zorder=-32)
                        ls = axs[graph_count].plot(s_div, y_div, 'o',
                                                   color=tableau20[color], label=shot,
                                                   alpha=0.5, zorder=-32, markerfacecolor='none', markeredgecolor=tableau20[color])
                        color += 1

                        # Prevent going past the size of tableau20 array.
                        if color > 23:
                            color -= 24

                    # Now plot the average value + std. dev. of the TS measurements
                    # with requested binning performed.
                    if oedge_data == 'Te':
                        y = ring_df['Te (eV)']
                    if oedge_data == 'ne':
                        y = ring_df['ne (m-3)']

                    # Bin the s values, get the mean and std of the corresponding y values.
                    if len(ring_df.index) == 0:
                        continue
                    s = ring_df['S (m)']
                    #print('Ring: ' + str(ring))
                    #print('smin: ' + str(s.min()))
                    #print('smax: ' + str(s.max()))

                    # If there is essentially just one bin.
                    if np.abs(s.min() - s.max()) < bin_width:
                        bin_means   = y.mean()
                        bin_stds    = y.std()
                        bin_centers = s.mean()

                    # Else do the binning procedure.
                    else:
                        bins = np.arange(s.min(), s.max()+bin_width, bin_width)
                        digi = np.digitize(s, bins)
                        bin_means = [y[digi == i].mean() for i in range(1, len(bins))]
                        bin_stds = [y[digi == i].std() for i in range(1, len(bins))]
                        #bin_centers = bins[:-1] + bin_width / 2.0
                        bin_centers = [s[digi == i].mean() for i in range(1, len(bins))]

                    line1 = axs[graph_count].errorbar(bin_centers, bin_means,
                                                      yerr=bin_stds, fmt='k.',
                                                      capthick=1.5, capsize=3)

                    # Find the rings from the OEDGE simulation now. Ring-1 for indexing.
                    oedge_s  = self.nc['KSS'][:][ring-1].data

                    # Te plots.
                    if oedge_data == 'Te':
                        oedge_dat = self.nc['KTEBS'][:][ring-1].data
                        axs[graph_count].set_ylabel('Te (eV)')
                        #axs[graph_count].set_ylim([0, 200])

                    # ne plots.
                    if oedge_data == 'ne':
                        oedge_dat = self.nc['KNBS'][:][ring-1].data
                        axs[graph_count].set_ylabel('ne (m-3)')
                        #axs[graph_count].set_ylim([0, 1.5e20])

                    # Drop some extra zeros that sneak in.
                    keep_idx = np.where(oedge_s != 0.0)[0]
                    oedge_s  = oedge_s[keep_idx]
                    oedge_dat = oedge_dat[keep_idx]

                    # Restrict to the outer half since that's where TS data is.
                    #half_idx = np.where(oedge_s > oedge_s.max()/2.0)[0][0]
                    #oedge_s  = oedge_s[half_idx:]
                    #oedge_dat = oedge_dat[half_idx:]

                    # Adjust limits.
                    try:
                        max_bin = max(bin_means)
                    except:
                        max_bin = 0
                    max_val = max([max_bin, max(oedge_dat)]) * 1.5
                    #print("{}, {}: {}, {} ({})".format(oedge_data, ring, max_bin, max(oedge_dat), max_val))
                    axs[graph_count].set_ylim([0, max_val])

                    # Add OEDGE data to plot. Note comma after 'line2', because plot
                    # returns two items, but errorbar doesn't.
                    line2, = axs[graph_count].plot(oedge_s, oedge_dat, 'k-')

                    # Extra plot commands.
                    try:
                        psin = self.nc.variables['PSIFL'][:][ring-1]
                        psin = psin[psin != 0].mean()
                        ring_label = r'Ring {} ($\phi_n$ = {:.3f})'.format(ring, psin)
                    except:
                        ring_label = 'Ring {}'.format(ring)
                    axs[graph_count].text(0.6, 0.8, ring_label,
                                          transform=axs[graph_count].transAxes)
                    axs[graph_count].set_xlabel('S (m)')
                    if show_legend == 'all':
                        axs[graph_count].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.2), fancybox=True, shadow=True)
                    if show_legend == 'short':
                        axs[graph_count].legend((line1, line2), ('Thomson', 'OEDGE'),
                                                loc='upper center', ncol=3,
                                                bbox_to_anchor=(0.5, 1.1),
                                                fancybox=True, shadow=True)

                    # Reset graph_count once it's time for a new page.
                    if graph_count == ngraphs-1:
                        graph_count = 0
                        fig.tight_layout()
                        #fig.legend(ls, shots)
                        # papertype depreciated.
                        #pdf.savefig(fig, papertype='letter')
                        pdf.savefig(fig)
                        plt.close()
                    else:
                        graph_count += 1

                # Save the last of figures remaining.
                if graph_count != 0:
                    fig.tight_layout()
                    # papertype depreciated.
                    #pdf.savefig(fig, papertype='letter')
                    pdf.savefig(fig)
                    plt.close()

    def check_ts(self, ts_filename, time_cutoff=None):
        """
        This helper function just plots the locations of the Thomson scattering
        points after they have been mapped to a common flux surface. Just to
        make sure everything looks alright. Note if there is a strike point
        sweep it will make it look like the core is sweeping a range as well
        when it really isn't. Chose time_cutoff to indicate start of sweep.

        ts_filename : The file created from create_ts.
        time_cutoff : Fix a small bug in how the core data is plotted when
                       there's a strike point sweep. This would be the start of
                       the sweep time so you don't use core data that is goofed up.
        """

        # Load up a plot. Just do Te, it doesn't really matter what data.
        fig = self.plot_contour_polygon('KTEBS', normtype='log')

        # Get our TS file that has the (R, Z) locations to plot over this fig.
        df = pd.read_excel(ts_filename)

        # Potentially a big problem, but unsure yet. Sweeping allows filling out
        # a 2D divertor space, but the mapping process used in these Thomson
        # files may goof up core measurements during a sweep. This allows
        # excluding core data after a certain time.
        if time_cutoff is not None:
            drop_idx = np.logical_and(df['System'] == 'core', df['Time'] > time_cutoff).values
            df = df.iloc[~drop_idx]

        rs = df['R (m)'].values
        zs = df['Z (m)'].values

        fig.axes[0].plot(rs, zs, 'k.', ms=1)
        fig.show()

    def find_knot_on_ring(self, ring, r, z):
        """
        Helper function to find the closest knot on a ring to an (R, Z).

        Input
        ring : Ring to find knot on.
        r : R location of the point on the grid.
        z : Z location of the point on the grid.

        Output
        closest_knot : The knot on this ring closest to r, z.
        """

        dist = np.sqrt((r - self.rs[ring])**2 + (z - self.zs[ring])**2)
        closest_knot = np.where(dist == dist.min())
        return closest_knot[0][0]

    def find_ring_knot(self, r, z, return_cell=False, verbal=True):
        """
        This function will return the ring and the knot on that ring of the cell
        closest to the input (r, z). Outputs as (ring, knot) in a tuple.

        Input
        r : R location of the point on the grid.
        z : Z location of the point on the grid.

        Output
        ring, knot : The ring and knot of which cell this point is in. The
            values are 1-indexed!!!
        """

        from shapely.geometry import Point
        from shapely.geometry.polygon import Polygon

        # Store Polygon objects when they're created to avoid having to create
        # them each call.
        try:
            self.polygons
        except:
            self.polygons = {}

        # Updated to a more robust, although more time-consuming, algorithm.
        point = Point(r, z)
        count = 0
        for cell in self.mesh:

            if count == len(self.mesh):
                if verbal:
                    print("Warning: Point not in mesh!")
                if return_cell:
                    return None
                else:
                    return (None, None)

            # Check if Polygon already exists.
            try:
                poly = self.polygons[count]
            except:
                poly = Polygon(cell)
                self.polygons[count] = poly

            if poly.contains(point):
                break
            count += 1
        cell = cell.copy()
        cell.append(cell[0])
        cell = np.array(cell)

        # Check if the rings and knots have been stored already to avoid loading
        # them each time.
        try:
            self.rings
            self.knots
        except:

            # This return the rings as 1-indexed.
            self.rings = self.read_data_2d("KTEBS", scaling="Ring")
            self.knots = self.read_data_2d("KTEBS", scaling="Knot")
        if count == len(self.rings):
            if count == len(self.mesh):
                if verbal:
                    print("Warning: Point not in mesh!")
                if return_cell:
                    return None
                else:
                    return (None, None)
        ring = int(self.rings[count])
        knot = int(self.knots[count])

        # I think I ditched this method because if you were outside the grid it
        # still returned something, which made it confusing. The above is able
        # to tell if the point is in the grid or not.
        #dist = np.sqrt(np.square(r - self.rs) + np.square(z - self.zs))
        #closest_cell = np.where(dist == dist.min())

        if return_cell:
            return cell
        else:
            return (ring, knot)

    def find_ring_from_psin(self, psin):
        """
        BUG IN THIS FUNCTION: When finding the ring for a point in the core, it
        may in fact return the ring with a matching psin in the PFZ instead.
        This is not the biggest issue since we don't really compare against
        points in those regions, but should be fixed somehow. It can also be
        an issue with extended grids, where you may have multiple rings on a
        very similar psin (i.e. maybe it made it tiny little ring tucked away
        in a baffle). This latter bug in effect will just "steal" some of the
        data, not an issue if you have enough of it.

        Helper function to find the ring with the closest average psin value to
        the input psin in.

        Input
        psin : Psin value of which to find the closest ring for.

        Output
        close_ring : The closest ring to this psin.
        """

        # See if we've already calculated psin_dict.
        try:
            self.psin_dict
        except:
            self.psin_dict = {}
            psifl = np.array(self.nc.variables['PSIFL'][:])
            for ring in range(0, psifl.shape[0]):
                if np.count_nonzero(psifl[ring]) == 0:
                    continue
                else:

                    # Store the average psin value into our dictionary.
                    psin_avg = psifl[ring][psifl[ring] != 0].mean()
                    self.psin_dict[ring] = psin_avg

        # Elegant one-liner to find closest ring to psin.
        close_ring, _ = min(self.psin_dict.items(), key=lambda item: abs(item[1] - psin))
        return close_ring

    def fake_probe(self, r_start, r_end, z_start, z_end, data='Te', num_locs=100,
                   plot=None, show_plot=True, fontsize=16, charge="all",
                   verbal=True, rings_only=False, skip_t13=False):
        """
        Return data, and plot, of a mock probe. Plot is useful if the probe has
        a constant R or Z value, just choose the correct option for it.

        Input
        r_start   : R coordinate of the measurement starting point.
        r_end     : R coordinate of the measurement ending point.
        z_start   : Z coordinate of the measurement starting point.
        z_end     : Z coordinate of the measurement ending point.
        data      : One of 'Te', 'ne', 'Mach', 'Velocity', 'L OTF', 'L ITF',
                      'nz', 'ring' or 'neut_dens'.
        num_locs  :
        plot      : Either None, 'R' or 'Z' (or 'r' or 'z'), or 'psin'. If the probe is at a
                     constant R, then use 'R', likewise for 'Z'.
        show_plot : Show the plot or not (i.e. if you just want the data).
        fontsize  : Font size for the data labels.
        verbal    : Not sure this is even used.
        rings_only: This option makes it so values are returned at every ring
                      between the start and end coordinates, no duplicates.
                      Currently this takes way too long to run, so I need to
                      think if there's a way to speed things up.
        skip_t13  : Don't bother adding on any additional flow that may have
                      been specified.

        Output
        probe_dict : The data used in the plot that simulates, for example, a
            plunging Langmuir probe or something.
        """

        from tqdm import tqdm

        # Create rs and zs to get measurements at and the dictionary to store
        # things in. If rings_only, the r, z's will be stored at each new ring.
        if rings_only:

            # This assumes that the coverage is fine enought o capture all the
            # rings. If not then decrease.
            num_locs = int(5000)
            rs = np.linspace(r_start, r_end, num_locs)
            zs = np.linspace(z_start, z_end, num_locs)
            probe_dict = {"r":[], "z":[], "psin":[], data:[], "ring":[], "knot":[]}
        else:
            rs = np.linspace(r_start, r_end, num_locs)
            zs = np.linspace(z_start, z_end, num_locs)
            probe_dict = {"r":rs, "z":zs, "psin":[], data:[], "ring":[], "knot":[]}

        # If we want the Mach number (or speed), we need to do a little data
        # preprocessing first to see if the additional T13 drift option was on.
        if data in ['Mach', 'Velocity']:

            # Get the 2D data from the netCDF file.
            scaling  = 1.0 / self.qtim
            kvhs     = self.nc['KVHS'][:] * scaling

            # Array to hold data with T13 data added on (will be same as
            # kvhs if T13 was off).
            kvhs_adj = kvhs

            # See if T13 was on and additional values need to be added.
            if not skip_t13:
                try:
                    pol_opt = float(self.dat_file.split('POL DRIFT OPT')[1].split(':')[0])
                    if pol_opt == 0.0:
                        print('Poloidal drift option T13 was OFF.')

                    elif pol_opt == 1.0:
                        print('Poloidal drift option T13 was ON.')

                        # Get the relevant table for the extra drifts out of the .dat file.
                        add_data = self.dat_file.split('TABLE OF DRIFT REGION BY RING - RINGS ' + \
                                                        'WITHOUT FLOW ARE NOT LISTED\n')[1]. \
                                                        split('DRIFT')[0].split('\n')

                        # Split the data between the spaces, put into DataFrame.
                        add_data = [line.split() for line in add_data]
                        add_df = pd.DataFrame(add_data[1:-1], columns=['IR', 'Vdrift (m/s)',
                                              'S_START (m)', 'S_END (m)'], dtype=float). \
                                              set_index('IR')

                        # Loop through the KVHS data one cell at a time, and if
                        # the cell nas extra Mach flow, add it.
                        for ir in range(self.nrs):
                            for ik in range(self.nks[ir]):
                                if self.area[ir, ik] != 0.0:

                                    # If this ring has additional drifts to be added.
                                    if ir in add_df.index:

                                        # Then add the drift along the appropriate s (or knot) range.
                                        if self.kss[ir][ik] > add_df['S_START (m)'].loc[ir] and \
                                           self.kss[ir][ik] < add_df['S_END (m)'].loc[ir]:

                                           kvhs_adj[ir][ik] = kvhs[ir][ik] + add_df['Vdrift (m/s)'].loc[ir]

                except AttributeError:
                    print("Error: .dat file has not been loaded.")
                except IndexError:
                    print("Warning: Can't add on T13 data if DIVIMP is not run.")

        for i in tqdm(range(0, num_locs)):

            # Get the cell that has the data at this R, Z. If rings_only then
            # continue looping until we land in a new ring. The R and Z values
            # that are appended are instead the RS and ZS values from that location.
            ring, knot = self.find_ring_knot(rs[i], zs[i], verbal=False)
            if rings_only:
                if ring in probe_dict["ring"]:
                    continue
                if ring == None:
                    probe_dict["r"].append(np.nan)
                    probe_dict["z"].append(np.nan)
                else:
                    probe_dict["r"].append(float(self.nc["RS"][ring, knot]))
                    probe_dict["z"].append(float(self.nc["ZS"][ring, knot]))

            # Outside of grid condition.
            if ring == None:
                probe_dict["ring"].append(np.nan)
                probe_dict["knot"].append(np.nan)
                probe_dict["psin"].append(np.nan)
                probe_dict[data].append(np.nan)
                continue

            probe_dict["ring"].append(ring)
            probe_dict["knot"].append(knot)

            # find_ring_knot returns the rings as 1-indexed (annoying,
            # confusing, I know). We need to subtract 1 when using them to
            # index any arrays because of python's 0-indexing.
            ring = ring - 1
            knot = knot - 1
            psin = float(self.nc['PSIFL'][ring, knot])

            if data == 'Te':
                probe = self.nc['KTEBS'][ring, knot]
                ylabel = 'Te (eV)'

            elif data == 'ne':
                probe = self.nc['KNBS'][ring, knot]
                ylabel = 'ne (m-3)'

            elif data == 'Mach':
                # Need to calculate the sound speed to back out the Mach number.
                te = self.nc['KTEBS'][ring, knot]
                ti = self.nc['KTIBS'][ring, knot]
                mb = self.crmb * 931.49 * 10**6 / ((3*10**8)**2)
                cs = np.sqrt((te + ti) / mb)
                if cs == 0.0:
                    probe = np.nan
                else:
                    probe = kvhs_adj[ring, knot] / cs
                ylabel = 'Mach'

            elif data == 'Velocity':
                probe = kvhs_adj[ring, knot]
                ylabel = 'Velocity (m/s)'

            # Plot of the connection length to the inner target.
            elif data == 'L OTF':
                smax = self.nc['KSMAXS'][ring]
                s    = self.nc['KSS'][ring, knot]
                probe = smax - s
                ylabel = 'L ITF (m)'

            elif data == 'L ITF':
                s    = self.nc['KSS'][ring, knot]
                probe = s
                ylabel = 'L OTF (m)'

            elif data == "Erad":
                probe = self.nc["E_RAD"][ring, knot]
                ylabel = "Erad (V/m)"

            elif data == "Epol":
                probe = self.nc["E_POL"][ring, knot]
                ylabel = "Epol (V/m)"

            elif data == "nz":
                if charge == "all":
                    probe = self.nc["DDLIMS"][1:,:,:].sum(axis=0)[ring, knot]
                else:
                    probe = self.nc["DDLIMS"][charge+1:,:,:].sum(axis=0)[ring, knot]
                ylabel = "nz (m-3)"

            elif data == "ring":
                probe = ring + 1  # To go back with the 1-indexed convention.
                ylabel = "Ring"

            elif data == "neut_dens":
                probe = self.nc["PINATO"][ring, knot]
                ylabel = "Neutral Density (m-3)"

            probe_dict["psin"].append(psin)
            probe_dict[data].append(float(probe))

        # Make a plot of the data.
        if plot is not None:

            # Get correct X and Y arrays for plotting.
            if plot.lower() in ['r', 'rminrsep']:
                x = probe_dict["r"]
                xlabel = 'R (m)'

                # Not implemented.
                if plot.lower() == 'rminrsep':
                    pass

            elif plot.lower() == 'z':
                x = probe_dict["z"]
                xlabel = 'Z (m)'
            elif plot.lower() == 'psin':
                x = probe_dict["psin"]
                xlabel = 'Psin'
            y = probe_dict[data]

            if show_plot:
                fig = plt.figure()
                ax  = fig.add_subplot(111)
                ax.plot(x, y, lw=3, color='k')
                ax.set_xlabel(xlabel, fontsize=fontsize)
                ax.set_ylabel(ylabel, fontsize=fontsize)
                ax.tick_params(axis='both', labelsize=fontsize*0.75)
                fig.tight_layout()
                fig.show()

        return probe_dict

    def along_ring(self, ring, dataname, ylabel=None, charge=None, vz_mult=0.0,
        plot_it=True, remove_zeros=True):
        """
        Plot data along a specified ring. Will return the x, y data just in case
        you want it.

        Input
        ring     : The ring number to plot data for. This is 1-INDEXED. I.e.
                     the first ring = 1, not 0 (as in the DIVIMP source).
        dataname : The NetCDF variable you want the data along the ring for. Special
                    options include 'Mach' or 'Velocity' that can add on the additional
                    drift option T13 if it was on.
        ylabel   : Label for the Y-axis.
        charge   : Charge, if needed.
        vz_mult  : The multiplier to be used in FF calculations (see calculate forces).
        plot_it  : Show the plot or not.

        Output
        x, y : The data used to make the plot. Would be S, data.
        """

        # Load the 2D array of KVHS with the additional drifts added on, optional.
        def load_kvhs_adj_2d():

            # Get the 2D data from the netCDF file.
            scaling  = 1.0 / self.qtim
            kvhs     = self.nc['KVHS'][:] * scaling

            # Array to hold data with T13 data added on (will be same as
            # kvhs if T13 was off).
            kvhs_adj = kvhs

            try:
                pol_opt = float(self.dat_file.split('POL DRIFT OPT')[1].split(':')[0])
                if pol_opt == 0.0:
                    print('Poloidal drift option T13 was OFF.')

                elif pol_opt in [1.0, 2.0, 3.0]:
                    print('Poloidal drift option T13 was ON.')

                    # Get the relevant table for the extra drifts out of the .dat file.
                    add_data = self.dat_file.split('TABLE OF DRIFT REGION BY RING - RINGS ' + \
                                                    'WITHOUT FLOW ARE NOT LISTED\n')[1]. \
                                                    split('DRIFT')[0].split('\n')

                    # Split the data between the spaces, put into DataFrame.
                    add_data = [line.split() for line in add_data]
                    add_df = pd.DataFrame(add_data[1:-1], columns=['IR', 'Vdrift (m/s)',
                                          'S_START (m)', 'S_END (m)'], dtype=float). \
                                          set_index('IR')

                    # Loop through the KVHS data one cell at a time, and if
                    # the cell has extra Mach flow, add it.
                    for ir in range(self.nrs):
                        for ik in range(self.nks[ir]):
                            if self.area[ir, ik] != 0.0:

                                # If this ring has additional drifts to be added.
                                if ir+1 in add_df.index:  # ir is zero-indexed

                                    # Then add the drift along the appropriate s (or knot) range.
                                    if self.kss[ir][ik] > add_df['S_START (m)'].loc[ir+1] and \
                                       self.kss[ir][ik] < add_df['S_END (m)'].loc[ir+1]:

                                       kvhs_adj[ir][ik] = kvhs[ir][ik] + add_df['Vdrift (m/s)'].loc[ir+1]
                return kvhs_adj

            except AttributeError:
                print("Error: .dat file has not been loaded.")

            except IndexError:

                # Happens if DIVIMP not run, and since T13 is a DIVIMP only option,
                # it's irrelevant when just creating background.
                return kvhs_adj

        # Get the parallel to B coordinate. Not sure why this "1" is needed, but
        # these is a like extra point in this KSB data to be ignored for the test
        # cases I work with.
        x = self.nc['KSB'][:][ring-1][1:].data

        # Some translations.
        if dataname in ['KVHS - Mach', 'KVHSimp - Mach']:
            dataname = 'Mach'
        elif dataname in ['KVHS', 'KVHSimp']:
            dataname = 'Velocity'

        # If we want the Mach number (or speed), we need to do a little data
        # preprocessing first to see if the additional T13 drift option was on.
        if dataname in ['Mach', 'Velocity']:

            kvhs_adj = load_kvhs_adj_2d()

            """
            # Get the 2D data from the netCDF file.
            scaling  = 1.0 / self.qtim
            kvhs     = self.nc['KVHS'][:] * scaling

            # Array to hold data with T13 data added on (will be same as
            # kvhs if T13 was off).
            kvhs_adj = kvhs

            try:
                pol_opt = float(self.dat_file.split('POL DRIFT OPT')[1].split(':')[0])
                if pol_opt == 0.0:
                    print('Poloidal drift option T13 was OFF.')

                elif pol_opt in [1.0, 2.0, 3.0]:
                    print('Poloidal drift option T13 was ON.')

                    # Get the relevant table for the extra drifts out of the .dat file.
                    add_data = self.dat_file.split('TABLE OF DRIFT REGION BY RING - RINGS ' + \
                                                    'WITHOUT FLOW ARE NOT LISTED\n')[1]. \
                                                    split('DRIFT')[0].split('\n')

                    # Split the data between the spaces, put into DataFrame.
                    add_data = [line.split() for line in add_data]
                    add_df = pd.DataFrame(add_data[1:-1], columns=['IR', 'Vdrift (m/s)',
                                          'S_START (m)', 'S_END (m)'], dtype=float64). \
                                          set_index('IR')

                    # Loop through the KVHS data one cell at a time, and if
                    # the cell nas extra Mach flow, add it.
                    for ir in range(self.nrs):
                        for ik in range(self.nks[ir]):
                            if self.area[ir, ik] != 0.0:

                                # If this ring has additional drifts to be added.
                                if ir in add_df.index:

                                    # Then add the drift along the appropriate s (or knot) range.
                                    if self.kss[ir][ik] > add_df['S_START (m)'].loc[ir] and \
                                       self.kss[ir][ik] < add_df['S_END (m)'].loc[ir]:

                                       kvhs_adj[ir][ik] = kvhs[ir][ik] + add_df['Vdrift (m/s)'].loc[ir]

            except AttributeError:
                print("Error: .dat file has not been loaded.")

            except IndexError:

                # Happens if DIVIMP not run, and since T13 is a DIVIMP only option,
                # it's irrelevant when jsut creating background.
                pass
            """

            # Finally put it into the y value of the ring we want.
            if dataname == 'Mach':

                # Need to calculate the sound speed to back out the Mach number.
                te = self.nc['KTEBS'][:][ring-1]
                ti = self.nc['KTIBS'][:][ring-1]
                mb = self.crmb * 931.49 * 10**6 / ((3*10**8)**2)
                cs = np.sqrt((te + ti) / mb)
                y  = kvhs_adj[ring-1] / cs
                #ylabel = 'Mach'

            elif dataname == 'Velocity':
                y  = kvhs_adj[ring-1]
                #ylabel = 'Velocity (m/s)'

        # If you want the impurity density, make sure you can handle the sum
        # across all the charge states or a specific charge state. Note we do
        # not do charge - 1 here since the zero index is actually neutrals.
        elif dataname =='DDLIMS':
            scaling = self.absfac
            if charge == 'all':
                #y = self.nc[dataname][:-1].sum(axis=0)[ring-1] * scaling
                y = self.nc[dataname][1:].sum(axis=0)[ring-1] * scaling
            else:
                #y = self.nc[dataname][charge][ring-1] * scaling
                y = self.nc[dataname][1+charge][ring-1] * scaling

        # Pull from the forces data if you want a force plot.
        elif dataname.lower() in ['ff', 'fig', 'feg', 'fpg', 'fe', 'fnet', 'ff']:

            # Some constants and required factors.
            col_log = 15
            qe      = 1.602E-19
            amu_kg  = 1.66E-27
            emi = 1.602192E-19 / 1.672614E-27
            #fact = np.power(self.qtim, 2) * qe / amu_kg
            fact = np.power(self.qtim, 2) * emi / self.crmi

            if dataname.lower() in ['fig', 'fnet']:

                # Need charge to calculate FiG.
                if charge == None:
                    print("Error: Must supply a charge state to calculate FiG")
                    return None

                mu = self.crmi / (self.crmi + self.crmb)
                beta = 3 * (mu + 5 * np.sqrt(2) * np.power(charge, 2) * \
                       (1.1 * np.power(mu, 5/2)- 0.35 * np.power(mu, 3/2)) - 1) / \
                       (2.6 - 2 * mu + 5.4 * np.power(mu, 2))
                print("FiG: Beta = {:.2f}".format(beta))

                # Again, must scale the gradient with this scaling factor.
                #kfigs = self.read_data_2d('KFIGS', scaling = qe / fact)

                kfigs = self.nc['KFIGS'][:][ring-1].data * qe / fact
                #kfigs = self.nc['KFIGS'][:][ring].data / fact

                # Calculate the force.
                fig = beta * kfigs
                y = np.array(fig, dtype=float64)

            if dataname.lower() in ['ff', 'fnet']:

                # Need charge to calculate FF.
                if charge == None:
                    print("Error: Must supply a charge state to calculate FF")
                    return None

                ti = self.nc['KTIBS'][:][ring-1].data
                ne = self.nc['KNBS'][:][ring-1].data

                # Slowing down time.
                tau_s = 1.47E13 * self.crmi * ti * np.sqrt(ti / self.crmb) / \
                        ((1 + self.crmb / self.crmi) * ne * np.power(charge, 2) * col_log)

                # TODO: Currently will assume impurity velocity is some fraction of the
                # background velocity (default zero), though this is obviously not true
                # and a better implementation would use the real values wherever they are.

                kvhs_adj = load_kvhs_adj_2d()

                # Don't forget KVHS is scaled by 1/QTIM to get m/s.
                #scaling = 1.0 / self.qtim
                #vi = self.read_data_2d('KVHS', scaling=scaling)
                #vi = self.nc['KVHS'][:][ring].data * scaling
                vi = kvhs_adj[:][ring-1]

                vz = vz_mult * vi

                # Calculate the force.
                ff = self.crmi * amu_kg * (vi - vz) / tau_s
                y = np.array(ff, dtype=float64)

            if dataname.lower() in ['fe', 'fnet']:

                # This also gets scaled by a factor as well in Jake's code.
                #e_pol = self.read_data_2d('E_POL', scaling = qe / fact
                e_pol = self.nc['E_POL'][:][ring-1].data * qe / fact
                fe = charge * qe * e_pol
                y = np.array(fe, dtype=float64)

            if dataname.lower() in ['feg', 'fnet']:

                alpha = 0.71 * np.power(charge, 2)

                # The electron temperature gradient from the code. I don't understand
                # this scaling factor, but Jake uses it to get into presumably eV/m.
                #kfegs = self.read_data_2d('KFEGS', scaling = qe / fact)
                kfegs = self.nc['KFEGS'][:][ring-1] * qe / fact

                # Calculate the force.
                feg = alpha * kfegs
                y = np.array(feg, dtype=float64)

            if dataname.lower() in ['fpg', 'fnet']:

                # Parallel collisional diffusion time.
                tau_par = 1.47E13 * self.crmi * ti * np.sqrt(ti / self.crmb) / \
                          (ne * np.power(charge, 2) * col_log)

                print("Warning: FPG not implemented.")
                fpg = np.full(len(x), 0)
                y = fpg

            if dataname.lower() == 'fnet':
                y = fig + ff + feg + fe + fpg

        else:
            # Get the data for this ring.
            if charge == None:
                y = np.array(self.nc[dataname][:][ring-1].data, dtype=float64)
            else:
                y = np.array(self.nc[dataname][:][charge-1][ring-1].data, dtype=float64)

        # Remove any (0, 0) data points that may occur due to fortran being fortran.
        drop_idx = np.array([], dtype=int)
        for i in range(0, len(x)):
             #print("{}: {} {}".format(i, x[i]==0.0, y[i]==0.0))
             if x[i] == 0.0 and y[i] == 0.0:
                 drop_idx = np.append(drop_idx, i)

        #print("{} drop_idx: {}".format(dataname, drop_idx))
        if remove_zeros:
            x = np.delete(x, drop_idx)
            y = np.delete(y, drop_idx)

        if plot_it:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x, y, 'k-')
            ax.set_xlabel('S (m)', fontsize=18)
            ax.set_ylabel(ylabel, fontsize=18)
            ax.tick_params(axis='both', labelsize=18*0.75)
            fig.tight_layout()
            fig.show()

        return x, y

    def input_for_midpoint_shift(self, knot):
        """
        This function will print out data that you can input into the DIVIMP
        input file for the shifted midpoint option G55. Knot is what knot you
        want the midpoint to be shifted to. Can rerun this function multiple
        times to see the knot until you find one you're happy with.

        knot : Knot you want the midpoint shifted to for each ring in just
                the main SOL.
        """

        # First plot the knot to make sure it's the one you want the midpoint
        # to be at.
        self.plot_contour_polygon('Knot', vmin=knot-1, vmax=knot+1, lut=3)

        # Load the Snorm data in 2D format.
        snorm = self.nc['KSB'][:] / self.ksmaxs[:, None]
        snorm = np.abs(snorm + snorm / 1.0 - 1.0)

        if self.nc["ZXP"][:].data > 0:
            soffset_mult = -1.0
            print("Plasma is USN... double check output after running DIVIMP...")
        else:
            soffset_mult = 1.0
            print("Plasma is LSN...")

        print("Lines for input: {}".format(self.irwall-self.irsep+1))
        print("Copy/paste under option G55 in input file:")

        # Print out the offset for each ring one at a time.
        for ring in range(self.irsep, self.irwall+1):
            soffset = soffset_mult * snorm[ring-1][knot-1] / 2.0  # Ring and knot are 1-indexed here
            print("{:8}{:8}{:8}{:10.6f}".format(1, ring, ring, soffset))


    def collector_probe(self, r1, z1, r2, z2, showplot=True, numlocs=20):
        """
        In the same vein as the OUT routine, try and approximate what a
        collector probe would see from a DIVIMP simulation. This is a considered
        subpar to a full-fledged 3DLIM workflow.

        Inputs
        r1, z1, r2, z2 (float): The start and end coordinates of the probe.
        """

        # Find the rings that the start and end location span. That's how many
        # data points we will have.
        #ir1, ik1 = self.find_ring_knot(r1, z1)
        #ir2, ik2 = self.find_ring_knot(r2, z2)

        numlocs = numlocs
        cprs = np.linspace(r1, r2, numlocs)
        cpzs = np.linspace(z1, z2, numlocs)

        cp = {"r":[], "z":[], "ir":[], "ik":[], "nz_avg1":[], "nz_avg2":[],
            "te":[], "ti":[], "cs":[], "flux1":[], "flux2":[], "rmrsomp":[]}
        for i in range(0, numlocs):
            print(i)

            # First, find the intersection ring, knots of each CP location.
            # If we've already computed for this location, no point in doing it
            # again so skip.
            ir, ik = self.find_ring_knot(cprs[i], cpzs[i])
            if ir in cp["ir"] and ik in cp["ik"]:
                continue
            if type(ik) == type(None):
                continue
            cp["r"].append(cprs[i])
            cp["z"].append(cpzs[i])
            cp["ir"].append(ir)
            cp["ik"].append(ik)

            # Store this location's R-Rsep @ OMP value. The values are !-indexed,
            # but we need to make it 0-indexed for indexing arrays.
            rmrsomp = float(self.nc["MIDIST"][1][ir-1])
            cp["rmrsomp"].append(rmrsomp)

            # The CP will collect from half the field line (the whole Lcoll
            # thing from Peter's 1985 paper can be added later). The average
            # density each side probe "sees" is then that over half the field
            # line in each direction.
            s = self.nc["KSS"][ir-1].data
            #mask = s > 0
            #s = s[mask]
            smax = max(s)
            sint = s[ik-1]
            scoll1 = sint + (smax-sint) / 2
            scoll2 = sint - (smax-sint) / 2

            # If either are negative, just means we exceeded the field line
            # length and to default to half the distance.
            if scoll1 < 0:
                scoll1 = sint / 2
            if scoll2 < 0:
                scoll2 = sint / 2

            # For side 1 (in the direction of increasing s), we integrate the
            # impurity density from sint to scoll1. For side 2 (direction of
            # negative s), we integrate from scoll2 to sint.
            scoll1_ik = np.argmin(np.abs(s - scoll1))
            scoll2_ik = np.argmin(np.abs(s - scoll2))

            # along_ring accepts the ring number as 1-indexed, we are working
            # with 0-indexed right now. If the scoll index = len(nz), then we are
            # at the end of the ring so there is no impurity to be collected.
            _, nz = self.along_ring(ir, "DDLIMS", charge="all", plot_it=False)
            if scoll1_ik != len(nz):
                nz_avg1 = np.trapz(nz[ik-1:scoll1_ik], s[ik-1:scoll1_ik]) / (scoll1 - sint)
            else:
                nz_avg1 = 0.0
            if scoll2_ik != len(nz):
                nz_avg2 = np.trapz(nz[scoll2_ik-1:ik], s[scoll2_ik-1:ik]) / (sint - scoll2)
            else:
                nz_avg2 = 0.0
                print(" Doink")
            cp["nz_avg1"].append(nz_avg1)
            cp["nz_avg2"].append(nz_avg2)

            # For the sound speed, which is needed for the flux (nz * cs), we
            # just use the sound speed at the probe location. Bit ad-hoc.
            te = float(self.nc["KTEBS"][ir-1, ik-1])
            ti = float(self.nc["KTIBS"][ir-1, ik-1])
            mb = self.crmb * 931.49 * 10**6 / ((3*10**8)**2)
            cs = np.sqrt((te + ti) / mb)
            cp["te"].append(te)
            cp["ti"].append(ti)
            cp["cs"].append(cs)

            cp["flux1"].append(nz_avg1 * cs)
            cp["flux2"].append(nz_avg2 * cs)

        if showplot:
            fig, ax = plt.subplots()
            ax.plot(cp["r"], cp["flux1"], label=1, color="tab:red", lw=2)
            ax.plot(cp["r"], cp["flux2"], label=2, color="tab:purple", lw=2)
            ax.legend()
            ax.set_xlabel("R (m)")
            ax.set_ylabel("Impurity Flux (m-2 s-1)")
            fig.tight_layout()
            fig.show()

        return cp


# ------------------------
# End of OedgePlots class.
#-------------------------


def combine_imp_plot(ncpaths=None, **args):
    """
    This standalone function will create a plot of DDLIMS by combining multiple
    runs together. This is useful where if it's quicker to run multiple small
    cases instead of one huge one.

    Inputs
    ncpaths : If None will prompt you to select the paths.
    **args : Anything else that goes into the plot parameters in plot_contour_polygon.

    Outputs
    ncpaths : Can pass back into function to save time redoing the plot.
    """

    # Needed for the file prompt.
    import tkinter as tk
    from tkinter import filedialog
    root = tk.Tk(); root.withdraw()

    # First ask the user to select all their runs to be combined.
    if ncpaths == None:
        print("Select runs to combine (CTRL + Click them)...")
        ncpaths = []
        try:

            # Shanw's run directory because he wrote the code and can do this.
            initialdir = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/'
            ncpaths = tk.filedialog.askopenfilenames(
                  filetypes=(('NetCDF files', '*.nc'),),
                  initialdir=initialdir)
        except:
            ncpaths = tk.filedialog.askopenfilenames(
                  filetypes=(('NetCDF files', '*.nc'),))

    # Load OedgePlots objects of each of these, and grab the impurity data.
    ops = []
    ddlims = []
    for i in range(0, len(ncpaths)):
        ops.append(OedgePlots(ncpaths[i]))
        ddlims.append(ops[i].read_data_2d('DDLIMS', charge='all', scaling=ops[i].absfac))

    # Find the average of the impurity data, ignoring zeros by converting them
    # to nans, and then back again to 0's at the end.
    ddlims = np.array(ddlims)
    ddlims[ddlims == 0] = np.nan
    ddlims_avg = np.nanmean(ddlims, axis=0)
    ddlims_avg[np.isnan(ddlims_avg)] = 0

    # Supply to the plot function. Doesn't matter which op we use since we're
    # using the own_data option. Supply random dataname since its required, but
    # it gets overwritten by our own_data anyways.
    ops[0].plot_contour_polygon('KTIBS', own_data=ddlims_avg, normtype='log',
                cbar_label='Tungsten Density (m-3)', cmap='nipy_spectral', **args)

    return ncpaths, ops
