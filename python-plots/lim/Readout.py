import numpy              as np
import pandas             as pd
import matplotlib.pyplot  as plt
import matplotlib.patches as patches
import matplotlib         as mpl
from matplotlib       import colors
from collections      import OrderedDict
from tkinter          import filedialog, Tk
from scipy.optimize   import curve_fit
import netCDF4


# Some plot properties to make them a bit nicer.
plt.ion()
plt.rcParams['font.family'] = 'serif'
fontsize = 12
ms = 2
lw = 5
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the tableau20 RGBs to numbers between (0,1) since this is how mpl accepts them.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)


class Readout:
        """
        Class object that holds a figure of a grid of plots from the netCDF
        output of 3DLIM. Example usage in controlling script lim_readout.py.
        """

        def __init__(self, netcdf_file=None, dat_file=None, lim_file=None,
                     figsize=(15,10), grid_shape=(3,3)):
            """
            netcdf_file: Path to the 3DLIM netCDF file. If None is entered then
                         it will use 'colprobe-test-m2.nc' as a default, which
                         is convienent for testing.
            dat_file:    Path to the casename.dat file. This has things like the
                         duration of the run in it.
            lim_file:    Path to the casename.lim file. This has the forces in it,
                         among many, many other things.
            figsize:     Size of the figure to hold all the plots. The default of
                         (15, 10) is a good size.
            grid_shape:  The shape of the grid of plots. Change if you want to
                         add more plots or whatever.
            """

            # Create the master figure to hold all the plots.
            self.master_fig = plt.figure(figsize=figsize)

            # If no netCDF file given, just use this test one and the dat file.
            if netcdf_file is None:
                self.netcdf = netCDF4.Dataset('colprobe-test-m2.nc')
                with open('colprobe-test-m2.dat') as f:
                    self.dat = f.read()
            else:
                self.netcdf = netCDF4.Dataset(netcdf_file)

            # Get the .dat file info as well, if supplied, otherwise it's None.
            if dat_file:
                with open(dat_file) as f:
                    self.dat = f.read()
            else:
                self.dat = None

            # Same with .lim file.
            if lim_file:
                with open(lim_file) as f:
                    self.lim = f.read()
            else:
                self.lim = None

            # Create figure with array of empty plots.
            for plot_num in range(1, grid_shape[0] * grid_shape[1] + 1):
                self.master_fig.add_subplot(grid_shape[0], grid_shape[1], plot_num)

        def print_readout(self):
            """
            Output a table with relevant info from the netcdf file.
            """

            # Let's just put everything we want into a dict so printing is easy.
            output = OrderedDict()
            output['3DLIM Version'] = self.netcdf['VERSION'][:].data.tostring().decode()
            output['Title'] = self.netcdf['TITLE'][:].data.tostring().decode()
            output['File'] = self.netcdf['JOB'][:].data.tostring().decode().split(' ')[0]
            #output['Particles'] = format(self.netcdf['MAXIMP'][:].data, ',')
            output['Conn. Length'] = self.netcdf['CL'][:].data

            if self.dat:
                # Get the total CPU time used.
                try:
                    time = int(self.dat.split('TOTAL CPU TIME USED     (S)')[1].split('\n')[0])
                    output['Time'] = str(time) + 's (' + format(time/3600, '.2f') + ' hours)'
                except:
                    pass

                try:
                    num = int(self.dat.split('NO OF IMPURITY IONS TO FOLLOW')[1].split('\n')[0])
                    output['No. Imp. Ions'] = "{:,}".format(num)
                except:
                    pass

            # Find longest output for formatting.
            pad = 0
            for val in output.values():
                if len(str(val)) > pad:
                    pad = len(str(val))

            # Printing commands.
            num_stars = 2 + 15 + 2 + pad
            print("\n" + "*"*num_stars)
            for key, val in output.items():
                print("* {:15}{:<{pad}} *".format(key, val, pad=pad))
            print("*"*num_stars)

            # Also while we're here put the figure title as the filename.
            self.master_fig.subplots_adjust(top=0.60)
            self.master_fig.suptitle(output['Title'], fontsize=26)

        def centerline(self, plot_num, mult_runs=False, log=False, fit_exp=False):
            """
            Plot the ITF and OTF deposition along the centerlines.

            plot_num: Location in grid to place this plot. I.e. if the grid_shape
                      is (3,3), then enter a number between 0-8, where the locations
                      are labelled left to right.
            """

            #The deposition array.
            #dep_arr = np.array(self.netcdf.variables['NERODS3'][0] * -1)
            dep_arr = self.get_dep_array(mult_runs)

            # Location of each P bin, and its width. Currently they all have the same width,
            # but it may end up such that there are custom widths so we leave it like this.
            ps     = np.array(self.netcdf.variables['PS'][:].data)
            pwids  = np.array(self.netcdf.variables['PWIDS'][:].data)

            # Array of poloidal locations (i.e. the center of each P bin).
            pol_locs = ps - pwids/2.0

            # Drop last row since it's garbage.
            dep_arr = dep_arr[:-1, :]
            pol_locs = pol_locs[:-1]

            # Distance cell centers along surface (i.e. the radial locations).
            rad_locs = np.array(self.netcdf.variables['ODOUTS'][:].data)

            # Get the centerline index (or closest to it).
            cline = np.abs(pol_locs).min()

            # Index the deposition array at the centerline for plotting.
            itf_x = rad_locs[np.where(rad_locs > 0.0)[0]]
            itf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
            otf_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
            otf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]

            # Plotting commands.
            ax = self.master_fig.axes[plot_num]
            if log:
                ax.semilogy(itf_x*100, itf_y, '-', label='ITF', ms=ms, color=tableau20[6])
                ax.semilogy(otf_x*100, otf_y, '-', label='OTF', ms=ms, color=tableau20[8])
            else:
                ax.plot(itf_x*100, itf_y, '-', label='ITF', ms=ms, color=tableau20[6])
                ax.plot(otf_x*100, otf_y, '-', label='OTF', ms=ms, color=tableau20[8])
            ax.legend(fontsize=fontsize)
            ax.set_xlabel('Distance along probe (cm)', fontsize=fontsize)
            ax.set_ylabel('Deposition (arbitrary units)', fontsize=fontsize)
            ax.set_xlim([0, 10])
            ax.set_ylim([0,None])

            # Option to perform an exponential fit to the data.
            if fit_exp:
                def exp_fit(x, a, b):
                    return a * np.exp(-b * x)

                popt_itf, pcov_itf = curve_fit(exp_fit, itf_x, itf_y, maxfev=5000)
                popt_otf, pcov_otf = curve_fit(exp_fit, otf_x, otf_y, maxfev=5000)

                fitx = np.linspace(0, 0.1, 100)
                fity_itf = exp_fit(fitx, *popt_itf)
                fity_otf = exp_fit(fitx, *popt_otf)

                if log:
                    ax.semilogy(fitx*100, fity_itf, '--', ms=ms, color=tableau20[6])
                    ax.semilogy(fitx*100, fity_otf, '--', ms=ms, color=tableau20[8])
                else:
                    ax.plot(fitx*100, fity_itf, '--', ms=ms, color=tableau20[6])
                    ax.plot(fitx*100, fity_otf, '--', ms=ms, color=tableau20[8])

                print("Lambdas")
                print("  ITF = {:.2f}".format(1/popt_itf[1]*100))
                print("  OTF = {:.2f}".format(1/popt_otf[1]*100))

            #print("Max ITF/OTF: {:.2f}".format(itf_y.max()/otf_y.max()))
            print("Total ITF/OTF: {:.2f}".format(itf_y.sum()/otf_y.sum()))

            return {'itf_x':itf_x, 'itf_y':itf_y, 'otf_x':otf_x, 'otf_y':otf_y}

        def deposition_contour(self, plot_num, side, probe_width=0.015, rad_cutoff=0.1, mult_runs=False):
            """
            Plot the 2D tungsten distribution across the face.

            plot_num:    Location in grid to place this plot. I.e. if the grid_shape
                         is (3,3), then enter a number between 0-8, where the locations
                         are labelled left to right.
            side: Either 'ITF' or 'OTF'.
            probe_width: The half-width of the collector probe (the variable CPCO).
                         A = 0.015, B = 0.005, C = 0.0025
            rad_cutoff:  Only plot data from the tip down to rad_cutoff. Useful
                         if we want to compare to LAMS since those scans only go
                         down a certain length of the probe.

            *** To-do ***
            - Instead of entering the width, pull out CPCO(?) from the netcdf file.
                Need to figure out the points being deposited outside the expected
                probe width first though.

            """

            #The deposition array.
            #dep_arr = np.array(self.netcdf.variables['NERODS3'][0] * -1)
            dep_arr = self.get_dep_array(mult_runs)

            # Location of each P bin, and its width. Currently they all have the same width,
            # but it may end up such that there are custom widths so we leave it like this.
            ps     = np.array(self.netcdf.variables['PS'][:].data)
            pwids  = np.array(self.netcdf.variables['PWIDS'][:].data)

            # Array of poloidal locations (i.e. the center of each P bin).
            pol_locs = ps - pwids/2.0

            # Drop last row since it's garbage.
            dep_arr = dep_arr[:-1, :]
            pol_locs = pol_locs[:-1]

            # Distance cell centers along surface (i.e. the radial locations).
            rad_locs = np.array(self.netcdf.variables['ODOUTS'][:].data)

            # Remove data beyond rad_cutoff.
            idx = np.where(np.abs(rad_locs)<rad_cutoff)[0]
            rad_locs = rad_locs[idx]
            dep_arr = dep_arr[:, idx]

            # Get only positive values of rad_locs for ITF...
            idx = np.where(rad_locs > 0.0)[0]
            X_itf, Y_itf = np.meshgrid(rad_locs[idx], pol_locs)
            Z_itf = dep_arr[:, idx]

            # ... negative for OTF.
            idx = np.where(rad_locs < 0.0)[0]
            X_otf, Y_otf = np.meshgrid(np.abs(rad_locs[idx][::-1]), pol_locs)
            Z_otf = dep_arr[:, idx][:, ::-1]

            # Make the levels for the contour plot out of whichever side has the max deposition.
            if Z_itf.max() > Z_otf.max():
                levels = np.linspace(0, Z_itf.max(), 15)
            else:
                levels = np.linspace(0, Z_otf.max(), 15)

            # Plotting commands.
            if side == 'ITF':
                X = X_itf; Y = Y_itf; Z = Z_itf
            else:
                X = X_otf; Y = Y_otf; Z = Z_otf

            ax = self.master_fig.axes[plot_num]
            ax.contourf(X*100, Y*100, Z, levels=levels, cmap='Reds')
            ax.set_xlabel('Distance along probe (cm)', fontsize=fontsize)
            ax.set_ylabel('Z location (cm)', fontsize=fontsize)
            ax.set_ylim([-probe_width*100, probe_width*100])
            props = dict(facecolor='white')
            ax.text(0.75, 0.85, side, bbox=props, fontsize=fontsize*1.5, transform=ax.transAxes)

            # Print out the total amount collected.
            #print("Total W Deposited ({:}): {:.2f}".format(side, Z.sum()))

        def velocity_contour_pol(self, pol_slice=0):
            """
            Plot the 2D distribution of the (tungsten? plasma?) velocity at a
            poloidal slice.

            pol_slice: The poloidal coordinate to get a velocity plot in (R, B) space.
            """
            pass

        def velocity_contour_par(self, par_slice=0):
            """
            Plot the 2D distribution of the (tungsten? plasma?) velocity at a
            parallel (to B) slice.

            par_slice: The parallel coordinate to get a velocity plot in (R, P) space.
            """
            pass

        def te_plot(self):
            """
            Plot the input Te (which is at the midplane?).
            """
            pass

        def ne_plot(self):
            """
            Plot the input ne (which is at the midplane?).
            """
            pass

        def te_contour(self, plot_num):
            """
            Plot the 2D background electron plasma temperature.

            plot_num: Location in grid to place this plot. I.e. if the grid_shape
                      is (3,3), then enter a number between 0-8, where the locations
                      are labelled left to right.
            """

            # Get the connection length to restrict the plot between the two absorbing surfaces.
            cl = float(self.netcdf['CL'][:].data)

            # Same with the location of the plasma center (the top of the box).
            ca = float(self.netcdf['CA'][:].data)

            # Get the X and Y grid data.
            x = self.netcdf.variables['XOUTS'][:].data
            y = self.netcdf.variables['YOUTS'][:].data

            # 2D grid of the temperature data.
            Z = self.netcdf.variables['CTEMBS'][:].data

            # Trim the zeros from the edges of the x and y arrays, and the associated
            # data points as well. This is done to stop this data from messing up
            # the contours in the contour plot.
            xkeep_min = np.nonzero(x)[0].min()
            xkeep_max = np.nonzero(x)[0].max()
            ykeep_min = np.nonzero(y)[0].min()
            ykeep_max = np.nonzero(y)[0].max()
            x = x[xkeep_min:xkeep_max]
            y = y[ykeep_min:ykeep_max]
            Z = Z[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

            # Furthermore, trim the data off that is beyond CL.
            ykeep_cl = np.where(np.abs(y) < cl)[0]
            y = y[ykeep_cl]
            Z = Z[ykeep_cl, :]

            # Replace zeros in Z with just the smallest density value. Again to
            # stop all these zeros from messing up the contour levels.
            try:
                Zmin = np.partition(np.unique(Z), 1)[1]
                Z = np.clip(Z, Zmin, None)
            except:
                pass

            # Create grid for plotting. Note we swap definitions for x and y since
            # we want the x-axis in the plot to be the parallel direction (it just
            # looks better that way).
            Y, X = np.meshgrid(x, y)

            # Plotting commands.
            ax = self.master_fig.axes[plot_num]
            cont = ax.contourf(X, Y, Z, cmap='magma', levels=10)
            ax.set_xlim([-cl, cl])
            #ax.set_ylim([None, ca])
            ax.set_ylim([None, 0.01])  # Contour weird near edge.
            ax.set_xlabel('Parallel (m)', fontsize=fontsize)
            ax.set_ylabel('Radial (m)', fontsize=fontsize)
            cbar = self.master_fig.colorbar(cont, ax=ax)
            cbar.set_label('Background Te (eV)')

        def ne_contour(self, plot_num):
            """
            Plot the 2D background plasma density.

            plot_num: Location in grid to place this plot. I.e. if the grid_shape
                      is (3,3), then enter a number between 0-8, where the locations
                      are labelled left to right.
            """
            # Get the connection length to restrict the plot between the two absorbing surfaces.
            cl = float(self.netcdf['CL'][:].data)
            # Same with the location of the plasma center (the top of the box)
            ca = float(self.netcdf['CA'][:].data)

            # Get the X and Y grid data.
            x = self.netcdf.variables['XOUTS'][:].data
            y = self.netcdf.variables['YOUTS'][:].data

            # 2D grid of the temperature data.
            Z = self.netcdf.variables['CRNBS'][:].data

            # Trim the zeros from the edges of the x and y arrays, and the associated
            # data points as well. This is done to stop this data from messing up
            # the contours in the contour plot.
            xkeep_min = np.nonzero(x)[0].min()
            xkeep_max = np.nonzero(x)[0].max()
            ykeep_min = np.nonzero(y)[0].min()
            ykeep_max = np.nonzero(y)[0].max()
            x = x[xkeep_min:xkeep_max]
            y = y[ykeep_min:ykeep_max]
            Z = Z[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

            # Furthermore, trim the data off that is beyond CL.
            ykeep_cl = np.where(np.abs(y) < cl)[0]
            y = y[ykeep_cl]
            Z = Z[ykeep_cl, :]

            # Replace zeros in Z with just the smallest density value. Again to
            # stop all these zeros from messing up the contour levels.
            Zmin = np.partition(np.unique(Z), 1)[1]
            Z = np.clip(Z, Zmin, None)

            # Create grid for plotting. Note we swap definitions for x and y since
            # we want the x-axis in the plot to be the parallel direction (it just
            # looks better that way).
            Y, X = np.meshgrid(x, y)

            # Plotting commands.
            ax = self.master_fig.axes[plot_num]

            # Create our own levels since the automatic ones are bad.
            lev_exp = np.arange(np.floor(np.log10(Z.min())-1), np.ceil(np.log10(Z.max())+1), 0.25)
            levs = np.power(10, lev_exp)

            cont = ax.contourf(X, Y, Z, cmap='magma', levels=levs, norm=colors.LogNorm())
            ax.set_xlim([-cl, cl])
            #ax.set_ylim([None, ca])
            ax.set_ylim([None, 0.01])  # Contour weird near edge.
            ax.set_xlabel('Parallel (m)', fontsize=fontsize)
            ax.set_ylabel('Radial (m)', fontsize=fontsize)
            cbar = self.master_fig.colorbar(cont, ax=ax)
            cbar.set_label('Background ne (m-3)')

        def avg_imp_vely(self, plot_num):
            """
            SVYBAR: Average impurity velocity at X coordinates in QXS.

            plot_num: Location in grid to place this plot. I.e. if the grid_shape
                      is (3,3), then enter a number between 0-8, where the locations
                      are labelled left to right.
            """

            # Grab the data.
            x = self.netcdf.variables['QXS'][:].data
            y = self.netcdf.variables['SVYBAR'][:].data

            # Plotting commands.
            ax = self.master_fig.axes[plot_num]
            ax.plot(x, y, '.', ms=ms, color=tableau20[6])
            ax.set_xlabel('Radial coordinates (m)', fontsize=fontsize)
            ax.set_ylabel('Average Y imp. vel. (m/s)', fontsize=fontsize)

        def avg_pol_profiles(self, plot_num, probe_width=0.015, rad_cutoff=0.5):
            """
            Plot the average poloidal profiles for each side. Mainly to see if
            deposition peaks on the edges.

            plot_num:    Location in grid to place this plot. I.e. if the grid_shape
                         is (3,3), then enter a number between 0-8, where the locations
                         are labelled left to right.
            probe_width: The half-width of the collector probe (the variable CPCO).
                         A = 0.015, B = 0.005, C = 0.0025
            rad_cutoff:  Only plot data from the tip down to rad_cutoff. Useful
                         if we want to compare to LAMS since those scans only go
                         down a certain length of the probe.
            """

            # Code copied from above function, deposition_contour. See for comments.
            dep_arr = np.array(self.netcdf.variables['NERODS3'][0] * -1)
            ps     = np.array(self.netcdf.variables['PS'][:].data)
            pwids  = np.array(self.netcdf.variables['PWIDS'][:].data)
            pol_locs = ps - pwids/2.0
            dep_arr = dep_arr[:-1, :]
            pol_locs = pol_locs[:-1]
            rad_locs = np.array(self.netcdf.variables['ODOUTS'][:].data)
            idx = np.where(np.abs(rad_locs)<rad_cutoff)[0]
            rad_locs = rad_locs[idx]
            dep_arr = dep_arr[:, idx]
            idx = np.where(rad_locs > 0.0)[0]
            X_itf, Y_itf = np.meshgrid(rad_locs[idx], pol_locs)
            Z_itf = dep_arr[:, idx]
            idx = np.where(rad_locs < 0.0)[0]
            X_otf, Y_otf = np.meshgrid(np.abs(rad_locs[idx][::-1]), pol_locs)
            Z_otf = dep_arr[:, idx][:, ::-1]

            # Average along the radial direction.
            avg_pol_itf = np.mean(Z_itf, 1)
            avg_pol_otf = np.mean(Z_otf, 1)

            # Get the centerline index (or closest to it).
            cline = np.abs(pol_locs).min()
            cline_idx = np.where(pol_locs == cline)[0][0]

            # Get average peaking factor for each side.
            peak1 = avg_pol_itf[:cline_idx].max() / avg_pol_itf[cline_idx]
            peak2 = avg_pol_itf[cline_idx:].max() / avg_pol_itf[cline_idx]
            itf_peak = (peak1 + peak2) / 2.0
            peak1 = avg_pol_otf[:cline_idx].max() / avg_pol_otf[cline_idx]
            peak2 = avg_pol_otf[cline_idx:].max() / avg_pol_otf[cline_idx]
            otf_peak = (peak1 + peak2) / 2.0

            print("OTF/ITF Peaking Ratio: {:.2f}".format(otf_peak/itf_peak))

            # Plotting commands.
            ax = self.master_fig.axes[plot_num]
            ax.plot(pol_locs, avg_pol_itf/avg_pol_itf.max(), label='ITF', color=tableau20[6])
            ax.plot(pol_locs, avg_pol_otf/avg_pol_otf.max(), label='OTF', color=tableau20[8])
            ax.legend(fontsize=fontsize)
            ax.set_xlabel('Poloidal (m)', fontsize=fontsize)
            ax.set_ylabel('Deposition (normalized)', fontsize=fontsize)
            ax.set_xlim([-probe_width, probe_width])

        def imp_contour_plot(self, plot_num, rmin=-0.005, rmax=0, iz_state=5):

            # Get positions of the center of each bin.
            xs       = self.netcdf.variables['XS'][:].data
            xwids    = self.netcdf.variables['XWIDS'][:].data
            rad_locs = xs-xwids/2.0
            ps       = self.netcdf.variables['PS'][:].data
            pwids    = self.netcdf.variables['PWIDS'][:].data
            pol_locs = ps-pwids/2.0
            ys       = self.netcdf.variables['YS'][:].data
            ywids    = self.netcdf.variables['YWIDS'][:].data
            par_locs = ys-ywids/2.0

            # Also mirror the par_locs to cover both sides (i.e. from -L to L instead of 0 to L).
            # Need to add a zero in the as the middle point, hence two appends.
            par_locs = np.append(np.append(-par_locs[::-1], 0), par_locs)

            # Load ddlim3 variable array of the specific ionization state.
            if type(iz_state) is list:
                # Add capability to do a range of ionization states.
                pass
            else:
                ddlim3 =self.netcdf.variables['DDLIM3'][:, iz_state, :, :].data

            # Sum over the radial range to create a 2D plot.
            sum_range = np.where(np.logical_and(rad_locs>rmin, rad_locs<rmax))[0]
            summed_ddlim3 = ddlim3[:,:,sum_range].sum(axis=2)

            # Plotting commands.
            X, Y = np.meshgrid(par_locs, pol_locs)
            Z    = summed_ddlim3
            ax   = self.master_fig.axes[plot_num]
            cont = ax.contourf(X, Y, Z)
            cbar = self.master_fig.colorbar(cont, ax=ax)
            cl   = float(self.netcdf['CL'][:].data)
            ax.set_xlim([-cl, cl])
            ax.set_xlabel('Parallel (m)', fontsize=fontsize)
            ax.set_ylabel('Poloidal (m)', fontsize=fontsize)
            cp = patches.Rectangle((-0.2,-0.015), width=0.4, height=0.03, color='k')
            ax.add_patch(cp)
            textstr = r'Integration region:' + \
                      r'\n$\mathrm{R_min}$ = ' + str(rmin) + \
                      r'\n$\mathrm{R_max}$ = ' + str(rmax)
            props = dict(facecolor='white')
            #ax.text(0.05, 0.95, textstr, bbox=props)

        def imp_contour_plot_radial(self, plot_num, pmin=-0.005, pmax=0, iz_state=5):

            # Get positions of the center of each bin.
            xs       = self.netcdf.variables['XS'][:].data
            xwids    = self.netcdf.variables['XWIDS'][:].data
            rad_locs = xs-xwids/2.0
            ps       = self.netcdf.variables['PS'][:].data
            pwids    = self.netcdf.variables['PWIDS'][:].data
            pol_locs = ps-pwids/2.0
            ys       = self.netcdf.variables['YS'][:].data
            ywids    = self.netcdf.variables['YWIDS'][:].data
            par_locs = ys-ywids/2.0

            # Also mirror the par_locs to cover both sides (i.e. from -L to L instead of 0 to L).
            # Need to add a zero in the as the middle point, hence two appends.
            par_locs = np.append(np.append(-par_locs[::-1], 0), par_locs)

            # Load ddlim3 variable array of the specific ionization state.
            if type(iz_state) is list:
                # Add capability to do a range of ionization states.
                pass
            else:
                ddlim3 =self.netcdf.variables['DDLIM3'][:, iz_state, :, :].data

            # Sum over the radial range to create a 2D plot.
            sum_range = np.where(np.logical_and(pol_locs>pmin, pol_locs<pmax))[0]
            summed_ddlim3 = ddlim3[sum_range,:,:].sum(axis=0)

            # Plotting commands.
            X, Y = np.meshgrid(par_locs, rad_locs)
            Z    = summed_ddlim3
            ax   = self.master_fig.axes[plot_num]
            cont = ax.contourf(X, Y, Z.T)
            cbar = self.master_fig.colorbar(cont, ax=ax)
            cl   = float(self.netcdf['CL'][:].data)
            ax.set_xlim([-cl, cl])
            ax.set_xlabel('Parallel (m)', fontsize=fontsize)
            ax.set_ylabel('Radial (m)', fontsize=fontsize)
            #cp = patches.Rectangle((-0.2,-0.015), width=0.4, height=0.03, color='k')
            #ax.add_patch(cp)
            textstr = r'Integration region:' + \
                      r'\n$\mathrm{P_min}$ = ' + str(pmin) + \
                      r'\n$\mathrm{P_max}$ = ' + str(pmax)
            props = dict(facecolor='white')
            #ax.text(0.05, 0.95, textstr, bbox=props)

        def force_plots(self, plot_num, rad_loc=-0.01, cl=9.9, separate_plot=True):

            # First, grab that while big force table, splitting it at the start
            # of each force table for each radial location (i.e. ix location).
            lim_sfs = self.lim.split('Static forces')[1:]

            # Fix the last element so it doesn't catch everything after ([:-2]
            # to ignore an extra \n and space that bugs the rest up).
            lim_sfs[-1] = lim_sfs[-1].split('***')[0][:-2]

            # Column names for the dataframe.
            col_names = ['IX', 'IY', 'XOUT', 'YOUT', 'FEG', 'FIG', 'FF', 'FE',
                         'FVH', 'FF2', 'FE2', 'FVH2', 'FTOT1', 'FTOT2', 'TEGS',
                         'TIGS', 'CFSS', 'CFVHXS', 'VP1', 'VP2', 'FFB', 'FEB',
                         'CVHYS', 'CEYS', 'TE', 'TI', 'NE', 'VELB']

            # List to hold all the force dataframes.
            dflist = []
            for sf in lim_sfs:

                # Split up the rows for this radial location.
                foo = sf.split('\n')[1:]

                # Split up each entry in each row (each row is one big string
                # at this point).
                foo = [bar.split() for bar in foo]

                # Put into dataframe and append to list with all of them.
                df = pd.DataFrame(foo, columns=col_names)
                dflist.append(df)

            # Create a single multidimensional df out of these so it's all easier
            # to work with.
            big_df = pd.concat(dflist, keys=np.arange(1, len(dflist)))

            # Plot at a location near the tip of the probe (say R=-0.01).
            for idx in np.unique(big_df.index.get_level_values(0).values):
                if np.float(big_df.loc[idx]['XOUT'][0]) > rad_loc:
                    break

            # Get values from dataframe for plotting.
            x  = np.array(big_df.loc[idx]['YOUT'].values, dtype=np.float64)
            y1 = np.array(big_df.loc[idx]['FTOT1'].values, dtype=np.float64)
            y2 = np.array(big_df.loc[idx]['FTOT2'].values, dtype=np.float64)

            # Remove nans.
            x = x[~np.isnan(x)]
            y1 = y1[~np.isnan(y1)]
            y2 = y2[~np.isnan(y2)]

            # Only want values between -cl and cl.
            valid_idx = np.where(np.logical_and(x > -cl, x < cl))
            x  = x[valid_idx]
            y1 = y1[valid_idx]
            y2 = y2[valid_idx]

            if plot_num == 99:
                pass
            else:
                ax = self.master_fig.axes[plot_num]
                ax.plot(x, y1, '-', color=tableau20[6], label='FTOT1')
                ax.plot(x, y2, '-', color=tableau20[8], label='FTOT2')
                ax.set_xlabel('Parallel (m)', fontsize=fontsize)
                ax.set_ylabel('Force (N?)', fontsize=fontsize)
                ax.legend(fontsize=fontsize)
                ax.axhline(0, linestyle='--', color='k')

            # If you want a separate plot made with all the forces, more detailed.
            if separate_plot:

                x     = np.array(big_df.loc[idx]['YOUT'].values,  dtype=np.float64)[:-1]
                valid_idx = np.where(np.logical_and(x > -cl, x < cl))
                x     = x[valid_idx]

                ftot1 = np.array(big_df.loc[idx]['FTOT1'].values, dtype=np.float64)[:-1][valid_idx]
                ftot2 = np.array(big_df.loc[idx]['FTOT2'].values, dtype=np.float64)[:-1][valid_idx]
                ff1   = np.array(big_df.loc[idx]['FF'].values,    dtype=np.float64)[:-1][valid_idx]
                ff2   = np.array(big_df.loc[idx]['FF2'].values,   dtype=np.float64)[:-1][valid_idx]
                feg   = np.array(big_df.loc[idx]['FEG'].values,   dtype=np.float64)[:-1][valid_idx]
                figf  = np.array(big_df.loc[idx]['FIG'].values,   dtype=np.float64)[:-1][valid_idx]
                fe    = np.array(big_df.loc[idx]['FE'].values,    dtype=np.float64)[:-1][valid_idx]
                fvh1  = np.array(big_df.loc[idx]['FVH'].values,   dtype=np.float64)[:-1][valid_idx]
                fvh2  = np.array(big_df.loc[idx]['FVH2'].values,  dtype=np.float64)[:-1][valid_idx]

                fig = plt.figure(figsize=(7,5))
                ax = fig.add_subplot(111)
                ax.plot(x, ftot1, '-',  color=tableau20[2],  label='FTOT1')
                ax.plot(x, ftot2, '--', color=tableau20[2],  label='FTOT2')
                ax.plot(x, ff1, '-',    color=tableau20[4],  label='FF1')
                ax.plot(x, ff2, '--',   color=tableau20[4],  label='FF2')
                ax.plot(x, feg, '-',    color=tableau20[6],  label='FEG')
                ax.plot(x, figf, '-',    color=tableau20[8],  label='FIG')
                ax.plot(x, fe, '-',     color=tableau20[10], label='FE')
                ax.plot(x, fvh1, '-',   color=tableau20[12], label='FVH1')
                ax.plot(x, fvh2, '--',  color=tableau20[12], label='FVH2')
                ax.legend(fontsize=fontsize)
                ax.set_xlabel('Parallel (m)')
                ax.set_ylabel('Force (N?)')
                fig.tight_layout()
                fig.show()

        def vel_plots(self, plot_num, vp, cl=10.0, separate_plot=True, vmin1=None, vmax1=None, vmin2=None, vmax2=None, clip=False):
            """
            TODO

            vp: Either 'vp1' (background plasma) or 'vp2' (disturbed plasma from CP).
            """

            # These lines are copied from the above force_plots. See them for comments.
            lim_sfs = self.lim.split('Static forces')[1:]
            lim_sfs[-1] = lim_sfs[-1].split('***')[0][:-2]
            col_names = ['IX', 'IY', 'XOUT', 'YOUT', 'FEG', 'FIG', 'FF', 'FE',
                         'FVH', 'FF2', 'FE2', 'FVH2', 'FTOT1', 'FTOT2', 'TEGS',
                         'TIGS', 'CFSS', 'CFVHXS', 'VP1', 'VP2', 'FFB', 'FEB',
                         'CVHYS', 'CEYS', 'TE', 'TI', 'NE', 'VELB']
            dflist = []
            for sf in lim_sfs:
                foo = sf.split('\n')[1:]
                foo = [bar.split() for bar in foo]
                df = pd.DataFrame(foo, columns=col_names)
                dflist.append(df)
            big_df = pd.concat(dflist, keys=np.arange(1, len(dflist)))

            # Unstack the big_df into individual velocity dfs for easier plotting.
            vp1_df = big_df['VP1'].unstack()
            vp2_df = big_df['VP2'].unstack()

            # May be a better way for this, but last column is nans, so drop them.
            vp1_df = vp1_df[vp1_df.columns[:-1]]
            vp2_df = vp2_df[vp2_df.columns[:-1]]

            # Get the x, y values into numpy arrays. Remove nans.
            x = np.unique(np.array(big_df['XOUT'].values, dtype=np.float64))
            y = np.unique(np.array(big_df['YOUT'].values, dtype=np.float64))
            x = x[~np.isnan(x)]
            y = y[~np.isnan(y)]

            # Create 2D grids for the contour plot.
            Y, X = np.meshgrid(x, y)
            Z1 = np.array(vp1_df.values.T, dtype=np.float64)
            Z2 = np.array(vp2_df.values.T, dtype=np.float64)

            # Choose the relevant Z. Cast to float.
            if vp == 'vp1':
                Z = Z1
            elif vp =='vp2':
                Z = Z2
            else:
                print("Error: vp must be either vp1 or vp2.")

            if vmin1 is None:
                vmin1 = -Z1.max()
            if vmax1 is None:
                vmax1 = Z1.max()

            if vmin2 is None:
                vmin2 = -Z2.max()
            if vmax2 is None:
                vmax2 = Z2.max()

            if clip:
                Z1 = np.clip(Z1, vmin1, vmax1)
                Z2 = np.clip(Z2, vmin2, vmax2)

            # Plotting commands.
            ax   = self.master_fig.axes[plot_num]
            cont = ax.contourf(X, Y, Z, vmin=-Z.max(), vmax=Z.max(), cmap='coolwarm')
            cbar = self.master_fig.colorbar(cont, ax=ax)
            ax.set_xlim([-cl, cl])
            ax.set_xlabel('Parallel (m)', fontsize=fontsize)
            ax.set_ylabel('Radial (m)', fontsize=fontsize)
            cbar.set_label(vp.upper() + ' (m/s)')

            if separate_plot:

                # Bounds needed for the colorbar. Will just do 10 levels.
                bounds1 = np.linspace(vmin1, vmax1, 10)
                bounds2 = np.linspace(vmin2, vmax2, 10)

                fig = plt.figure()
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
                cont1 = ax1.contourf(X, Y, Z1, vmin=vmin1, vmax=vmax1, cmap='coolwarm')
                cont2 = ax2.contourf(X, Y, Z2, vmin=vmin2, vmax=vmax2, cmap='coolwarm')
                cbar1 = self.master_fig.colorbar(cont1, ax=ax1, ticks=bounds1)
                cbar2 = self.master_fig.colorbar(cont2, ax=ax2, ticks=bounds2)
                ax1.set_xlim([-cl, cl])
                ax2.set_xlim([-cl, cl])
                ax1.set_xlabel('Parallel (m)', fontsize=fontsize)
                ax2.set_xlabel('Parallel (m)', fontsize=fontsize)
                ax1.set_ylabel('Radial (m)', fontsize=fontsize)
                ax2.set_ylabel('Radial (m)', fontsize=fontsize)
                cbar1.set_label('VP1 (m/s)')
                cbar2.set_label('VP2 (m/s)')
                fig.tight_layout()
                fig.show()

        def force_plot_2d(self, vmin=None, vmax=None, force='FF', xlim=(-10,10)):

            # First, grab that while big force table, splitting it at the start
            # of each force table for each radial location (i.e. ix location).
            lim_sfs = self.lim.split('Static forces')[1:]

            # Fix the last element so it doesn't catch everything after ([:-2]
            # to ignore an extra \n and space that bugs the rest up).
            lim_sfs[-1] = lim_sfs[-1].split('***')[0][:-2]

            # Column names for the dataframe.
            col_names = ['IX', 'IY', 'XOUT', 'YOUT', 'FEG', 'FIG', 'FF', 'FE',
                         'FVH', 'FF2', 'FE2', 'FVH2', 'FTOT1', 'FTOT2', 'TEGS',
                         'TIGS', 'CFSS', 'CFVHXS', 'VP1', 'VP2', 'FFB', 'FEB',
                         'CVHYS', 'CEYS', 'TE', 'TI', 'NE', 'VELB']

            # List to hold all the force dataframes.
            dflist = []
            for sf in lim_sfs:

                # Split up the rows for this radial location.
                foo = sf.split('\n')[1:]

                # Split up each entry in each row (each row is one big string
                # at this point).
                foo = [bar.split() for bar in foo]

                # Put into dataframe and append to list with all of them.
                df = pd.DataFrame(foo, columns=col_names)
                dflist.append(df)

            # Create a single multidimensional df out of these so it's all easier
            # to work with.
            big_df = pd.concat(dflist, keys=np.arange(1, len(dflist)))

            # Our Z data will be our force, reshape to appropriate shape.
            rows = big_df.index.max()[0]
            cols = big_df.index.max()[1]
            Z = np.array(big_df[force].values.reshape(rows, cols+1), dtype=np.float)

            # Last column is Nones. Comment out if this changes ever.
            Z = Z[:, :-1]

            # The x-values for plotting minus the nans.
            xs = np.array(big_df['XOUT'].unique(), dtype=np.float)
            xs = xs[~np.isnan(xs)]

            # The y-values.
            ys = np.array(big_df['YOUT'].unique(), dtype=np.float)
            ys = ys[~np.isnan(ys)]

            if vmin is None:
                vmin = Z.min()
            if vmax is None:
                vmax = Z.max()

            # Clip array if needed.
            Z = np.clip(Z, vmin, vmax)

            Y, X = np.meshgrid(xs, ys)

            # Plotting commands.
            fig = plt.figure()
            """
            ax = fig.add_subplot(121)
            cax = fig.add_subplot(122)
            cmap = plt.get_cmap('coolwarm')
            bounds = np.linspace(vmin, vmax, 10)
            norm = colors.BoundaryNorm(bounds, cmap.N)
            im = ax.imshow(Z, cmap=cmap, norm=norm, aspect='auto')
            cbar = fig.colorbar(im, cax=cax, cmap=cmap, norm=norm, boundaries=bounds)
            """
            ax = fig.add_subplot(111)
            divnorm = colors.DivergingNorm(vmin=vmin, vmax=vmax, vcenter=0)
            cont = ax.contourf(X, Y, Z.T, cmap='coolwarm', norm=divnorm, levels=np.linspace(vmin, vmax, 11))
            ax.set_xlim(xlim)
            cbar = fig.colorbar(cont)
            ax.set_xlabel('Parallel (m)')
            ax.set_ylabel('Radial (m)')
            cbar.ax.set_ylabel(force + ' (N)')


        def get_dep_array(self, mult_runs):

            # Only load it once. Keep track if it's already been loaded by trying
            # to see if it's been defined yet.
            try:
                self.dep_arr

            # Not defined, so create it.
            except AttributeError:

                # Create the deposition array for the initial file.
                dep_arr = np.array(self.netcdf.variables['NERODS3'][0] * -1)

                # Add the contributions from multiple runs together.
                if mult_runs:
                    file_count = 1
                    while True:
                        try:

                            # Get location to each file from user. Only .nc needed.
                            print('Choose file #{:} (press cancel to continue)'.format(file_count+1))
                            root = Tk(); root.withdraw()
                            netcdf_path = filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),))
                            if netcdf_path == '':
                                # Force an error to exit, otherwise there's a seg fault (this is a lazy fix).
                                fail
                            add_netcdf  = netCDF4.Dataset(netcdf_path)
                            add_dep     = np.array(add_netcdf.variables['NERODS3'][0] * -1)
                            dep_arr     = dep_arr + add_dep
                            file_count  = file_count + 1
                        except:
                            print('{:} files entered.'.format(file_count))
                            break

                # Define dep_arr so next time you won't have to choose all the file
                # locations.
                self.dep_arr = dep_arr

            return self.dep_arr

        def show_fig(self):
            """
            Accessor for showing the master_fig and cleaning up the layout.
            """
            self.master_fig.tight_layout()
            self.master_fig.show()
