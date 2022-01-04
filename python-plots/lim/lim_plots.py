"""
Author: Shawn Zamperini
Email:  zamp@utk.edu

This is a collection of plotting functions for plotting the input/output data
from a 3DLIM run.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import netCDF4
from scipy.optimize import curve_fit
from matplotlib     import colors
from matplotlib     import patches
from matplotlib     import ticker


# Some plot properties to make them a bit nicer.
plt.ion()
#plt.rcParams['font.family'] = 'serif'
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


class LimPlots:
    """
    Class object to store data and plotting routines for a 3DLIM run.
    """

    def __init__(self, ncpath, combine_repeat_runs=True):
        """
        Just intialize with the netCDF file. Will assume the lim and dat file
        are in the same folder.

        ncpath: Path to the NetCDF file.
        combine_repeat_runs: Combine the NERODS3 arrays from multiple repeat
          runs together for improved statistics. This assumes the following
          naming convention: ncfile, ncfile2, ncfile3, ...
        """

        # Save the repeat runs flags for later.
        self.combine_repeat_runs = combine_repeat_runs
        self.ncpath = ncpath

        # Just a default file for testing.
        if ncpath == 'test':
            ncpath = 'colprobe-z1-001e.nc'

        # Load in netcdf file.
        self.nc = netCDF4.Dataset(ncpath)

        # Load in the lim file.
        #limpath = ncpath[:-2] + 'lim'
        #with open(limpath) as f:
        #    self.lim = f.read()

        # load in the dat file.
        #datpath = ncpath[:-2] + 'dat'
        #with open(datpath) as f:
        #    self.dat = f.read()

        # Save file name.
        self.file = ncpath.split('/')[-1][:-3]

        # Save case name.
        self.case = ''
        for b in self.nc.variables['TITLE'][:].data:
            self.case = self.case + b.decode('UTF-8')

    def __repr__(self):

        message = 'LimPlots Object\n' + \
                  '  Case: ' + self.case + '\n' + \
                  '  File: ' + self.file + '\n'

        return message

    def summary(self):

        # Output dictionary to print out results easier.
        output = dict()

        # Case info.
        output['Case'] = self.case
        output['File'] = self.file

        # Time for run in hours.
        time = int(self.dat.split('TOTAL CPU TIME USED     (S)')[1].split('\n')[0])
        output['Time'] = str(time) + 's (' + format(time/3600, '.2f') + ' hours)'

        # Number of impurities followed.
        num = int(self.dat.split('NO OF IMPURITY IONS TO FOLLOW')[1].split('\n')[0])
        output['Ions Followed'] = "{:,}".format(num)

        # Find longest output for formatting.
        pad = 0
        for val in output.values():
            if len(str(val)) > pad:
                pad = len(str(val))

        # Printing commands.
        num_stars = 2 + 15 + 2 + pad
        print("\n" + "*" * num_stars)
        for key, val in output.items():
            print("* {:15}{:<{pad}} *".format(key, val, pad=pad))
        print("*" * num_stars)

    def get_dep_array(self, num_runs=399):
        """
        Load the deposition arrays for the collector probes.

        To-Do
        - Add option to combine the arrays of multiple repeat runs.
        """

        # Only load it once. Keep track if it's already been loaded by trying
        # to see if it's been defined yet.
        try:
            self.dep_arr

        # Not defined, so load it.
        except AttributeError:

            if self.combine_repeat_runs:

                # Try and get the dep_arr from the base case. If it doesn't exist,
                # that means nothing landed on the probe and 3DLIM won't save
                # and array of all zeros apparently. So just create the dep_arr
                # of all zeros.
                try:
                    dep_arr = np.array(self.nc.variables['NERODS3'][0] * -1)
                except:
                    dep_arr = np.zeros((6, 2*self.nc.variables['MAXNPS'][:]+1, self.nc.variables['MAXOS'][:]))
                    print("  No NERODS3.")

                # Add on contributions from repeat runs.
                for i in range(1, num_runs):
                    try:
                        ncpath_add = self.ncpath.split('.nc')[0] + str(i) + '.nc'
                        #print('Looking for {}...'.format(ncpath_add))
                        nc = netCDF4.Dataset(ncpath_add)
                        print("Found additional run: {}".format(ncpath_add))
                        try:
                            dep_arr = dep_arr + np.array(nc.variables['NERODS3'][0] * -1)
                        except KeyError:
                            print("  No NERODS3.")
                    except:
                        pass
            else:

                # Create the deposition array for the initial file.
                dep_arr = np.array(self.nc.variables['NERODS3'][0] * -1)

            # Define dep_arr so next time you won't have to choose all the file
            # locations.
            self.dep_arr = dep_arr

        return self.dep_arr


    def centerline(self, log=False, fit_exp=False, plotnum=0, show_plot=True):

        """
        Plot the ITF and OTF deposition along the centerlines on the same plot.

        log:       Option to make y axis a log scale.
        fit_exp:   Do an exponential fit onto the data and get the lambdas.

        To-Do
        - Add option so ITF/OTF is only over, say, first 5 cm.
        """

        #The deposition array.
        dep_arr = self.get_dep_array()

        # Location of each P bin, and its width.
        ps     = np.array(self.nc.variables['PS'][:].data)
        pwids  = np.array(self.nc.variables['PWIDS'][:].data)

        # Array of poloidal locations (i.e. the center of each P bin).
        pol_locs = ps - pwids/2.0

        # Distance cell centers along surface (i.e. the radial locations).
        rad_locs = np.array(self.nc.variables['ODOUTS'][:].data)

        # Get the centerline index (or closest to it).
        cline = np.abs(pol_locs).min()

        # Index the deposition array at the centerline for plotting.
        itf_x = rad_locs[np.where(rad_locs > 0.0)[0]]
        itf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
        otf_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
        otf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]

        # Plotting commands.
        if plotnum == 0:
            if show_plot:
                fig = plt.figure()
                ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum-1]

        # Option for a log axis.
        if show_plot:
            if log:
                ax.semilogy(itf_x*100, itf_y, '-', label='ITF', ms=ms, color=tableau20[6])
                ax.semilogy(otf_x*100, otf_y, '-', label='OTF', ms=ms, color=tableau20[8])
            else:
                ax.plot(itf_x*100, itf_y, '-', label='ITF', ms=ms, color=tableau20[6])
                ax.plot(otf_x*100, otf_y, '-', label='OTF', ms=ms, color=tableau20[8])

            # The tips of the probes can have spuriously high data points, so set the
            # may for the ylim just the second highest number in the datasets.
            ymax = sorted(np.concatenate((otf_y, itf_y)))[-2]

            ax.legend(fontsize=fontsize)
            ax.set_xlabel('Distance along probe (cm)', fontsize=fontsize)
            ax.set_ylabel('Deposition (arbitrary units)', fontsize=fontsize)
            #ax.set_xlim([0, 10])
            #ax.set_ylim([0, ymax])

        # Option to perform an exponential fit to the data.
        if fit_exp:
            def exp_fit(x, a, b):
                return a * np.exp(-b * x)

            popt_itf, pcov_itf = curve_fit(exp_fit, itf_x, itf_y, maxfev=5000)
            popt_otf, pcov_otf = curve_fit(exp_fit, otf_x, otf_y, maxfev=5000)

            fitx = np.linspace(0, 0.1, 100)
            fity_itf = exp_fit(fitx, *popt_itf)
            fity_otf = exp_fit(fitx, *popt_otf)

            if show_plot:
                if log:
                    ax.semilogy(fitx*100, fity_itf, '--', ms=ms, color=tableau20[6])
                    ax.semilogy(fitx*100, fity_otf, '--', ms=ms, color=tableau20[8])
                else:
                    ax.plot(fitx*100, fity_itf, '--', ms=ms, color=tableau20[6])
                    ax.plot(fitx*100, fity_otf, '--', ms=ms, color=tableau20[8])

            print("Lambdas")
            print("  ITF = {:.2f}".format(1/popt_itf[1]*100))
            print("  OTF = {:.2f}".format(1/popt_otf[1]*100))

        if plotnum ==0:
            if show_plot:
                fig.tight_layout()
                fig.show()

        print("Center ITF/OTF: {:.2f}".format(itf_y.sum()/otf_y.sum()))

        return {'itf_x':itf_x, 'itf_y':itf_y, 'otf_x':otf_x, 'otf_y':otf_y}

    def deposition_contour(self, side, probe_width=0.015, rad_cutoff=0.05,
                           plotnum=0, vmax=None, print_ratio=True):

        """
        Plot the 2D tungsten distribution across the face.

        side: Either 'ITF' or 'OTF'.
        probe_width: The half-width of the collector probe (the variable CPCO).
                     A = 0.015, B = 0.005, C = 0.0025
        rad_cutoff:  Only plot data from the tip down to rad_cutoff. Useful
                     if we want to compare to LAMS since those scans only go
                     down a certain length of the probe.

        *** To-Do ***
        - Instead of entering the width, pull out CPCO(?) from the netcdf file.
            Need to figure out the points being deposited outside the expected
            probe width first though.
        - Print out the ITF/OTF ratio from this analysis.
        """

        #The deposition array.
        dep_arr = self.get_dep_array()

        # Location of each P bin, and its width. Currently they all have the same width,
        # but it may end up such that there are custom widths so we leave it like this.
        ps     = np.array(self.nc.variables['PS'][:].data)
        pwids  = np.array(self.nc.variables['PWIDS'][:].data)

        # Array of poloidal locations (i.e. the center of each P bin).
        pol_locs = ps - pwids/2.0

        # Distance cell centers along surface (i.e. the radial locations).
        rad_locs = np.array(self.nc.variables['ODOUTS'][:].data)

        # Remove data beyond rad_cutoff.
        idx = np.where(np.abs(rad_locs)<rad_cutoff)[0]
        rad_locs = rad_locs[idx]
        dep_arr = dep_arr[:, idx]

        # Seems a junk number can sneak into PS and PWIDS. Clean that up.
        idx = np.where(ps < 9999)[0]
        pol_locs = pol_locs[idx]
        dep_arr = dep_arr[idx]

        # Get only positive values of rad_locs for ITF...
        idx = np.where(rad_locs > 0.0)[0]
        X_itf, Y_itf = np.meshgrid(rad_locs[idx], pol_locs)
        Z_itf = dep_arr[:, idx]

        # ... negative for OTF.
        idx = np.where(rad_locs < 0.0)[0]
        X_otf, Y_otf = np.meshgrid(np.abs(rad_locs[idx][::-1]), pol_locs)
        Z_otf = dep_arr[:, idx][:, ::-1]

        # Make the levels for the contour plot out of whichever side has the max deposition.
        if vmax == None:
            if Z_itf.max() > Z_otf.max():
                levels = np.linspace(0, Z_itf.max(), 15)
            else:
                levels = np.linspace(0, Z_otf.max(), 15)
        else:
            levels = np.linspace(0, vmax, 15)

        # Plotting commands.
        if side == 'ITF':
            X = X_itf; Y = Y_itf; Z = Z_itf
        else:
            X = X_otf; Y = Y_otf; Z = Z_otf

        if plotnum == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum-1]

        ax.contourf(X*100, Y*100, Z, levels=levels, cmap='Reds')
        ax.set_xlabel('Distance along probe (cm)', fontsize=fontsize)
        ax.set_ylabel('Z location (cm)', fontsize=fontsize)
        ax.set_ylim([-probe_width*100, probe_width*100])
        props = dict(facecolor='white')
        ax.text(0.75, 0.85, side, bbox=props, fontsize=fontsize*1.5, transform=ax.transAxes)

        if plotnum == 0:
            fig.tight_layout()
            fig.show()

        if print_ratio:
            print('Total ITF/OTF (0-{} cm): {:.2f}'.format(rad_cutoff*100, Z_itf.sum()/Z_otf.sum()))

    def avg_pol_profiles(self, probe_width=0.015, rad_cutoff=0.5, plotnum=0):

        """
        Plot the average poloidal profiles for each side. Mainly to see if
        deposition peaks on the edges.

        probe_width: The half-width of the collector probe (the variable CPCO).
                     A = 0.015, B = 0.005, C = 0.0025
        rad_cutoff:  Only plot data from the tip down to rad_cutoff. Useful
                     if we want to compare to LAMS since those scans only go
                     down a certain length of the probe.
        """

        # Code copied from above function, deposition_contour. See for comments.
        dep_arr = np.array(self.nc.variables['NERODS3'][0] * -1)
        ps     = np.array(self.nc.variables['PS'][:].data)
        pwids  = np.array(self.nc.variables['PWIDS'][:].data)
        pol_locs = ps - pwids/2.0
        dep_arr = dep_arr[:-1, :]
        pol_locs = pol_locs[:-1]
        rad_locs = np.array(self.nc.variables['ODOUTS'][:].data)
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

        # Plotting commands.
        if plotnum ==0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum-1]

        ax.plot(pol_locs, avg_pol_itf/avg_pol_itf.max(), label='ITF', color=tableau20[6])
        ax.plot(pol_locs, avg_pol_otf/avg_pol_otf.max(), label='OTF', color=tableau20[8])
        ax.legend(fontsize=fontsize)
        ax.set_xlabel('Poloidal (m)', fontsize=fontsize)
        ax.set_ylabel('Deposition (normalized)', fontsize=fontsize)
        ax.set_xlim([-probe_width, probe_width])

        if plotnum==0:
            fig.tight_layout()
            fig.show()

        # Print and then return the message for the GUI to use.
        message = "OTF/ITF Peaking Ratio: {:.2f}".format(otf_peak/itf_peak)
        print(message)
        #return message
        return {'pol_locs':pol_locs, 'avg_pol_itf':avg_pol_itf, 'avg_pol_otf':avg_pol_otf}

    def te_contour(self, plotnum=0):
        """
        Plot the 2D background electron plasma temperature.

        plot_num: Location in grid to place this plot. I.e. if the grid_shape
                  is (3,3), then enter a number between 0-8, where the locations
                  are labelled left to right.
        """

        # Get the connection length to restrict the plot between the two absorbing surfaces.
        cl = float(self.nc['CL'][:].data)

        # Same with the location of the plasma center (the top of the box).
        ca = float(self.nc['CA'][:].data)
        caw = float(self.nc['CAW'][:].data)

        # Get the X and Y grid data.
        x = self.nc.variables['XOUTS'][:].data
        y = self.nc.variables['YOUTS'][:].data

        # 2D grid of the temperature data.
        Z = self.nc.variables['CTEMBS'][:].data
        #Z = self.nc.variables['velplasma'][1].data

        # Trim leading and trailing zeros in the data.
        #x = x[:len(x)-np.argmin(x[::-1]==0)]
        #y = y[np.argmin(y==0):len(y)-np.argmin(y[::-1]==0)]

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
        #try:
        #    Zmin = np.partition(np.unique(Z), 1)[1]
        #    Z = np.clip(Z, Zmin, None)
        #except:
        #    pass

        # Create grid for plotting. Note we swap definitions for x and y since
        # we want the x-axis in the plot to be the parallel direction (it just
        # looks better that way).
        Y, X = np.meshgrid(x, y)

        # Plotting commands.
        if plotnum == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum-1]

        cont = ax.contourf(X, Y, Z, cmap='magma', levels=10)

        # Add a grey box to show any possible step in the boundary.
        # Swap x and y bc that's how it plotted.
        if 'xabsorb1a_step' in self.nc.variables.keys():
            y = self.nc.variables['xabsorb1a_step'][:].data
            x = self.nc.variables['yabsorb1a_step'][:].data
            rect = patches.Rectangle((x, y), 999, 999, angle=180, color='grey')
            ax.add_patch(rect)
            y = self.nc.variables['xabsorb2a_step'][:].data
            x = self.nc.variables['yabsorb2a_step'][:].data
            rect = patches.Rectangle((x, y), 999, 999, angle=270, color='grey')
            ax.add_patch(rect)

        # Set the limits. Some older runs won't have the absorbing boundaries
        # in the netcdf file, so in those cases use the connection length (which
        # isn't actually the connection length but it's something).
        if 'yabsorb1a' in self.nc.variables.keys():
            ax.set_xlim([self.nc.variables['yabsorb2a'][:], self.nc.variables['yabsorb1a'][:]])
        else:
            ax.set_xlim([-cl, cl])

        ax.set_ylim([caw, ca])
        # ax.set_ylim([None, ca])
        #ax.set_ylim([None, 0.01])  # Contour weird near edge.
        ax.set_xlabel('Parallel (m)', fontsize=fontsize)
        ax.set_ylabel('Radial (m)', fontsize=fontsize)
        if plotnum == 0:
            cbar = fig.colorbar(cont, ax=ax)
        else:
            cbar = self.master_fig.colorbar(cont, ax=ax)
        cbar.set_label('Background Te (eV)')
        if plotnum==0:
            fig.tight_layout()
            fig.show()

    def ne_contour(self, plotnum=0):
        """
        Plot the 2D background plasma density.

        plot_num: Location in grid to place this plot. I.e. if the grid_shape
                  is (3,3), then enter a number between 0-8, where the locations
                  are labelled left to right.
        """
        # Get the connection length to restrict the plot between the two absorbing surfaces.
        cl = float(self.nc['CL'][:].data)

        # Same with the location of the plasma center (the top of the box)
        ca  = float(self.nc['CA'][:].data)
        caw = float(self.nc['CAW'][:].data)

        # Get the X and Y grid data.
        x = self.nc.variables['XOUTS'][:].data
        y = self.nc.variables['YOUTS'][:].data

        # 2D grid of the temperature data.
        Z = self.nc.variables['CRNBS'][:].data

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

        if plotnum ==0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum-1]

        # Create our own levels since the automatic ones are bad.
        lev_exp = np.arange(np.floor(np.log10(Z.min()) - 1), np.ceil(np.log10(Z.max()) + 1), 0.25)
        levs = np.power(10, lev_exp)

        cont = ax.contourf(X, Y, Z, cmap='magma', levels=levs, norm=colors.LogNorm())

        # Add a grey box to show any possible step in the boundary.
        # Swap x and y bc that's how it plotted.
        if 'xabsorb1a_step' in self.nc.variables.keys():
            y = self.nc.variables['xabsorb1a_step'][:].data
            x = self.nc.variables['yabsorb1a_step'][:].data
            rect = patches.Rectangle((x, y), 999, 999, angle=180, color='grey')
            ax.add_patch(rect)
            y = self.nc.variables['xabsorb2a_step'][:].data
            x = self.nc.variables['yabsorb2a_step'][:].data
            rect = patches.Rectangle((x, y), 999, 999, angle=270, color='grey')
            ax.add_patch(rect)

        # Set the limits. Some older runs won't have the absorbing boundaries
        # in the netcdf file, so in those cases use the connection length (which
        # isn't actually the connection length but it's something).
        if 'yabsorb1a' in self.nc.variables.keys():
            ax.set_xlim([self.nc.variables['yabsorb2a'][:], self.nc.variables['yabsorb1a'][:]])
        else:
            ax.set_xlim([-cl, cl])

        ax.set_ylim([caw, ca])  # Contour weird near edge.
        ax.set_xlabel('Parallel (m)', fontsize=fontsize)
        ax.set_ylabel('Radial (m)', fontsize=fontsize)
        if plotnum == 0:
            cbar = fig.colorbar(cont,ax=ax)
        else:
            cbar = self.master_fig.colorbar(cont, ax=ax)
        cbar.set_label('Background ne (m-3)')

        if plotnum==0:
            fig.tight_layout()
            fig.show()

    def avg_imp_vely(self, plotnum=0):
        """
        SVYBAR: Average impurity velocity at X coordinates in QXS.

        plotnum: Location in grid to place this plot. I.e. if the grid_shape
                  is (3,3), then enter a number between 0-8, where the locations
                  are labelled left to right.
        """

        # Grab the data.
        x = self.nc.variables['QXS'][:].data
        y = self.nc.variables['SVYBAR'][:].data

        # Plotting commands.
        if plotnum == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum - 1]

        ax.plot(x, y, '.', ms=ms, color=tableau20[6])
        ax.set_xlabel('Radial coordinates (m)', fontsize=fontsize)
        ax.set_ylabel('Average Y imp. vel. (m/s)', fontsize=fontsize)

        if plotnum==0:
            fig.tight_layout()
            fig.show()

    def imp_contour_plot(self, plotnum=0, rmin=-0.005, rmax=0, iz_state=10, show_steps=True):

        # Get positions of the center of each bin.
        xs = self.nc.variables['XS'][:].data
        xwids = self.nc.variables['XWIDS'][:].data
        rad_locs = xs - xwids / 2.0
        ps = self.nc.variables['PS'][:].data
        pwids = self.nc.variables['PWIDS'][:].data
        pol_locs = ps - pwids / 2.0
        ys = self.nc.variables['YS'][:].data
        ywids = self.nc.variables['YWIDS'][:].data
        par_locs = ys - ywids / 2.0

        # Also mirror the par_locs to cover both sides (i.e. from -L to L instead of 0 to L).
        # Need to add a zero in the as the middle point, hence two appends.
        par_locs = np.append(np.append(-par_locs[::-1], 0), par_locs)

        # Load ddlim3 variable array of the specific ionization state.
        if type(iz_state) is list:
            # Add capability to do a range of ionization states.
            pass
        if iz_state == "all":
            ddlim3 = self.nc.variables['DDLIM3'][:].data.sum(axis=1)
        else:
            ddlim3 = self.nc.variables['DDLIM3'][:, iz_state, :, :].data

        # Sum over the radial range to create a 2D plot.
        sum_range = np.where(np.logical_and(rad_locs > rmin, rad_locs < rmax))[0]
        summed_ddlim3 = ddlim3[:, :, sum_range].sum(axis=2)

        # Plotting commands.
        X, Y = np.meshgrid(par_locs, pol_locs)
        Z = summed_ddlim3

        if plotnum == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum - 1]

        cont = ax.contourf(X, Y, Z)
        if plotnum == 0:
            cbar = fig.colorbar(cont, ax=ax)
        else:
            cbar = self.master_fig.colorbar(cont, ax=ax)
        cl = float(self.nc['CL'][:].data)
        ax.set_xlim([-cl, cl])
        ax.set_xlabel('Parallel (m)', fontsize=fontsize)
        ax.set_ylabel('Poloidal (m)', fontsize=fontsize)
        #cp = patches.Rectangle((-0.2, -0.015), width=0.4, height=0.03, color='k')
        cp = patches.Rectangle((0, 0), width=0.4, height=0.03, color='k')
        ax.add_patch(cp)
        textstr = r'Integration region:' + \
                  r'\n$\mathrm{R_min}$ = ' + str(rmin) + \
                  r'\n$\mathrm{R_max}$ = ' + str(rmax)
        props = dict(facecolor='white')

        # Show where the step in the absorbing boundary occurs.
        if show_steps:

            # Left step.
            if rmin < self.nc.variables['xabsorb1a_step'][:].data:
                #ax.axvline(self.nc.variables['yabsorb1a_step'][:].data, color='k', linestyle='--')
                rect = patches.Rectangle((self.nc.variables['yabsorb1a_step'][:].data, 999), 9999, 9999, color='grey', angle=180, alpha=0.6)
                ax.add_patch(rect)

            # Right step.
            if rmin < self.nc.variables['xabsorb2a_step'][:].data:
                #ax.axvline(self.nc.variables['yabsorb2a_step'][:].data, color='k', linestyle='--')
                rect = patches.Rectangle((self.nc.variables['yabsorb2a_step'][:].data, -999), 9999, 9999, color='grey', alpha=0.6)
                ax.add_patch(rect)

        # ax.text(0.05, 0.95, textstr, bbox=props)

        if plotnum==0:
            fig.tight_layout()
            fig.show()

    def imp_contour_plot_radial(self, plotnum=0, pmin=-0.005, pmax=0,
                                iz_state=10, show_steps=True, vmax=None):
        """
        Create a plot of the impurity density in the (parallel, radial) plane.

        plotnum:
        pmin:
        pmax:
        iz_state:
        separate_plot:
        show_steps:
        vmax:
        """

        # Get positions of the center of each bin.
        xs = self.nc.variables['XS'][:].data
        xwids = self.nc.variables['XWIDS'][:].data
        rad_locs = xs - xwids / 2.0
        ps = self.nc.variables['PS'][:].data
        pwids = self.nc.variables['PWIDS'][:].data
        pol_locs = ps - pwids / 2.0
        ys = self.nc.variables['YS'][:].data
        ywids = self.nc.variables['YWIDS'][:].data
        par_locs = ys - ywids / 2.0

        # Also mirror the par_locs to cover both sides (i.e. from -L to L instead of 0 to L).
        # Need to add a zero in the as the middle point, hence two appends.
        par_locs = np.append(np.append(-par_locs[::-1], 0), par_locs)

        # Load ddlim3 variable array of the specific ionization state.
        if type(iz_state) is list:
            # Add capability to do a range of ionization states.
            pass
        elif iz_state == "all":
            ddlim3 = self.nc.variables['DDLIM3'][:].data.sum(axis=1)
        else:
            ddlim3 = self.nc.variables['DDLIM3'][:, iz_state, :, :].data

        # This is incorrect. We need to integrate, not just sum, by multiplying by the bin widths.
        # Sum over the radial range to create a 2D plot.
        #sum_range = np.where(np.logical_and(pol_locs >= pmin, pol_locs <= pmax))[0]
        sum_range = np.where(np.logical_and(ps >= pmin, ps <= pmax))[0]

        # To truly get the range from pmin to pmax, we probably will need to reach
        # into the bin where the limit actually is.
        #                   sum_range[0]          sum_range[1]
        #                       |                     |
        #       0          1    |     2          3    |     4
        # --------------------------------------------------------
        # |          |     *****|**********|**********|******    |
        # |          |     *****|**********|**********|******    |
        # --------------------------------------------------------
        #                 ^                                 ^
        #                 |                                 |
        #               pmin                              pmax
        #
        # Since the value for ddlim3 is essentially an average for that bin, we
        # can just take abs(ps[1]-pmin) / pwidth[1] * ddlim3[1] of that bin and
        # add it on to our summed_ddlim3.

        # This isn't elegant but matrix math is confusing.
        summed_ddlim3 = np.zeros((len(par_locs), len(rad_locs)))
        for i in range(0, len(sum_range)):
            summed_ddlim3 += ddlim3[sum_range,:,:][i] * abs(ps[sum_range[i]])

        # Now add on the contributions from the edge cells it may or may not
        # spill into. If the poloidal boundary of the "left-most" ps cell does not
        # equal pmin, then that means we're spilling into a cell to the left.
        if abs(ps[sum_range[0]] - pmin) > 1e-5:
            left_cell = sum_range[0] - 1
            frac_of_left = 1 - abs(ps[left_cell] - pmin) / pwids[left_cell]
            print("ps[{}] = {:.4f} but pmin = {:.4f}. Adding on {:.2f} from left cell, ps[{}] = {:.3f}.".format(sum_range[0], ps[sum_range[0]], pmin, frac_of_left, left_cell, ps[left_cell]))
            summed_ddlim3 += ddlim3[left_cell, :, :] * frac_of_left
        if abs(ps[sum_range[-1]] - pmax) > 1e-5:
            right_cell = sum_range[-1] + 1
            frac_of_right = 1 - abs(ps[right_cell] - pmax) / pwids[right_cell]
            print("ps[{}] = {:.4f} but pmax = {:.4f}. Adding on {:.2f} from right cell, ps[{}] = {:.3f}.".format(sum_range[-1], ps[sum_range[-1]], pmax, frac_of_right, right_cell, ps[right_cell]))
            summed_ddlim3 += ddlim3[right_cell, :, :] * frac_of_right

        #summed_ddlim3 = ddlim3[sum_range, :, :].sum(axis=0)

        if vmax != None:
            summed_ddlim3 = np.clip(summed_ddlim3, 0, vmax)

        # Plotting commands.
        X, Y = np.meshgrid(par_locs, rad_locs)
        Z = summed_ddlim3

        if plotnum == 0:
            fig = plt.figure(figsize=(10,8))
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum - 1]

        def fmt(x, pos):
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)

        Z_low = sorted(np.unique(Z))[1]
        cont = ax.contourf(X, Y, Z.T, norm=colors.LogNorm(vmin=Z_low, vmax=Z.max()))
        if plotnum == 0:
            cbar = fig.colorbar(cont, ax=ax, format='%.2e')
        else:
            cbar = self.master_fig.colorbar(cont, ax=ax, format='%.2e')
        cl = float(self.nc['CL'][:].data)
        ax.set_xlim([-cl, cl])
        ax.set_ylim([rad_locs.min(), rad_locs.max()])
        ax.set_xlabel('Parallel (m)', fontsize=fontsize)
        ax.set_ylabel('Radial (m)', fontsize=fontsize)
        cbar.set_label('Imp. Density W{}+'.format(iz_state))
        # cp = patches.Rectangle((-0.2,-0.015), width=0.4, height=0.03, color='k')
        # ax.add_patch(cp)
        textstr = r'Integration region:' + '\n' + \
                  r'$\mathrm{P_{min}}$ = ' + str(pmin) + '\n' + \
                  r'$\mathrm{P_{max}}$ = ' + str(pmax)
        props = dict(facecolor='white')

        if show_steps:

            try:
                # Swap x and y bc that's how it plotted.
                y = self.nc.variables['xabsorb1a_step'][:].data
                x = self.nc.variables['yabsorb1a_step'][:].data
                rect = patches.Rectangle((x, y), 999, 999, angle=180, color='grey')
                ax.add_patch(rect)
                y = self.nc.variables['xabsorb2a_step'][:].data
                x = self.nc.variables['yabsorb2a_step'][:].data
                rect = patches.Rectangle((x, y), 999, 999, angle=270, color='grey')
                ax.add_patch(rect)
            except:
                pass

        #ax.text(0.6, 0.05, textstr, fontsize=fontsize, bbox=props, transform=ax.transAxes)
        ax.text(0.1, 0.1, textstr, fontsize=10, bbox=props, transform=ax.transAxes)

        if plotnum==0:
            fig.tight_layout()
            fig.show()

    def force_plots(self, plotnum=0, rad_loc=-0.01, cl=9.9, separate_plot=False, show_steps=True):

        # First, grab that while big force table, splitting it at the start
        # of each force table for each radial location (i.e. ix location).
        lim_sfs = self.lim.split('Static forces')[1:]

        # Seems at some point the forces stopped being recorded in the lim
        # file?
        if lim_sfs == []:
            print("Error: Forces not in .lim file. Skipped.")
            return None

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
        x  = np.array(big_df.loc[idx]['YOUT'].values,  dtype=np.float64)
        y1 = np.array(big_df.loc[idx]['FTOT1'].values, dtype=np.float64)
        y2 = np.array(big_df.loc[idx]['FTOT2'].values, dtype=np.float64)

        # Remove nans.
        x = x[~np.isnan(x)]
        y1 = y1[~np.isnan(y1)]
        y2 = y2[~np.isnan(y2)]

        # Only want values between -cl and cl.
        valid_idx = np.where(np.logical_and(x > -cl, x < cl))
        x = x[valid_idx]
        y1 = y1[valid_idx]
        y2 = y2[valid_idx]

        # If you want a separate plot made with all the forces, more detailed.
        if separate_plot:
            x = np.array(big_df.loc[idx]['YOUT'].values, dtype=np.float64)[:-1]
            valid_idx = np.where(np.logical_and(x > -cl, x < cl))
            x = x[valid_idx]

            ftot1 = np.array(big_df.loc[idx]['FTOT1'].values, dtype=np.float64)[:-1][valid_idx]
            ftot2 = np.array(big_df.loc[idx]['FTOT2'].values, dtype=np.float64)[:-1][valid_idx]
            ff1 = np.array(big_df.loc[idx]['FF'].values, dtype=np.float64)[:-1][valid_idx]
            ff2 = np.array(big_df.loc[idx]['FF2'].values, dtype=np.float64)[:-1][valid_idx]
            feg = np.array(big_df.loc[idx]['FEG'].values, dtype=np.float64)[:-1][valid_idx]
            figf = np.array(big_df.loc[idx]['FIG'].values, dtype=np.float64)[:-1][valid_idx]
            fe = np.array(big_df.loc[idx]['FE'].values, dtype=np.float64)[:-1][valid_idx]
            #fvh1 = np.array(big_df.loc[idx]['FVH'].values, dtype=np.float64)[:-1][valid_idx]
            #fvh2 = np.array(big_df.loc[idx]['FVH2'].values, dtype=np.float64)[:-1][valid_idx]

            fig = plt.figure(figsize=(7, 5))
            ax = fig.add_subplot(111)
            ax.plot(x, ftot1, '-', color=tableau20[2], label='FTOT1')
            ax.plot(x, ftot2, '--', color=tableau20[2], label='FTOT2')
            ax.plot(x, ff1, '-', color=tableau20[4], label='FF1')
            ax.plot(x, ff2, '--', color=tableau20[4], label='FF2')
            ax.plot(x, feg, '-', color=tableau20[6], label='FEG')
            ax.plot(x, figf, '-', color=tableau20[8], label='FIG')
            ax.plot(x, fe, '-', color=tableau20[10], label='FE')
            #ax.plot(x, fvh1, '-', color=tableau20[12], label='FVH1')
            #ax.plot(x, fvh2, '--', color=tableau20[12], label='FVH2')
            ax.legend(fontsize=fontsize)
            ax.set_xlabel('Parallel (m)')
            ax.set_ylabel('Force (N?)')

            # Set the limits. Some older runs won't have the absorbing boundaries
            # in the netcdf file, so in those cases use the connection length (which
            # isn't actually the connection length but it's something).
            if 'yabsorb1a' in self.nc.variables.keys():
                ax.set_xlim([self.nc.variables['yabsorb2a'][:], self.nc.variables['yabsorb1a'][:]])
            else:
                ax.set_xlim([-cl, cl])

            # Show where the step in the absorbing boundary occurs.
            if show_steps:

                # Left step.
                if rad_loc < self.nc.variables['xabsorb1a_step'][:].data:
                    #ax.axvline(self.nc.variables['yabsorb1a_step'][:].data, color='k', linestyle='--')
                    rect = patches.Rectangle((self.nc.variables['yabsorb1a_step'][:].data, 999), 9999, 9999, color='grey', angle=180, alpha=0.6)
                    ax.add_patch(rect)

                # Right step.
                if rad_loc < self.nc.variables['xabsorb2a_step'][:].data:
                    #ax.axvline(self.nc.variables['yabsorb2a_step'][:].data, color='k', linestyle='--')
                    rect = patches.Rectangle((self.nc.variables['yabsorb2a_step'][:].data, -999), 9999, 9999, color='grey', alpha=0.6)
                    ax.add_patch(rect)


            fig.tight_layout()
            fig.show()

        else:

            if plotnum == 99:
                pass
            else:
                if plotnum == 0:
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                else:
                    ax = self.master_fig.axes[plotnum - 1]
                ax.plot(x, y1, '-', color=tableau20[6], label='FTOT1')
                ax.plot(x, y2, '-', color=tableau20[8], label='FTOT2')
                ax.set_xlabel('Parallel (m)', fontsize=fontsize)
                ax.set_ylabel('Force (N?)', fontsize=fontsize)
                ax.legend(fontsize=fontsize)
                ax.axhline(0, linestyle='--', color='k')

                # Set the limits. Some older runs won't have the absorbing boundaries
                # in the netcdf file, so in those cases use the connection length (which
                # isn't actually the connection length but it's something).
                if 'yabsorb1a' in self.nc.variables.keys():
                    ax.set_xlim([self.nc.variables['yabsorb2a'][:], self.nc.variables['yabsorb1a'][:]])
                else:
                    ax.set_xlim([-cl, cl])

                # Show where the step in the absorbing boundary occurs.
                if show_steps:

                    try:
                        # Left step.
                        if rad_loc < self.nc.variables['xabsorb1a_step'][:].data:
                            #ax.axvline(self.nc.variables['yabsorb1a_step'][:].data, color='k', linestyle='--')
                            rect = patches.Rectangle((self.nc.variables['yabsorb1a_step'][:].data, 999), 9999, 9999, color='grey', angle=180, alpha=0.6)
                            ax.add_patch(rect)

                        # Right step.
                        if rad_loc < self.nc.variables['xabsorb2a_step'][:].data:
                            #ax.axvline(self.nc.variables['yabsorb2a_step'][:].data, color='k', linestyle='--')
                            rect = patches.Rectangle((self.nc.variables['yabsorb2a_step'][:].data, -999), 9999, 9999, color='grey', alpha=0.6)
                            ax.add_patch(rect)
                    except KeyError:
                        # Using an old file or just one without the step data.
                        pass

                if plotnum == 0:
                    fig.tight_layout()
                    fig.show()


    def vel_plots(self, vp, plotnum=0, cl=10.0, separate_plot=False, vmin1=None,
                  vmax1=None, vmin2=None, vmax2=None, clip=False, show_steps=True):
        """
        TODO

        vp: Either 'vp1' (background plasma) or 'vp2' (disturbed plasma from CP).
        """

        # Get the connection length to restrict the plot between the two absorbing surfaces.
        cl = float(self.nc['CL'][:].data)

        # Same with the location of the plasma center (the top of the box)
        ca  = float(self.nc['CA'][:].data)
        caw = float(self.nc['CAW'][:].data)

        # These lines are copied from the above force_plots. See them for comments.
        lim_sfs = self.lim.split('Static forces')[1:]
        if lim_sfs == []:
            print("Error: Forces not in .lim file. Skipped.")
            return None
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
        elif vp == 'vp2':
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

        if separate_plot:

            # Bounds needed for the colorbar. Will just do 10 levels.
            bounds1 = np.linspace(vmin1, vmax1, 10)
            bounds2 = np.linspace(vmin2, vmax2, 10)

            # Clip the data.
            Z1 = np.clip(Z1, vmin1, vmax1)
            Z2 = np.clip(Z2, vmin2, vmax2)

            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            cont1 = ax1.contourf(X, Y, Z1, vmin=vmin1, vmax=vmax1, cmap='coolwarm', levels=bounds1)
            cont2 = ax2.contourf(X, Y, Z2, vmin=vmin2, vmax=vmax2, cmap='coolwarm', levels=bounds2)
            cbar1 = fig.colorbar(cont1, ax=ax1, ticks=bounds1, boundaries=bounds1)
            cbar2 = fig.colorbar(cont2, ax=ax2, ticks=bounds2, boundaries=bounds2)
            ax1.set_xlabel('Parallel (m)', fontsize=fontsize)
            ax2.set_xlabel('Parallel (m)', fontsize=fontsize)
            ax1.set_ylabel('Radial (m)', fontsize=fontsize)
            ax2.set_ylabel('Radial (m)', fontsize=fontsize)
            cbar1.set_label('VP1 (m/s)')
            cbar2.set_label('VP2 (m/s)')

            if show_steps:

                try:
                    # Swap x and y bc that's how it plotted.
                    y = self.nc.variables['xabsorb1a_step'][:].data
                    x = self.nc.variables['yabsorb1a_step'][:].data
                    rect1 = patches.Rectangle((x, y), 999, 999, angle=180, color='grey')
                    rect2 = patches.Rectangle((x, y), 999, 999, angle=180, color='grey')
                    ax1.add_patch(rect1)
                    ax2.add_patch(rect2)
                    y = self.nc.variables['xabsorb2a_step'][:].data
                    x = self.nc.variables['yabsorb2a_step'][:].data
                    rect1 = patches.Rectangle((x, y), 999, 999, angle=270, color='grey')
                    rect2 = patches.Rectangle((x, y), 999, 999, angle=270, color='grey')
                    ax1.add_patch(rect1)
                    ax2.add_patch(rect2)
                except KeyError:
                    # Step info not in output file.
                    pass

            ax1.set_xlim([-cl, cl])
            ax2.set_xlim([-cl, cl])
            ax1.set_ylim([Y.min(), Y.max()])
            ax2.set_ylim([Y.min(), Y.max()])
            fig.tight_layout()
            fig.show()

        else:

            # Plotting commands.
            if plotnum == 0:
                fig = plt.figure()
                ax = fig.add_subplot(111)
            else:
                ax = self.master_fig.axes[plotnum - 1]

            cont = ax.contourf(X, Y, Z, vmin=-Z.max(), vmax=Z.max(), cmap='coolwarm')
            cbar = self.master_fig.colorbar(cont, ax=ax)
            ax.set_xlim([-cl, cl])
            ax.set_xlabel('Parallel (m)', fontsize=fontsize)
            ax.set_ylabel('Radial (m)', fontsize=fontsize)
            cbar.set_label(vp.upper() + ' (m/s)')

            # Set the limits. Some older runs won't have the absorbing boundaries
            # in the netcdf file, so in those cases use the connection length (which
            # isn't actually the connection length but it's something).
            if 'yabsorb1a' in self.nc.variables.keys():
                ax.set_xlim([self.nc.variables['yabsorb2a'][:], self.nc.variables['yabsorb1a'][:]])
            else:
                ax.set_xlim([-cl, cl])

            ax.set_ylim([caw, ca])

            if show_steps:

                try:
                    # Swap x and y bc that's how it plotted.
                    y = self.nc.variables['xabsorb1a_step'][:].data
                    x = self.nc.variables['yabsorb1a_step'][:].data
                    rect = patches.Rectangle((x, y), 999, 999, angle=180, color='grey')
                    ax.add_patch(rect)
                    y = self.nc.variables['xabsorb2a_step'][:].data
                    x = self.nc.variables['yabsorb2a_step'][:].data
                    rect = patches.Rectangle((x, y), 999, 999, angle=270, color='grey')
                    ax.add_patch(rect)
                except KeyError:
                    # Step info not in putput file.
                    pass

            if plotnum==0:
                fig.tight_layout()
                fig.show()


    def force_plot_2d(self, plotnum=0, vmin=None, vmax=None, force='FF', xlim=(-10, 10), show_steps=True):

        # First, grab that while big force table, splitting it at the start
        # of each force table for each radial location (i.e. ix location).
        lim_sfs = self.lim.split('Static forces')[1:]

        # Seems at some point the forces stopped being recorded in the lim
        # file?
        if lim_sfs == []:
            print("Error: Forces not in .lim file. Skipped.")
            return None

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
        Z = np.array(big_df[force].values.reshape(rows, cols + 1), dtype=np.float)

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

        if show_steps:

            # Swap x and y bc that's how it plotted.
            y = self.nc.variables['xabsorb1a_step'][:].data
            x = self.nc.variables['yabsorb1a_step'][:].data
            rect = patches.Rectangle((x, y), 999, 999, angle=180, color='grey')
            ax.add_patch(rect)
            y = self.nc.variables['xabsorb2a_step'][:].data
            x = self.nc.variables['yabsorb2a_step'][:].data
            rect = patches.Rectangle((x, y), 999, 999, angle=270, color='grey')
            ax.add_patch(rect)

    def overviewplot(self, iz_state=5):

        self.master_fig = plt.figure(figsize=(15,10))
        for x in range(1, 10):
            self.master_fig.add_subplot(3, 3, x)

        self.centerline(plotnum=1, log=True)
        self.deposition_contour(side='ITF', plotnum=2)
        self.deposition_contour(side='OTF', plotnum=3, print_ratio=False)
        self.te_contour(plotnum=4)
        self.ne_contour(plotnum=5)
        self.avg_pol_profiles(plotnum=6)
        #self.imp_contour_plot_radial(plotnum=7, iz_state=iz_state)
        self.force_plots(plotnum=8)
        self.vel_plots(vp='vp2', plotnum=9)

        self.master_fig.tight_layout()
        self.master_fig.show()

    def multiplot_start(self):
        """
        self.master_fig = plt.figure(figsize=(12, 6))
        for x in range(1, 10):
            self.master_fig.add_subplot(3, 3, x)
        """

    def multiplot_end(self):

        self.master_fig.tight_layout()
        self.master_fig.show()

    def compare_rcp(self, rcp_csv_path):
        """
        3DLIM automatically uses the same target data for each target. This means
        you should input just the RCP data, as it is upstream data! This script
        simply compares the upstream 3DLIM data at Y=0, which experimentally is where
        the CP and RCP data is.

        rcp_csv_path (str): Path to a csv file with the RCP data mapped to 3DLIM
          coordinate (Y=0 at the tip of the CP, negative for CP values and
          positive for beyond the CP).
        """

        rcp_df = pd.read_csv(rcp_csv_path)

        # Get the upstream data at the RCP location, assumed to be at the
        # middle of the Y data. This assumes the Y grid is symmetric about
        # Y = 0, which has always been the case here.
        mid = int(self.nc.variables['CTEMBS'][:].data.shape[0] / 2)
        te_up = self.nc.variables['CTEMBS'][:].data[mid]
        ne_up = self.nc.variables['CRNBS'][:].data[mid]
        rad = self.nc.variables['XOUTS'][:].data
        par = self.nc.variables['YOUTS'][:].data
        print("Y location of radial profile: {:.3f} m".format(par[mid]))

        # Remove zeros.
        keep = te_up != 0
        te_up = te_up[keep]
        ne_up = ne_up[keep]
        rad = rad[keep]

        # Parallel profiles at the tip as well.
        tip = np.abs(rad).argmin()
        #te_tip = self.nc.variables['CTEMBS'][:].data[:, tip]
        te_tip = self.nc.variables['velplasma'][1].data[:, tip]
        ne_tip = self.nc.variables['CRNBS'][:].data[:, tip]

        # Remove zeros.
        keep = te_tip != 0
        par = par[keep]
        te_tip = te_tip[keep]
        ne_tip = ne_tip[keep]

        # Closest ne, Te data point from the RCP.
        tip_idx = np.abs(rcp_df["X (m)"]).argmin()
        rcp_tip_te = rcp_df["Te (eV)"][tip_idx]
        rcp_tip_ne = rcp_df["ne (m-3)"][tip_idx]

        fontsize = 14

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8))

        ax1.plot(rad, te_up, color="tab:red")
        ax1.scatter(rcp_df["X (m)"], rcp_df["Te (eV)"])
        ax1.set_xlabel("X (m)", fontsize=fontsize)
        ax1.set_ylabel("Te (eV)", fontsize=fontsize)

        ax2.plot(rad, ne_up, color="tab:purple")
        ax2.scatter(rcp_df["X (m)"], rcp_df["ne (m-3)"])
        ax2.set_xlabel("X (m)", fontsize=fontsize)
        ax2.set_ylabel("ne (m-3)", fontsize=fontsize)

        ax3.plot(par, te_tip, color="tab:red")
        ax3.scatter([0], rcp_tip_te)
        ax3.set_xlabel("Y (m)", fontsize=fontsize)
        ax3.set_ylabel("Te (eV)", fontsize=fontsize)

        ax4.plot(par, ne_tip, color="tab:purple")
        ax4.scatter([0], rcp_tip_ne)
        ax4.set_xlabel("Y (m)", fontsize=fontsize)
        ax4.set_ylabel("ne (m-3)", fontsize=fontsize)

        fig.tight_layout()
        fig.show()
