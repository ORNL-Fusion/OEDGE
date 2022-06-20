"""
Author  : Shawn Zamperini
Email   : zamp@utk.edu, zamperinis@fusion.gat.com
Updated : 2/8/22

This contain a class object capable of importing, manipulating and plotting
results from a 3DLIM run. 3DLIM data is often output in a confusing format
and is not always very intuitive, so this script can help sort out the mess.

List of functions:

__init__
"""
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys


class LimPlots:

    def __init__(self, ncpath):
        """
        Create a LimPlots object by loading in the NetCDF file from the run.

        ncpath (str): Path to NetCDF file from 3DLIM run.
        """

        self.nc = netCDF4.Dataset(ncpath)

    def get_dep_array(self, additional_files=None):
        """
        Load the deposition array, typically these represent a collector probe.

        --- Inputs ---
        additional_files (list): One can pass in a list of additional netCDF
          files, from which the NERODS3 data will be added onto the current
          run. This enables running multiple identical runs instead of waiting
          for one large one to finish.

        --- Outputs ---
        dep_arr (array): The 2D array of the NERODS3 data.
        """

        # Only load it once. Keep track if it's already been loaded by trying
        # to see if it's been defined yet.
        try:
            self.dep_arr

        # Not defined, so load it.
        except AttributeError:

            # Try and get the dep_arr from the base case. If it doesn't exist,
            # that means nothing landed on the probe and 3DLIM won't save
            # an array of all zeros apparently. So just create the dep_arr
            # of all zeros.
            try:
                dep_arr = np.array(self.nc.variables['NERODS3'][0] * -1)
            except:
                dep_arr = np.zeros((6, 2*self.nc.variables['MAXNPS'][:]+1, self.nc.variables['MAXOS'][:]))
                print("  No NERODS3.")

            if type(additional_files) != type(None):

                # Add on contributions from repeat runs.
                for i in range(1, len(additional_files)):
                    try:
                        #ncpath_add = self.ncpath.split('.nc')[0] + str(i) + '.nc'
                        #print('Looking for {}...'.format(ncpath_add))
                        #nc = netCDF4.Dataset(ncpath_add)
                        nc = netCDF4.Dataset(additional_files[i])
                        print("Found additional run: {}".format(ncpath_add))
                        try:
                            dep_arr = dep_arr + np.array(nc.variables['NERODS3'][0] * -1)
                        except KeyError:
                            print("  No NERODS3.")
                    except:
                        pass
            else:

                # Create the deposition array for the initial file.
                try:
                    dep_arr = np.array(self.nc.variables['NERODS3'][0] * -1)
                except KeyError:
                    print("Error! No deposition on probe faces.")
                    return None

            # Define dep_arr so next time you won't have to choose all the file
            # locations.
            self.dep_arr = dep_arr

        return self.dep_arr

    def centerline(self, log=False, showplot=True):

        """
        Plot the deposition along the centerlines for each side of the simulated
        collector probe on the same plot. Whichever side corresponds to the
        ITF and OTF is on the user to figure out.

        --- Input ---
        log (bool): Option to make y axis a log scale.
        showplot (bool): Whether to show the plot or not.

        --- Output ---
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
        x1 = rad_locs[np.where(rad_locs > 0.0)[0]]
        y1 = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
        x2 = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
        y2 = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]

        if showplot:
            fig, ax = plt.subplots()
            ax.plot(x1, y1, label="Side 1")
            ax.plot(x2, y2, label="Side 2")
            ax.set_xlabel("Distance along probe (m)")
            ax.set_ylabel("Deposition")
            ax.legend()
            if log:
                ax.set_yscale("log")
            fig.tight_layout()
            fig.show()

        return {"x1":x1, "y1":y1, "x2":x2, "y2":y2}

    def plot_par_rad(self, dataname, pol_idx, charge="all", cbar_levels=None,
      vmin=None, vmax=None, showplot=True, show_limiter=False):
        """
        Make a 2D plot of output data in the parallel, radial plane.

        --- Input ---
        dataname (str): One of the following:
          nz: Impurity density, maps to DDLIM3
          velocity1: Background velocity used when not in a collector probe
            region, maps to first index of velplasma or the velplasma1 array
            output when running with a 2D bound.
          velocity2: As above, but for the regions within the collector probe.
        pol_idx (int): The poloidal index at which to plot the parallel, radial
          plot at.
        charge (int): For options that can be broken down into charge-specific
          results. Defaults to "all", where the results from all charges are
          summed up. 1 is the neutrals, 2 is, e.g., C1+, 3 is C2+, ... The
          reason the index is 1 greater than the charge is because in 3DLIM
          DDLIM3 goes from -1 to maxizs (e.g. maxizs=6 for C). -1 is "primary"
          neutrals, 0 is "total" neutrals, 1 is C1+, etc. When going to python
          which is 0-indexed, this means 0 is the primary neutrals, 1 is the
          total neutrals, 2 is C1+, etc.
        cbar_levels(int): Number of levels to break the colorbar into. Defaults
          to None (continuous colorscale).
        vmin, vmax (float): Sometimes the default scales for the colorbar are
          not great, so you can manually set them here.
        showplot (bool): Whether to show the plot or not. Useful if you just
          want the data.
        show_limiter (bool): Whether to show the limiter on the plot or not.
          Not implemented yet, but not really relevant for collector probe
          simulations since that is a slab limiter of zero thickness.

        --- Output ---
        Returns dictionary with the 2D data used to make the plot. X = the
          parallel coordinates, Y = the radial coordinates, Z = the data.
        """

        # Extract all relevant data up front.
        nc    = self.nc
        cl    = float(nc['CL'][:].data)
        ca    = float(nc['CA'][:].data)
        caw   = float(nc['CAW'][:].data)
        cion  = int(nc["CION"][:].data)
        xouts = nc['XOUTS'][:].data
        youts = nc['YOUTS'][:].data
        xwids = nc["XWIDS"][:].data
        yabsorb1a = float(nc["yabsorb1a"][:].data)
        yabsorb2a = float(nc["yabsorb2a"][:].data)

        # Drop bins with zero widths (arrays are made larger than necessary).
        mask = xwids != 0
        xouts = xouts[mask]
        xwids = xwids[mask]

        # The bounds will only be available if we ran with Z01 = 1.
        # To plot the boundary need to apply shift by half the bin width.
        step_y = xouts + xwids / 2
        #step_y = np.append(step_y, xouts[-1] + ca-xouts[-1])
        try:
            bounds1a = nc["bounds_1a"][:].data
            bounds2a = nc["bounds_2a"][:].data
            pos_bound = bounds1a[pol_idx]
            neg_bound = bounds2a[pol_idx]
            vary_2d_bound = True
        except:
            pos_bound = np.full(len(step_y), yabsorb1a)
            neg_bound = np.full(len(step_y), yabsorb2a)
            vary_2d_bound = False
            pass

        try:
            absfac = float(nc["ABSFAC"][:].data)
        except:
            absfac = 1.0
        print("ABSFAC = {:.2e}".format(absfac))

        # 2D representation of the data in the parallel, radial plane.
        if dataname == "nz":

            # See description for charge parameter above as to why we start
            # at index 1 for the charge.
            if charge == "all":
                print("Summing across all charge states...")
                Z = nc["DDLIM3"][:].data[pol_idx, 1:, :, :].sum(axis=0)
            else:
                Z = nc["DDLIM3"][:].data[pol_idx, charge, :, :]
            Z = Z * absfac

        # If vary_2d_bound was on (Z01 = 1), then a different velocity array
        # is used.
        elif dataname == "velocity1":
            if vary_2d_bound:
                Z = nc["velplasma_4d_1"][:, :, pol_idx].data
            else:
                Z = nc["velplasma"][0, :, :].data
        elif dataname == "velocity2":
            if vary_2d_bound:
                Z = nc["velplasma_4d_2"][:, :, pol_idx].data
            else:
                Z = nc["velplasma"][1, :, :].data

        elif dataname == "te":
            if vary_2d_bound:
                Z = nc["ctembs_3d"][:, :, pol_idx].data
            else:
                Z = nc["CTEMBS"][:].data
        elif dataname == "ti":
            if vary_2d_bound:
                Z = nc["ctembsi_3d"][:, :, pol_idx].data
            else:
                Z = nc["CTEMBSI"][:].data
        elif dataname == "ne":
            if vary_2d_bound:
                Z = nc["crnbs_3d"][:, :, pol_idx].data
            else:
                Z = nc["CRNBS"][:].data
        else:
            print("Error! Unrecognized dataname: {}".format(dataname))
            sys.exit()

        # Trim the zeros from the edges of the x and y arrays, and the
        # associated data points as well. This is done to stop this data from
        # messing up the contours in the contour plot.
        xkeep_min = np.nonzero(xouts)[0].min()
        xkeep_max = np.nonzero(xouts)[0].max()
        ykeep_min = np.nonzero(youts)[0].min()
        ykeep_max = np.nonzero(youts)[0].max()
        xouts = xouts[xkeep_min:xkeep_max]
        youts = youts[ykeep_min:ykeep_max]
        xwids = xwids[xkeep_min:xkeep_max]
        Z = Z[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

        # Furthermore, trim the data off that is beyond the abosrbing
        # boundaries.
        ykeep = np.where(np.logical_and(youts>=yabsorb2a, youts<=yabsorb1a))[0]
        youts = youts[ykeep]
        Z = Z[ykeep, :]

        # Mask zeros, if it makes sense to, and create appropriate normalization
        # for colormap.
        if dataname in ["ne", "te", "ti"]:
            Z = np.ma.masked_where(Z<=0, Z)
            if type(vmin) == type(None):
                vmin = Z.min()
            if type(vmax) == type(None):
                vmax = Z.max()
            if type(cbar_levels) == type(None):
                cmap = "inferno"
                norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            else:
                bounds = np.linspace(vmin, vmax, cbar_levels)
                cmap = plt.cm.get_cmap("inferno", cbar_levels)
                norm = mpl.colors.BoundaryNorm(boundaries=bounds,
                  ncolors=cbar_levels)
            if dataname == "te":
                cbar_label = "Te (eV)"
            elif dataname == "ti":
                cbar_label = "Ti (eV)"
            elif dataname == "ne":
                cbar_label = "ne (m-3)"

        elif dataname in ["velocity1", "velocity2"]:
            if type(vmin) == type(None):
                vmin = -np.abs(Z).max()
            if type(vmax) == type(None):
                vmax = np.abs(Z).max()
            if type(cbar_levels) == type(None):
                cmap = "coolwarm"
                norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            else:
                bounds = np.linspace(vmin, vmax, cbar_levels)
                cmap = plt.cm.get_cmap('coolwarm', cbar_levels)
                norm = mpl.colors.BoundaryNorm(boundaries=bounds,
                  ncolors=cbar_levels)
            cbar_label = "Velocity (m/s)"

        elif dataname == "nz":
            Z = np.ma.masked_where(Z<=0, Z)
            cmap = "inferno"
            if type(vmin) == type(None):
                vmin = Z[Z != 0].min() * 10000
            if type(vmax) == type(None):
                vmax = Z.max()
            norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)

            # Build the colorbar label based on the ion characteristics.
            if cion == 6:
                ion = "C"
            elif cion == 12:
                ion = "Si"
            elif cion == 74:
                ion = "W"
            else:
                ion = "?"
            if charge == "all":
                cbar_label = "{} Density (m-3)".format(ion)
            else:
                cbar_label = "{}{}+ Density (m-3)".format(ion, charge-1)

        Y, X = np.meshgrid(xouts, youts)

        # Plotting commands.
        if showplot:
            fig, ax = plt.subplots()
            cont = ax.pcolormesh(X, Y, Z, shading="auto", norm=norm, cmap=cmap)
            cbar = fig.colorbar(cont)
            ax.step(pos_bound, step_y, color="k", where="post")
            ax.step(neg_bound, step_y, color="k", where="post")
            ax.set_xlabel("Parallel to B (m)", fontsize=14)
            ax.set_ylabel("Radial (m)", fontsize=14)
            cbar.set_label(cbar_label, fontsize=14)
            fig.tight_layout()
            fig.show()

        # Return all data for whatever use.
        return {"X":X, "Y":Y, "Z":Z, "pos_bound_x":pos_bound,
          "neg_bound_x":neg_bound, "bound_y":step_y}

    def sim_reerosion(self, side, strength=1.0, launch_energy=5, charge_ion=1,
      ti_mult=1.0, delta_t=1e-8, cdf_opt=1, rad_vel=None, te_mult=1.0,
      cdf_te_cutoff=5, cdf_opt_range=[0, 2], gauss_width=10, ff_mod=1.0,
      nparts=1000, cdf_exp_falloff=0.05, nruns=1, cdf_x_cutoff=999):
        """
        This is a postprocessor that is actually a Monte Carlo simulator. It
        takes a number of different inputs, most of which are experimental
        unknowns (e.g. interpretive modelling), and applies the effects of
        reerosion post-mortem. This is subpar compared to full-scale codes
        such as GITR or WallDYN, but it's fun in a pinch.

        --- Inputs ---
        side (int): Either side 1 or 2, as labeled in the centerlien function.
        strength (float): A number from 0-1 that determines how many of the
            particles are eroded. E.g. a strength of 0.5 means only half of the
            particles in the specified region (which depends on the cdf_opt) are
            reeroded and deposited.
        launch_energy (float): The energy of the reeroded ions in eV.
        charge_ion (int): Charge of the reeroded ion. Factors into the friction
            force calculation.
        ti_mult (float): The factor for Ti = x * Te.
        delta_t (float): Time step for the simulation.
        cdf_opt (int): Determines the cumulative distribution function from
            which the locations of the rerroded particles start from:
            0: Rectangular CDF where particles are reeroded uniformly between
                the bounds specified by cdf_opt_range.
            1: The CDF is constructed from the Te data at the surface of the
                probe. The reerosion can be specified to not occur below a
                certain Te specified by cdf_te_cutoff.
            2: The CDF is an exponential with the maximum
                probability at the tip of the probe. The 1/e falloff length of
                the CDF is specified by cdf_exp_falloff.
        rad_vel (float): If None, the radial velocity will use whatever the
            convective velocity was in the 3DLIM simulation (CVPOUT). Override
            by giving a value to rad_vel.
        te_mult (float): Multiplier to the Te data at the probe surface.
        cdf_te_cutoff (float): See documentation for cdf_opt = 1.
        cdf_opt_range (list): See documentation for cdf_opt = 0.
        gauss_width (float): The radial velocity of the reeroded particles are
            chosen randomly from a Gaussian centered on rad_vel with a FWHM
            specified by this input.
        ff_mod (float): Values above 1 increase the strength of the friction
            force. This acts to decrease the radial distance at which reerode
            particles will travel, and vice-versa.
        nparts (int): The number of reerodes particles to *track*. Larger values
            just give better statistics and will NOT change how many particles
            are reeroded. You can consider this just a "how good do you want
            your statistics to be" parameter.
        cdf_exp_falloff (float): See documentation for cdf_opt = 2.
        nruns (int):
        cdf_x_cutoff (float):
        """

        from scipy.interpolate import interp1d
        import random
        from tqdm import tqdm

        # Mass of impurity in kg.
        mass_ion = float(self.nc.variables["CRMI"][:].data)
        mass_ion_kg = mass_ion * 1.66E-27

        # This little model only works if the ions radially transport via
        # convection.
        if type(rad_vel) == type(None):
            try:
                rad_vel = -float(self.nc.variables["CVPOUT"][:].data)
                print("rad_vel = {}".format(rad_vel))
                if rad_vel == 0.0:
                    print("Error! Model only works if ions are assumed to radially " + \
                      "transport via convection (arbitrary pinch velocoty in 3DLIM " + \
                      "input file).")
                    sys.exit()
            except:
                print("Warning! CVPOUT not in netCDF file. Assuming rad_vel = 350 m/s")
                rad_vel = 350

        # Pull out the deposition data. This will form the basis for or little
        # mini simulations.
        center = self.centerline(showplot=False)
        if side == 1:
            dep_x = center["x1"]
            dep_y = center["y1"]
        elif side == 2:
            dep_x = center["x2"]
            dep_y = center["y2"]
        else:
            print("Error! Side must be either 1 or 2.")
            sys.exit()
        #print("dep_y sum before: {}".format(dep_y.sum()))

        # Pull the Te data along the limiter. Assumes same data on both
        # sides.
        lim_x = self.nc.variables["QDISTS"][:].data[0][1:]
        lim_te = self.nc.variables["QTEMBS"][:].data[0][:-2] * te_mult
        lim_ne = self.nc.variables["QRNBS"][:].data[0][:-2]
        f_te = interp1d(lim_x, lim_te)
        f_ne = interp1d(lim_x, lim_ne)
        te = f_te(dep_x)
        ti = te * ti_mult
        ne = f_ne(dep_x)

        # Normalize the deposition to a relative number of particles.
        norm_fact = nparts / dep_y.sum()
        dep_y = np.array(np.round(dep_y * norm_fact), dtype=int)
        net_y = np.copy(dep_y)
        n_reerode = int(dep_y.sum() * strength)

        # The cdf_opt decides the probabilities of erosion from each location.
        # Option 0: Define range, outside of which probability = 0.
        # Option 1: CDF is proportional to Te, and is zero for Te < cdf_te_cutoff.
        if cdf_opt == 0:
            pdf = np.zeros(len(dep_x))
            pdf[np.where(np.logical_and(dep_x>=cdf_opt_range[0], dep_x<=cdf_opt_range[1]))] = 1.0
        elif cdf_opt == 1:

            # Interpolate onto probe coordinates to create unnormalized PDF.
            pdf = te
            pdf[pdf < cdf_te_cutoff] = 0.0

        # Exponentially weighted pdf towards the tip.
        elif cdf_opt == 2:
            pdf = np.exp(-dep_x / cdf_exp_falloff)

        # PDF is the deposition (essentially every ion has an equal chance of
        # reeroding, so larger chance of reerosion occuring where there are more
        # atoms).
        elif cdf_opt == 3:
            pdf = net_y

        # Anything beyond cdf_x_cutoff is considered non-reerodable

        # Create CDF.
        cdf = np.cumsum(pdf / pdf.sum())
        #if cdf[-1] != 1.0:
        #    print("Error! CDF not correct.")
        #    print(np.trapz(pdf, dep_x))
        #    print(cdf)
        #    sys.exit()

        # The parallel force will just be a constant Ffric at a single te.
        # Background deuterium velocity assumed to just be the sound speed since
        # we're simulating a region close to the probe. Speed is negative since
        # it is a velocity towards the probe, which we designate as the negative
        # direction here.
        col_log = 15
        tau_s = 1.47E13 * mass_ion * ti * np.sqrt(ti / 2.0) / \
                ((1 + 2.0 / mass_ion) * ne * np.power(charge_ion, 2) * col_log)
        vi = -np.sqrt((te + ti) * 1.609e-19 / (2.0 * 1.66e-27))

        # Just launch at same velocity every time for now.
        launch_v = np.sqrt(2 * launch_energy * 1.609e-19 / mass_ion_kg)

        # Can perform multiple runs to do multiple steps of reerosion.
        for j in range(0, nruns):
            print("Run #{}".format(j+1))

            # Arrays for stat keeping.
            x_dists = np.zeros(n_reerode)
            erodes = np.zeros(len(dep_x))

            # Track one reeroded particle at a time.
            for i in tqdm(range(0, n_reerode)):

                # Location of reeroded particle chosen at random, find nearest location in
                # array. Count as reeroded (negative particle), don't go below zero.
                erode_idx = np.argmin(np.abs(cdf - random.random()))
                erode_x = dep_x[erode_idx]
                if net_y[erode_idx] == 0:
                    continue
                else:
                    net_y[erode_idx] -= 1
                erodes[erode_idx] += 1
                start_x = erode_x

                # Let's just use DIVIMP launch option #3 because it's easy.
                # Vel/angle flag 0 : theta =+/-asin($), $ in (0,1)
                launch_ang = (-1 + random.random() * 2) * np.arcsin(random.random())
                vy_int = launch_v * np.cos(launch_ang)

                # The X velocity is chosen from a Gaussian with indicated width.
                # If negative redraw.
                vx = 0
                while vx <= 0:
                    vx = random.gauss(rad_vel, gauss_width)

                # Track particle until it redeposits (y = 0).
                x = start_x
                y = 0.001  # Start 1 mm off of surface I guess.
                vy = vy_int
                count = 0
                while True:

                    # X and Y steps are calculated up front.
                    x_step = vx * delta_t
                    y_step = vy * delta_t

                    # If y has gone below zero then it's redepositied.
                    if y <= 0:
                        redep_idx = np.argmin(np.abs(dep_x - x))
                        lost = False
                        break

                    # Update values for next iteration.
                    x += x_step
                    y += y_step
                    idx = np.argmin(np.abs(dep_x - x))
                    ff = mass_ion_kg * (vi[idx] - vy) / tau_s[idx] * ff_mod
                    vy = vy + (ff / mass_ion_kg) * delta_t
                    count += 1

                    # If the new X is less than zero, then we should consider it
                    # lost for good.
                    if x < 0:
                        lost = True
                        break

                # Accounting for where it landed.
                if not lost:
                    net_y[redep_idx] += 1
                    x_dists[i] = x - start_x

            print("Average re-eroded distance: {:.2} cm".format(x_dists.mean() * 100))

            # If we are using cdf_opt 3, we need to reassign the CDF to the new
            # deposition data.
            if cdf_opt == 3:
                pdf = net_y
                cdf = np.cumsum(pdf / pdf.sum())


        # Return after converting dep_y back to whatever units it was before.
        dep_y = dep_y / norm_fact
        net_y = net_y / norm_fact
        #print("dep_y sum after: {}".format(dep_y.sum()))
        #print("net_y sum after: {}".format(net_y.sum()))
        return {"dep_x":dep_x, "dep_y":dep_y, "net_y":net_y, "cdf":cdf}
