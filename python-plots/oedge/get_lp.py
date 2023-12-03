# This script pulls the divertor langmuir probes and puts them into a dictionary,
# among other things. It is essentially a python translation of the Matlab script
# get_lp. Just import this script and run the function get_dict_of_lps(shot).
#
# Author: Shawn Zamperini
import MDSplus as mds
import numpy   as np
import pandas  as pd
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from scipy.optimize import curve_fit
from scipy.special import erfc


def get_mds_active_probes(shot, tunnel=True):
    """
    Get the probes that were active during the shot. Used in main function.

    shot: the shot you want
    """

    # MDSplus connection to atlas where the data is store on the "LANGMUIR" tree.
    if tunnel:
        conn = mds.Connection("localhost")
    else:
        conn = mds.Connection('atlas.gat.com')
    conn.openTree("LANGMUIR", shot)

    tmin = conn.get("\LANGMUIR::TOP.TMIN").data()
    tmax = conn.get("\LANGMUIR::TOP.TMAX").data()
    try:
        runid = conn.get("\LANGMUIR::TOP.RUNID").data()
    except:
        print("Error loading RUNID")
        runid = "N/A"

    mds_index = []
    found_probes = []
    for mds_pnum in range(1,85):

        # Make sure probe name is in correct formart: 001, 002, ... , 084, 085.
        if mds_pnum < 10:
            probe = "00" + str(mds_pnum)
        else:
            probe = "0" + str(mds_pnum)

        # The complete path name to the lp. PNUM is the probe number, which does
        # not match its number in mdsplus (001 - 085).
        pname = "\LANGMUIR::TOP.PROBE_" + probe + ".PNUM"

        # Get the actual probe number if it is there. Not all MDS probes are used.
        try:
            check_pnum = conn.get(pname).data()
        except:
            pass
            #print "No data in probe " + str(probe) + "."

        # It will be '0' or blank if the MDS entry isn;t used. Otherwise it will
        # have the actual probe number in it.
        if check_pnum > 0:
            #print("Probe " + str(check_pnum) + " is MDS probe " + str(mds_pnum))
            mds_index.append(mds_pnum)
            found_probes.append(check_pnum)

    number_of_probes = len(found_probes)
    print("Found data for " + str(number_of_probes) + " probes.")

    # Store in dictionary and return it.
    active = {}
    active["tmin"] = tmin
    active["tmax"] = tmax
    active["runid"] = runid
    active["probes"] = found_probes
    active["mds_index"] = mds_index

    return active

def get_mds_lp_data(shot, mds_index, tunnel=True):
    """
    Get LP data for a single probe. Used in main function.

    shot: the shot you want
    mds_index: a number 1-85 that corresponds to the mds node. These do not
      match the probe number (which is PNUM).
    """

    # MDS connection required through atlas tunnel.
    if tunnel:
        conn = mds.Connection("localhost")
    else:
        conn = mds.Connection("atlas.gat.com")
    conn.openTree("LANGMUIR", shot)

    # Use correct form of probe name.
    if mds_index < 10:
        probe = "00" + str(mds_index)
    else:
        probe = "0" + str(mds_index)

    pname = "\LANGMUIR::TOP.PROBE_" + probe

    # All the data stored in a dictionary. All the data is in the subtree
    # indicated in pname. Just specify the node and grab the data.
    lp_data = {}
    lp_data["time"]       = conn.get(pname + ":TIME").data()
    lp_data["rprobe"]     = conn.get(pname + ":R").data()
    lp_data["zprobe"]     = conn.get(pname + ":Z").data()
    lp_data["label"]      = conn.get(pname + ":LABEL").data()
    lp_data["ntimes"]     = conn.get(pname + ":NTIMES").data()
    lp_data["pnum"]       = conn.get(pname + ":PNUM").data()
    lp_data["isat"]       = conn.get(pname + ":ISAT").data()
    lp_data["jsat"]       = conn.get(pname + ":JSAT").data()
    lp_data["temp"]       = conn.get(pname + ":TEMP").data()
    lp_data["dens"]       = conn.get(pname + ":DENS").data()
    lp_data["pot"]        = conn.get(pname + ":POT").data()
    lp_data["psin"]       = conn.get(pname + ":PSIN").data()
    lp_data["angle"]      = conn.get(pname + ":ANGLE").data()
    lp_data["area"]       = conn.get(pname + ":AREA").data()
    lp_data["delrsepout"] = conn.get(pname + ":DELRSEPOUT").data()
    lp_data["delrsepin"]  = conn.get(pname + ":DELRSEPIN").data()
    lp_data["delzsepout"] = conn.get(pname + ":DELZSEPOUT").data()
    lp_data["delzsepin"]  = conn.get(pname + ":DELZSEPIN").data()
    lp_data["csq"]        = conn.get(pname + ":CSQ").data()
    lp_data["res_err"]    = conn.get(pname + ":RES_ERR").data()
    lp_data["heatflux"]   = conn.get(pname + ":HEATFLUX").data()
    lp_data["pnum"]       = conn.get(pname + ":PNUM").data()

    # Include an estimate of the ground or SOL current.
    with np.errstate(all='ignore'):
        lp_data["ground_j"]   = lp_data["jsat"] * (1 - np.exp(-lp_data["pot"]/lp_data["temp"]))

    #print "Data stored for probe " + str(lp_data["pnum"]) + " (MDS index " + str(mds_index) + ")."

    return lp_data


def get_dict_of_lps(shot, tunnel=True):
    """
    Run this function to get the Langmuir probe data in a dictionary
    of dictionaries. Each entry will be all the probe data in the form
    of a dictionary (it just sound confusing in word it isn't really
    that weird).

    shot: the shot you want the data for.
    """

    # Get a dictionary with the probe active during this shot.
    active = get_mds_active_probes(shot, tunnel=tunnel)
    print("")

    # Get a dictionary of each probe data, then store it all in one big dictionary.
    lps = {}
    for mds_index in active["mds_index"]:
        try:
            lp_data = get_mds_lp_data(shot, mds_index, tunnel=tunnel)
            probe_name = "probe " + str(lp_data["pnum"])
            lps[probe_name] = lp_data
            print("Data stored for " + str(probe_name) + " (MDS index " + str(mds_index) + ").")
        except Exception as e:
            print("Error loading data for mds index {}".format(mds_index))
            print("Exception: {}".format(e))

    return lps

def plot_lps(shot, tmin, tmax, xtype='rminrsep', xlim=None, filter='median',
             bins=5, tunnel=True, csv_path=None, up=False, showplot=True):
    """
    Plot LP data, with optional filtering applied.

    shot (int): Shot you want data for.
    tmin (float): Start time for data.
    tmax (float): End time for data.
    xtype (str): X-axis for plot. One of "rminrsep", "psin" or "time".
    xlim ([list, float]): X limits for the plot. Entered as list or tuple,
      e.g. (-0.99, 1.4).
    filter (str): One of "median" or "average". How to treat the binned LP data.
    bins (int): Number of bins to divide the LP data up into. Divides an LP
      signal in time up into number of bins and then performs the filtering on
      each one.
    tunnel (bool): Whether to tunnel through atlas or not. If True, require ssh
      linking atlas ot localhost.
    csv_path (str): Optional path to save data to as a csv file.
    up (bool): If True only plot upper probes.
    """

    # Load lp data.
    lps = get_dict_of_lps(shot, tunnel)

    # Output lists.
    pnames_filt   = []
    x_filt        = []
    te_filt       = []
    ne_filt       = []
    jsat_filt     = []
    heatflux_filt = []
    r_filt        = []
    z_filt        = []
    ground_filt   = []
    up_probe      = []
    labels_filt   = []
    t_filt        = []

    # Go through one probe at a time to get data for plotting.
    print("\nBinning and filtering data...")
    for key in lps.keys():

        # Get times and restrict to given time range.
        times = lps[key]['time']
        idx = np.logical_and(times > tmin, times < tmax)
        times = times[idx]

        # If we get here and there's no data, continue to the next one.
        if len(times) == 0:
            print("Error! Probe {} has no data.".format(key))
            continue

        # Get Te, ne and jsat. Also get R-Rsep out.
        te       = lps[key]['temp'][idx]
        ne       = lps[key]['dens'][idx]
        jsat     = lps[key]['jsat'][idx]
        heatflux = lps[key]['heatflux'][idx]
        ground   = lps[key]['ground_j'][idx]

        if xtype == 'rminrsep':
            x = lps[key]['delrsepout'][idx]
        elif xtype == 'psin':
            x = lps[key]['psin'][idx]
        elif xtype == 'time':
            x = lps[key]['time'][idx]
        else:
            print("Error in xtype entry. Must be either rminrsep or psin.")

        # Put into one array so we can sort them all by rminrsep, low -> high.
        probe_data = np.array((x, te, ne, jsat, heatflux, ground, times))
        probe_data = probe_data[:, np.argsort(probe_data[0])]

        # If one simply wants all the data (probably messy but maybe for
        # testing reasons), you can enter a massive number and it will just
        # default to a data point per bin.
        if bins >= len(times):
            bins = len(times)

        # Divide the data up into the number of 'bins', then take the filter of each bin.
        bin_size = len(probe_data[0]) / bins
        #print("Bin size = {:.2f} --> Using integer = {}".format(bin_size, int(bin_size)))
        bin_size = int(bin_size)

        for bin in range(0, bins):
            tmp_x        = probe_data[0][bin_size*bin:bin_size*(bin+1)]
            tmp_te       = probe_data[1][bin_size*bin:bin_size*(bin+1)]
            tmp_ne       = probe_data[2][bin_size*bin:bin_size*(bin+1)]
            tmp_jsat     = probe_data[3][bin_size*bin:bin_size*(bin+1)]
            tmp_heatflux = probe_data[4][bin_size*bin:bin_size*(bin+1)]
            tmp_ground   = probe_data[5][bin_size*bin:bin_size*(bin+1)]
            tmp_t        = probe_data[6][bin_size*bin:bin_size*(bin+1)]
            #print(tmp_t)

            # Apply the preferred filter.
            if filter == 'median':
                #filter = np.median
                filter = np.nanmedian
            elif filter == 'average':
                filter = np.mean
                #filter = np.nanmean

            # Filter and add to the output lists to be plotted.
            x_filt.append(filter(tmp_x))
            te_filt.append(filter(tmp_te))
            ne_filt.append(filter(tmp_ne))
            jsat_filt.append(filter(tmp_jsat))
            heatflux_filt.append(filter(tmp_heatflux))
            ground_filt.append(filter(tmp_ground))
            #t_filt.append(tmp_t.mean())

            # Assign probe names so we can identify these data points later.
            pnames_filt.append(key)
            r_filt.append(lps[key]['rprobe'])
            z_filt.append(lps[key]['zprobe'])
            labels_filt.append(lps[key]["label"])

            # Assign if this is a lower or upper probe.
            pname = int(key[6:])
            if np.logical_and(pname>=57, pname<=112):
                up_probe.append(True)
            else:
                up_probe.append(False)

    if showplot:

        # Enumerate the pnames so they can be used for color selection.
        pnames_enum = np.array(list(enumerate(np.unique(pnames_filt))))

        # General plotting commands. These are the "Tableau 20" colors as RGB.
        tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                     (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                     (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                     (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                     (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

        # A nice looking font.
        #plt.rcParams['font.family'] = 'serif'

        # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
        for i in range(len(tableau20)):
            r, g, b = tableau20[i]
            tableau20[i] = (r / 255., g / 255., b / 255.)

        # Create figure.
        fig = plt.figure(figsize=(15,7))

        # Function for plotting.
        def plot_ax(fig, x, y, ylabel, ax_num, high_y, xlabel, xlim, legend=False, up=False):
            ax = fig.add_subplot(ax_num)

            # For each data point assign correct color.
            for i in range(0, len(x)):
                for pnames_pair in pnames_enum:
                    if pnames_pair[1] == pnames_filt[i]:
                        color = int(pnames_pair[0])
                        label = pnames_pair[1]

                if up:
                    if up_probe[i]:
                        pass
                    else:
                        continue

                ax.plot(x[i], y[i], '^', ms=10, color=tableau20[color], label=label.title())
                ax.set_xlabel(xlabel, fontsize=18)
                ax.set_ylabel(ylabel, fontsize=18)
                ax.axvline(0.0, linestyle='--', color='k')
                ax.set_xlim(xlim)
                ax.set_ylim([0, high_y])

                # Process to remove duplicate legend entries.
                if legend:
                    handles, labels = plt.gca().get_legend_handles_labels()
                    newLabels, newHandles = [], []
                    for handle, label in zip(handles, labels):
                      if label not in newLabels:
                        newLabels.append(label)
                        newHandles.append(handle)
                    ax.legend(newHandles, newLabels, framealpha=0.5)

        # Assign plot limits, if specified.
        if xtype == 'rminrsep':
            xlabel = 'R-Rsep (m)'
            if xlim is None:
                xlim = [-0.05, 0.1]
        elif xtype == 'psin':
            xlabel = 'Psin'
            if xlim is None:
                xlim = [0.98, 1.1]
        elif xtype == 'time':
            xlabel = 'Time (ms)'

        plot_ax(fig, x_filt, te_filt, 'Te (eV)', 131, 50, xlabel, xlim, legend=True, up=up)
        plot_ax(fig, x_filt, ne_filt, 'ne (cm-3)', 132, 10e13, xlabel, xlim, up=up)
        plot_ax(fig, x_filt, jsat_filt, 'jsat (A/cm2)', 133, 100, xlabel, xlim, up=up)
        fig.tight_layout()
        fig.show()

    # Organize into dictionary for output.
    lp_dict = {xtype:x_filt, 'Te (eV)':te_filt, 'ne (cm-3)':ne_filt,
           'jsat (A/cm2)':jsat_filt, 'heatflux (W/cm2)':heatflux_filt,
           'ground_filt':ground_filt, 'pnames':pnames_filt, 'R':r_filt,
           'Z':z_filt, 'labels':labels_filt}

    # Output to a csv file.
    if csv_path != None:
        df = pd.DataFrame(lp_dict)
        df.to_csv(csv_path)

    return lp_dict

def plot_contours(lp_dict, tstart, tend, region="sas", filter_type="median",
    filter_window=201, te_vmin=None, te_vmax=None, show_locs=False):
    """
    Plot contours of the profiles with respect to the separatrix.
    """

    # Based on region select the right probes.
    data = {}
    if region == "sas":
        for pname in lp_dict.keys():
            label = lp_dict[pname]["label"].strip()
            print(label)
            if label in ["A15", "A16", "A17", "A18", "A19", "A20", "A21", "A22"]:
                continue
            elif label[0] == "A":
                data[pname] = lp_dict[pname]

    if filter_type == "median":
        from scipy.signal import medfilt
        filter = medfilt

    # Construct full arrays consisting of the coordinates time, psin and the
    # measurement.
    time = []
    psin = []
    te = []
    jsat = []
    labels = []
    pnames = []
    te_raw = []
    jsat_raw = []
    print("Filtering data...")
    for pname in data.keys():

        t = data[pname]["time"]
        keep = np.where(np.logical_and(t>=tstart, t<=tend))[0]
        time = np.append(time, data[pname]["time"][keep])
        psin = np.append(psin, data[pname]["psin"][keep])
        te = np.append(te, filter(data[pname]["temp"][keep], filter_window))
        te_raw = np.append(te_raw, data[pname]["temp"][keep])
        jsat = np.append(jsat, filter(data[pname]["jsat"][keep], filter_window))
        jsat_raw = np.append(jsat_raw, data[pname]["temp"][keep])
        labels = np.append(labels, np.full(len(keep), data[pname]["label"].strip()))
        pnames = np.append(pnames, np.full(len(keep), pname))
    time = np.array(time).flatten()
    psin = np.array(psin).flatten()
    te = np.array(te).flatten()
    jsat = np.array(jsat).flatten()
    labels = labels.flatten()
    pnames = pnames.flatten()

    if type(te_vmin) == type(None):
        te_vmin = te.min()
    if type(te_vmax) == type(None):
        te_vmax = te.max()
    te = np.clip(te, te_vmin, te_vmax)
    #print(te_vmin)
    #print(te_vmax)


    print("Creating triangle contour plots...")
    bounds = np.linspace(te_vmin, te_vmax, 20)
    norm = BoundaryNorm(bounds, len(bounds))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, sharey=True)
    ax1.set_facecolor("grey")
    ax2.set_facecolor("grey")

    cont1 = ax1.tricontourf(time, psin, te, cmap="hot", vmin=te_vmin, vmax=te_vmax)
    ax1.tricontour(time, psin, te, levels=[10], colors="k", linestyles="--", linewidths=2)
    ax1.axhline(1.0, color="k", lw=3)
    ax1.set_ylabel("Psin")
    cbar1 = fig.colorbar(cont1, ax=ax1, norm=norm, boundaries=bounds)
    cbar1.set_label("Te (eV)")
    if show_locs:
        ax1.plot(time, psin, color="k")

    cont2 = ax2.tricontourf(time, psin, jsat, cmap="winter")
    ax2.axhline(1.0, color="k", lw=3)
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylabel("Psin")
    cbar2 = fig.colorbar(cont2, ax=ax2)
    cbar2.set_label("jsat (A/cm2)")

    fig.tight_layout()
    fig.show()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    for pname in np.unique(pnames):
        mask = pnames == pname
        ax1.plot(time[mask], te[mask], label=labels[mask][0])
        ax2.plot(time[mask], jsat[mask], label=labels[mask][0])

    ax1.set_xlabel("Time (ms)")
    ax1.set_ylabel("Te (eV)")
    ax2.set_ylabel("jsat (A/cm2)")
    ax1.legend()

    fig.tight_layout()
    fig.show()

    return {"time":time, "psin":psin, "te":te, "jsat":jsat, "pnames":pnames,
        "labels":labels, "te_raw":te_raw, "jsat_raw":jsat_raw}



def fit_conv_gauss(lp_xl_path, lp_xl_sheet="Data Fixed", lp_xl_xdata="psin",
    lp_xl_ydata="jsat fixed (A/cm2)", gauss_range=[1.0, 1.04], ylabel=None,
    skiprows=0):
    """
    To get to the point where you would want to use this function probably
    requires a little manual labor. Supply here an Excel file where you have
    already manipulated the data to make it all agree (perhaps you needed to
    multiply some probes by constants to bring them into agreement with the
    other probes for example).

    Is this function flawproof? Absolutely not! But hopefully you can mess with
    it enough to get a good fit, which you can afterwards do one more round of
    manual tickering on where the exponential and gaussian fits don't overlap.

    lp_xl_path (str): Path to Excel file where you have manually organized
      the LP data to make sure it looks good for fitting.
    lp_xl_sheet (str): The sheet in the Excel file with your data you want to
      fit.
    lp_xl_ydata (str): The column name of either your jsat or Te data to be
      fitted in the Excel file.
    gauss_range (list, float): Between these psin values fit to a convoluted
      gaussian. Outside, fit to exponentials.
    ylabel (str): ylabel for the plot.
    skiprows (int): For the pd.read_excel function.
    """

    # The fitting functions.
    def exp_fit_left(x, a, b, c):
        return a * np.exp(b * (x - c))

    def exp_fit_right(x, a, b, c):
        return a * np.exp(-b * (x - c))

    def gauss_conv_exp_fit(s, width, lambda_n, n0, n_bg, s0):
        fx = 5
        return n0 / 2.0 * np.exp((width/(2*lambda_n*fx))**2 - (s-s0)/(lambda_n *
          fx)) * erfc(width/(2*lambda_n*fx) - (s-s0)/width)

    # Load Excel file and pull out the data for the fitting.
    df = pd.read_excel(lp_xl_path, sheet_name=lp_xl_sheet, skiprows=skiprows)
    try:
        psin = df[lp_xl_xdata].values
        y = df[lp_xl_ydata].values
    except KeyError as e:
        print(e)
        print("Options are: {}".format(df.columns))

    # Sort the data.
    sidx = np.argsort(psin)
    psin = psin[sidx]
    y    = y[sidx]

    # Get rid of nans if some slip in.
    kidx = ~np.isnan(y)
    psin = psin[kidx]
    y    = y[kidx]

    # Divide the data up into our three regions, i.e. the region with the
    # guassian fit (center) and the two exponential fits to the "left" and
    # "right" that surround it.
    left   = np.where(psin < gauss_range[0])[0]
    right  = np.where(psin > gauss_range[1])[0]
    center = np.where(np.logical_and(psin >= gauss_range[0], psin <= gauss_range[1]))[0]

    # Do the exponential fits first.
    if len(psin[left]) > 0:
        left_popt, left_pcov = curve_fit(exp_fit_left, psin[left], y[left],
          p0=(1, 10, 1), maxfev=5000)

    # If you want the gaussian fit to just go all the way and don't include any
    # data beyond it for an exponential.
    if len(psin[right]) > 0:
        right_popt, right_pcov = curve_fit(exp_fit_right, psin[right], y[right],
          p0=(1, 10, 1), maxfev=5000)

    # Now the convoluted gaussian fit.
    guess = (0.05, 0.02, y.max(), 0.0, 1.0)
    center_popt, center_pcov = curve_fit(gauss_conv_exp_fit, psin[center],
      y[center], p0=guess, maxfev=5000)

    # Plotting to see how it all looks.
    lw = 4
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(psin, y, "k.")

    if len(psin[left]) > 0:
        ax.plot(psin[left], exp_fit_left(psin[left], *left_popt), "-", color="springgreen", lw=lw)

    if len(psin[right]) > 0:
        ax.plot(psin[right], exp_fit_right(psin[right], *right_popt), "-", color="deepskyblue", lw=lw)

    ax.plot(psin[center], gauss_conv_exp_fit(psin[center], *center_popt), "-", color="crimson", lw=lw)

    ax.set_xlabel("Psin", fontsize=16)
    if ylabel == None:
        ylabel = lp_xl_ydata
    ax.set_ylabel(ylabel, fontsize=16)
    ax.grid()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.show()

    # Just return the results of the fit with points 0.0005 psin a.
    #psin_fit = np.arange(psin.min(), psin.max(), 0.0005)
    psin_fit = np.linspace(psin.min(), psin.max(), 100)
    y_fit = np.zeros(len(psin_fit))

    center_idx = np.where(np.logical_and(psin_fit > psin[center][0], psin_fit < psin[center][-1]))[0]
    y_fit[center_idx] = gauss_conv_exp_fit(psin_fit[center_idx], *center_popt)

    if len(psin[left]) > 0:
        left_idx = np.where(psin_fit < psin[left].max())[0]
        y_fit[left_idx] = exp_fit_left(psin_fit[left_idx], *left_popt)

        # The points in between the fits just put in a linear fit between the
        # last two points.
        l = left_idx[-1]
        c = center_idx[0] + 1
        m = (y_fit[c] - y_fit[l]) / (psin_fit[c] - psin_fit[l])
        y_fit[l:c] = m * (psin_fit[l:c] - psin_fit[l]) + y_fit[l]

    if len(psin[right]) > 0:
        right_idx = np.where(psin_fit > psin[right].min())[0]
        y_fit[right_idx] = exp_fit_right(psin_fit[right_idx], *right_popt)
        r = right_idx[0] + 1
        c = center_idx[-1]
        m = (y_fit[r] - y_fit[c]) / (psin_fit[r] - psin_fit[c])
        y_fit[c:r] = m * (psin_fit[c:r] - psin_fit[c]) + y_fit[c]

    return {"psin_fit":psin_fit, "y_fit":y_fit, "psin":psin, "y":y}
