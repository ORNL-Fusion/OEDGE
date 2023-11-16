# This script is used to pull the target Langmuir probe data and to output it in a format friendly to load into an
# Excel worksheet for further manual validation.
import sys
import argparse
import oedge_plots
import numpy as np
from scipy.signal import medfilt
import pandas as pd

# Import the get_lp script.
sys.path.append("../tools")
import get_lp

# Setup and parse arguments.
parser = argparse.ArgumentParser()
parser.add_argument("shot", help="shot number", type=int)
parser.add_argument("tmin", help="minimum time", type=float)
parser.add_argument("tmax", help="maximum time", type=float)
parser.add_argument("ncpath", help="path to the netcdf file from an OEDGE run", type=str)
parser.add_argument("-n", "--nvalues", help="number of data points needed to perform an average on", type=int,
                    default=10)
parser.add_argument("-w", "--window", help="window size of the running median filter for probe data", type=int,
                    default=51)
parser.add_argument("-o", "--outer", help="probe numbers assigned to the outer target", type=int, nargs="+")
parser.add_argument("-i", "--inner", help="probe numbers assigned to the inner target", type=int, nargs="+")
args = parser.parse_args()

# Load the Langmuir probe data.
lpdict = get_lp.get_dict_of_lps(args.shot, tunnel=False)

# Have the user specify which probe correspond to which target if it was not passed in.
if type(args.outer) == type(None):
    print("Target Specification")
    print("====================")
    print("The probe numbers and labels are printed below. You can")
    print("see where each probe is located by loading EFIT Viewer ")
    print("and turning on the Langmuir Probe diagnostic. Once you ")
    print("know which each probe is, rerun and pass them in directly")
    print("via the -i and -o flags, e.g.,")
    print("  -o 23 25 27 29 -i 51 53")
    for p, data in lpdict.items():
        print(" {:10} : {}".format(p, data["label"].strip()))

    # Exit so the user can rerun with the probe specifications passed in.
    sys.exit()

# For debugging, copy/paste this to create a sample arg object.
"""
lpdict = get_lp.get_dict_of_lps(167195, tunnel=False)
class debug_arg:

    def __init__(self):
        self.shot = 167195
        self.tmin = 4000
        self.tmax = 5000
        self.ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/oedge_tutorial/d3d-167196-osm-v1.nc"
        self.outer = [23, 25, 29, 31, 33, 35, 51, 53]
        self.inner = [131]
        self.window = 51
        self.nvalues = 10

args = debug_arg()
"""

# Load the OedgePlots object.
op = oedge_plots.OedgePlots(args.ncpath)

# For each probe...
mapped_data = {}
all_ring_data = {r: {"te_in": [], "jsat_in": [], "te_out": [], "jsat_out": []} for r in range(1, op.nrs + 1)}
for p, data in lpdict.items():
    pnum = int(p.split()[1])
    if pnum not in args.outer and pnum not in args.inner:
        continue
    elif pnum in args.outer:
        targ = "outer"
    elif pnum in args.inner:
        targ = "inner"
    print(p)

    # Create mask covering the time range of interest.
    mask = np.logical_and(data["time"] >= args.tmin, data["time"] <= args.tmax)

    # Perform a median filter to get filtered signals.
    te = medfilt(data["temp"][mask], args.window)
    jsat = medfilt(data["jsat"][mask], args.window)

    # For every data point, assign the ring number. core_or_pfz = "pfz" is used to indicate a psin of, e.g., 0.98
    # should return the ring in the PFZ and not the core.
    psin = data["psin"][mask]
    ring = np.array([op.find_ring_from_psin(p, core_or_pfz="pfz") for p in psin])

    # Store data for each ring together.
    for r in ring:
        ring_mask = ring == r
        if targ == "outer":
            all_ring_data[r]["te_out"].append(list(te[ring_mask]))
            all_ring_data[r]["jsat_out"].append(list(jsat[ring_mask]))
        elif targ == "inner":
            all_ring_data[r]["te_in"].append(list(te[ring_mask]))
            all_ring_data[r]["jsat_in"].append(list(jsat[ring_mask]))

    # Store data.
    mapped_data[p] = {"psin": psin, "ring": ring, "te": te, "jsat": jsat,
                      "r": np.full(len(psin), data["rprobe"]),
                      "z": np.full(len(psin), data["zprobe"]),
                      "label": np.full(len(psin), data["label"])}


def mean(l):
    flat_list = [item for sublist in l for item in sublist]
    return np.mean(flat_list)


def std(l):
    flat_list = [item for sublist in l for item in sublist]
    return np.std(flat_list)


# Consolidate that data for each ring where there are nvalues or more data points available.
target_data = {"outer": {}, "inner": {}}
for r in all_ring_data.keys():

    # Once if there is inner target data...
    if len(all_ring_data[r]["te_in"]) >= args.nvalues:
        target_data["inner"][r] = {"te": mean(all_ring_data[r]["te_in"]),
                                   "te_err": std(all_ring_data[r]["te_in"]),
                                   "jsat": mean(all_ring_data[r]["jsat_in"]) * 1e4,  # A/cm2 to A/m2
                                   "jsat_err": std(all_ring_data[r]["jsat_in"]) * 1e4,
                                   "target": "inner",
                                   "nvals": len(all_ring_data[r]["te_in"]),
                                   "psin": op.nc["PSIFL"][r][0]}

    # ... and again for any outer target data.
    if len(all_ring_data[r]["te_out"]) >= args.nvalues:
        target_data["outer"][r] = {"te": mean(all_ring_data[r]["te_out"]),
                                   "te_err": std(all_ring_data[r]["te_out"]),
                                   "jsat": mean(all_ring_data[r]["jsat_out"]) * 1e4,
                                   "jsat_err": std(all_ring_data[r]["jsat_out"]) * 1e4,
                                   "target": "outer",
                                   "nvals": len(all_ring_data[r]["te_out"]),
                                   "psin": op.nc["PSIFL"][r][0]}

# Combine into a single DataFrame, output to csv file.
outer_df = pd.DataFrame.from_dict(target_data["outer"], orient="index")
inner_df = pd.DataFrame.from_dict(target_data["inner"], orient="index")
df = pd.concat([outer_df, inner_df], axis=0)
fname = "{}_{}_{}.csv".format(args.shot, int(args.tmin), int(args.tmax))
df.to_csv(fname)
