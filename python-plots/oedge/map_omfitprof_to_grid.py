import netCDF4
import pandas as pd
import oedge_plots
import argparse
import numpy as np
from scipy.interpolate import interp1d


# Setup and parse arguments.
parser = argparse.ArgumentParser()
parser.add_argument("tmin", help="minimum time to average across", type=float)
parser.add_argument("tmax", help="maximum time to average across", type=float)
parser.add_argument("oedge_ncpath", help="path to the netcdf file from an OEDGE run", type=str)
parser.add_argument("omfit_ncpath", help="path to the netcdf file from OMFITprofiles", type=str)
args = parser.parse_args()

# For debugging copy/paste the below.
"""
class debug_arg:

    def __init__(self):
        self.tmin = 2500
        self.tmax = 5000
        self.oedge_ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/oedge_tutorial/d3d-167196-osm-v1.nc"
        self.omfit_ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/oedge_tutorial/OMFITprofiles_167196_FIT.nc"

args = debug_arg()
"""

# Load the OMFITprofiles data and average across the designated time range.
omfit = netCDF4.Dataset(args.omfit_ncpath)
omfit_time = omfit["time"][:].data
omfit_psin = omfit["psi_n"][:].data
time_idx = np.logical_and(omfit_time >= args.tmin, omfit_time <= args.tmax)
omfit_te = omfit["T_e"][time_idx].data.mean(axis=0)
omfit_ne = omfit["n_e"][time_idx].data.mean(axis=0)

# Use interpolation functions to get the ne, Te value at each core ring.
f_te = interp1d(omfit_psin, omfit_te)
f_ne = interp1d(omfit_psin, omfit_ne)
op = oedge_plots.OedgePlots(args.oedge_ncpath)
op_te = []
op_ne = []
op_psin = []
rings = list(range(1, op.irsep))
for ring in rings:
    psin = op.nc["PSIFL"][ring-1][0]  # 0-indexed
    op_psin.append(psin)
    op_te.append(f_te(psin))
    op_ne.append(f_ne(psin))

# Put into a DataFrame, consolidate by averaging each ring value.
d = {"ring":rings, "psin":op_psin, "te":op_te, "ne":op_ne}
df = pd.DataFrame.from_dict(d).set_index("ring")
df.to_csv("omfit_mapped_to_oedge.csv")

