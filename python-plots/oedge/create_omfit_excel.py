"""
Plug this script into an OMFIT session Command box. The steps for this workflow
are as follows.

1. Log onto iris. Type "module load omfit" and then "omfit" to bring up an
    OMFIT session.
2. Go to File > Import module > OMFITprofiles.
3. Double-click OMFITprofiles in the main OMFIT window to bring up the GUI.
4. Choose your shot and times. Times can be entered as a range, i.e.
    arange(2000, 4504, 4) to get every time between 2000-4500 in 4 ms increments
    (the time between every core TS measurement at 250 Hz).
5. In the Fetch tab select just the Thomson data, and click "Fetch and map data".
6. In the Filter tab, mess with the setting to filter out ELMs. Seems to work
    best on the core data, not so much the DTS data. Defaults seem to do alright.
    Can check plots in the Select tab.
7. Once you are happy with the filtered data, copy and paste the below script
    into the Command box back on the main OMFIT window and press Execute. A
    file called "omfit_excel_file.csv" will now be in your home directory.
8. Copy this file to your computer (WinSCP does well), and perform further
    filtering by hand. Rather simple to just delete the rows with unrealistic
    spikes in the data.
9. This file is then input into the create_ts_from_omfit function in oedge_plots.
    That function returns a file that is then used in the Thomson Input File in
    the OEDGE Plots GUI. Plan to turn this last step into a GUI only process.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from os.path import expanduser


def uncert_var(arr):
    val = np.array([a.n for a in arr])
    err = np.array([a.s for a in arr])
    return val, err

# Location of all TS data. Filtered times are in there.
ts = OMFIT['OMFITprofiles']['OUTPUTS']['RAW']['TS']

# Store all the DTS data in our DataFrame.
columns = ['subsystem', 'channel', 'time', 'psin', 'r', 'z', 'te', 'te_err', 'ne', 'ne_err']
df = pd.DataFrame(columns=columns)
for subsystem in ['core_r+1', 'core_r+0', 'divertor_r-1', 'tangential_r+0']:

    # Assign correct number of channels.
    if   subsystem == 'core_r+1':       channels = 33 + 1
    elif subsystem == 'core_r+0':       channels = 5 + 1
    elif subsystem == 'divertor_r-1':   channels = 7 + 1
    elif subsystem == 'tangential_r+0': channels = 5 + 1

    # Go through one channels at a time, appending to our main DataFrame.
    for channel in range(0, channels):
        ts_chan = ts[subsystem][channel]
        sub_df = pd.DataFrame(columns=columns)

        # Grab the Te, ne data and associated times and psin.
        te, te_err = uncert_var(ts_chan['T_e'].values.flatten())
        ne, ne_err = uncert_var(ts_chan['n_e'].values.flatten())
        time = ts_chan['time'].values.flatten()
        psin = ts_chan['psi_n'].values.flatten()
        r = ts_chan['R'].values.flatten()
        z = ts_chan['Z'].values.flatten()

        # Only include the filtered data.
        select_times = list(map(bool, ts_chan['selected_times'].values))
        te = te[select_times]; te_err = te_err[select_times]
        ne = ne[select_times]; ne_err = ne_err[select_times]
        time = time[select_times]
        psin = psin[select_times]
        subsys = np.full(len(te), subsystem)
        chan   = np.full(len(te), channel)
        r = r[select_times]; z = z[select_times]

        # Put into a DataFrame, then append onto our main DataFrame.
        sub_df['te']   = te;   sub_df['te_err'] = te_err
        sub_df['ne']   = ne;   sub_df['ne_err'] = ne_err
        sub_df['r']    = r;    sub_df['z']      = z
        sub_df['time'] = time; sub_df['psin']   = psin
        sub_df['subsystem'] = subsys
        sub_df['channel']   = chan

        df = df.append(sub_df)

# Save to a csv file (can try Excel if openpyxl is installed).
print(df)
home = expanduser("~")
#df.to_excel(home + '/omfit_excel_file.xlsx')
df.to_csv(home + '/omfit_excel_file.csv')
