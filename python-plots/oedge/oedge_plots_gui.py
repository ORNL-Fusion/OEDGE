import tkinter as tk
from tkinter import filedialog
from tkinter import font
import oedge_plots as oedge
import numpy as np


plot_opts = ['B Ratio', 'E Radial', 'E Poloidal', 'ExB Poloidal', 'ExB Radial',
             'Flow Velocity', 'Flow Velocity - Mach', 'Flow Velocity (with T13)',
             'Flow Velocity (with T13) - Mach', 'Force - Net', 'Force - FF', 'Force - FiG',
             'Force - FE', 'Force - FeG', 'Impurity Density',
             'Impurity Density - Charge', 'Impurity Ionization', 'Density',
             'Rings', 'S Coordinate', 'Electron Temperature',
             'Area of Cells', 'Ion Temperature', 'Psin']
plot_opts_cp = ['R-Rsep OMP vs. Flux - Midplane', 'R-Rsep OMP vs. Flux - Crown']
fake_opts = ['Mach', 'Te', 'ne', 'Velocity', 'nz', 'L ITF', 'L OTF']

# Spacing constants for padding.
padx = 3
pady = 4

class Window(tk.Frame):

    def __init__(self, master=None):

        # This line does something so you can like use the Tk.Frame methods or
        # something. All the examples have it at least.
        super().__init__(master)
        self.master = master
        self.master.title('OEDGE Plotting GUI')
        self.netcdf_loaded = False
        self.cp_cb_mid_var = tk.IntVar()
        self.cp_cb_top_var = tk.IntVar()
        self.cp_cb_dim_var = tk.IntVar()
        self.mr_cb_var     = tk.IntVar()
        self.core_cb_var   = tk.IntVar()
        self.charge_entry  = tk.Entry(self.master)
        self.charge_entry.insert(0, 30)
        self.vzmult_entry = tk.Entry(self.master)
        self.vzmult_entry.insert(0, 0)
        self.create_widgets()

    def add_message(self, text):
        """
        Just a small helper function to add text to the message box and scroll
        to the bottom.
        """

        self.message_box.insert(tk.END, text)
        self.message_box.see("end")

    def create_widgets(self):
        """
        Create all the buttons, message box, etc. and lay them all out. A quick
        disclaimer: I am not a pro at creating GUIs, and some of the numbers
        for like height, width, columns, etc. we landed upon through trial and
        error. This is no guarantee this looks good across different platforms.
        """

        # Variable to keep track of row number so we don't have to!
        row = 0

        # Color name to indicate background color of section.
        cname = 'gray75'
        cname2 = 'gray90'

        # Width of the columns.
        col1_width = 10
        col2_width = 30
        col3_width = 7

        # Add a message box to act as readout for errors or such.
        #self.message_box = tk.Text(self.master, height=33, width=55)
        self.message_box = tk.Text(self.master, height=7, width=65)
        self.message_box.grid(row=0, column=3, rowspan=5, padx=padx, sticky='NS', )
        self.add_message("Click 'Browse...' to load path to netCDF file.\n")

        # Add scrollbar to message box.
        self.scroll = tk.Scrollbar(self.master)
        self.scroll.grid(row=0, column=4, rowspan=5, pady=pady, sticky='NS')
        self.scroll.config(command=self.message_box.yview)
        self.message_box.config(yscrollcommand=self.scroll.set)

        # Create an Entry box for location of NetCDF file. Place this, and the associated
        # rows, into a frame that way we can distinguish sections by the color of the
        # frames they're put into.
        self.netcdf_frame = tk.Frame(self.master, bg=cname)
        self.netcdf_frame.grid(row=row, column=0, columnspan=3, sticky='WE')
        tk.Label(self.netcdf_frame, text='NetCDF File:', bg=cname, width=col1_width).grid(row=row, column=0, sticky='E')
        self.netcdf_entry = tk.Entry(self.netcdf_frame, width=col2_width)
        self.netcdf_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE')

        # Add a Browse button next to it.
        self.netcdf_button = tk.Button(self.netcdf_frame, text='Browse...', width=col3_width)
        self.netcdf_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.netcdf_button['command'] = self.browse_netcdf

        row += 1

        # Place Entry for .dat file into the same frame since this is all a part
        # of the 2D plotting.
        self.netcdf_frame.grid(row=row, column=0, columnspan=3, sticky='WE', padx=padx, pady=pady)
        tk.Label(self.netcdf_frame, text='Dat File:', bg=cname).grid(row=row, column=0, sticky='E')
        self.dat_entry = tk.Entry(self.netcdf_frame)
        self.dat_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE')

        # Add a Browse button next to it as well.
        self.dat_button = tk.Button(self.netcdf_frame, text='Browse...')
        self.dat_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.dat_button['command'] = self.browse_dat

        row += 1

        # The actual plot option. This is still in the same frame as the above NetCDF and
        # .dat file selections.
        self.netcdf_frame.grid(row=row, column=0, columnspan=3, sticky='WE', padx=padx, pady=pady)
        tk.Label(self.netcdf_frame, text='Plot:', bg=cname).grid(row=row, column=0, sticky='E')
        self.current_option = tk.StringVar(self.netcdf_frame)
        self.current_option.set(plot_opts[0])
        self.plot_option = tk.OptionMenu(self.netcdf_frame, self.current_option, *sorted(plot_opts))
        self.plot_option.grid(row=row, column=1, sticky='WE', padx=padx, pady=pady)
        self.plot_button = tk.Button(self.netcdf_frame, text='Plot')
        self.plot_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.plot_button['command'] = self.plot_command

        row += 1

        # Add button for extra plot options, still in the same frame as the NetCDF, .dat
        # and actual plot choices.
        self.extra_plot_button = tk.Button(self.netcdf_frame, text='Plot Options...')
        self.extra_plot_button.grid(row=row, column=1, columnspan=2, padx=padx, pady=pady*4)
        self.extra_plot_button['command'] = self.extra_plot_opts
        row += 1

        # Label for along ring plots. It looks fine without putting everything in a frame
        # since we alternate colors as we go down, so this section would just be the
        # same as the master background. But for the sake of consistency, we put it
        # in a frame anyways.
        self.along_frame = tk.Frame(self.master, bg=cname2)
        self.along_frame.grid(row=row, column=0, columnspan=3, sticky='WE', padx=padx, pady=pady)
        tk.Label(self.along_frame, text='Along Ring: ', width=col1_width, bg=cname2).grid(row=row, column=0, sticky='E', padx=padx, pady=pady)
        self.along_ring_entry = tk.Entry(self.along_frame, width=col2_width)
        self.along_ring_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE')
        self.plot_along_button = tk.Button(self.along_frame, text='Plot', width=col3_width)
        self.plot_along_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.plot_along_button['command'] = self.plot_along_command

        row += 1

        # Now some "fake" probe options (simulate a reciprocating probe).
        self.fake_frame = tk.Frame(self.master, bg=cname)
        self.fake_frame.grid(row=row, column=0, columnspan=3, sticky='WE', padx=padx, pady=pady)
        self.current_fake = tk.StringVar(self.fake_frame)
        self.current_fake.set(fake_opts[0])
        self.fake_option = tk.OptionMenu(self.fake_frame, self.current_fake, *fake_opts)
        self.fake_option.grid(column=0, row=0, sticky='EW', padx=padx, pady=pady)

        # Plot type is either constant R or Z, or along psin.
        self.current_fake_type = tk.StringVar(self.fake_frame)
        self.current_fake_type.set('Along R')
        self.fake_type = tk.OptionMenu(self.fake_frame, self.current_fake_type, *['Along R', 'Along Z', 'Psin'])
        self.fake_type.grid(column=0, row=1, sticky='EW', padx=padx, pady=pady)

        # A row for R coordinates, and row for Z coordinates.
        rz_entry_width = 7
        tk.Label(self.fake_frame, text='  R Start =', bg=cname).grid(row=0, column=1, sticky='WE')
        self.rstart_entry = tk.Entry(self.fake_frame, width=rz_entry_width)
        self.rstart_entry.grid(row=0, column=2, sticky='WE', padx=padx, pady=pady)
        tk.Label(self.fake_frame, text='  R End =', bg=cname).grid(row=0, column=3, sticky='WE')
        self.rend_entry = tk.Entry(self.fake_frame, width=rz_entry_width)
        self.rend_entry.grid(row=0, column=4, sticky='WE', padx=padx, pady=pady)
        tk.Label(self.fake_frame, text='  Z Start =', bg=cname).grid(row=1, column=1, sticky='WE')
        self.zstart_entry = tk.Entry(self.fake_frame, width=rz_entry_width)
        self.zstart_entry.grid(row=1, column=2, sticky='WE', padx=padx, pady=pady)
        tk.Label(self.fake_frame, text='  Z End =', bg=cname).grid(row=1, column=3, sticky='WE')
        self.zend_entry = tk.Entry(self.fake_frame, width=rz_entry_width)
        self.zend_entry.grid(row=1, column=4, sticky='WE', padx=padx, pady=pady)

        # Fill in some default values for a midplane probe.
        self.rstart_entry.insert(0, 2.22)
        self.rend_entry.insert(0, 2.37)
        self.zstart_entry.insert(0, -0.188)
        self.zend_entry.insert(0, -0.188)

        # Plot button at the end.
        self.fake_button = tk.Button(self.fake_frame, text='Plot', width=col3_width)
        self.fake_button.grid(row=0, column=5, rowspan=2, sticky='WE', padx=padx*2, pady=pady)
        self.fake_button['command'] = self.fake_plot
        row += 1

        # Entry for Thomson input file.
        self.ts_frame = tk.Frame(self.master, bg=cname2)
        self.ts_frame.grid(row=row, column=0, columnspan=3, sticky='WE', padx=padx, pady=pady)
        tk.Label(self.ts_frame, text='Thomson\n Input File: ', width=col1_width, bg=cname2).grid(row=row, column=0, sticky='E', padx=padx, pady=pady, rowspan=2)
        self.ts_entry = tk.Entry(self.ts_frame, width=col2_width)
        self.ts_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE', rowspan=2)
        self.ts_button = tk.Button(self.ts_frame, text='Browse...', width=col3_width)
        self.ts_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.ts_button['command'] = self.browse_ts
        row += 1

        # Undecided what will get here, but put a gray box as a filler.
        self.empty_frame = tk.Frame(self.master, bg=cname)
        self.empty_frame.grid(row=5, column=3, rowspan=3, sticky='NSWE', padx=padx, pady=pady)
        tk.Label(self.empty_frame, text='\n\n\n              Expansion Slot', bg=cname).grid(row=0, column=0)

        row += 1

        # Entry for Thomson file.
        tk.Label(self.ts_frame, text='Thomson\n Output File: ', bg=cname2).grid(row=row, column=0, sticky='E', padx=padx, pady=pady, rowspan=2)
        self.ts_out_entry = tk.Entry(self.ts_frame)
        self.ts_out_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE', rowspan=2)
        self.compare_button = tk.Button(self.ts_frame, text='Plot')
        self.compare_button.grid(row=row, column=2, padx=padx, pady=pady*4, sticky='WE')
        self.compare_button['command'] = self.compare_ts_command

        row += 1

        # Add a Help button.
        self.help_button = tk.Button(self.master, text='Help')
        self.help_button.grid(row=row, column=0, padx=padx, pady=pady*5, sticky='SE')
        self.help_button['command'] = self.help_command

        # Add a quit button.
        self.quit_button = tk.Button(self.master, text='Quit')
        self.quit_button.grid(row=row, column=1, padx=padx, pady=pady*5, columnspan=2, sticky='S')
        self.quit_button['command'] = self.quit_command

    def browse_netcdf(self):
        """
        Function linked to 'Browse...' button to get location of netCDF file.
        """
        self.message_box.insert(tk.END, 'Loading...\n')
        root = tk.Tk(); root.withdraw()

        # Path to DIVIMP files on Shawn's computer.
        try:
            initialdir = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/'
            netcdf_path = tk.filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),), initialdir=initialdir)
        except:
            netcdf_path = tk.filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),))
        self.netcdf_entry.delete(0, tk.END)
        self.netcdf_entry.insert(0, netcdf_path)
        self.op = oedge.OedgePlots(self.netcdf_entry.get())
        self.add_message('Loaded file: {}\n'.format(netcdf_path.split('/')[-1]))

        # Grab the .dat file while we're at it.
        dat_path = netcdf_path.split('.nc')[0] + '.dat'
        try:
            f = open(dat_path, 'r')
            self.dat_entry.delete(0, tk.END)
            self.dat_entry.insert(0, dat_path)
            self.op.add_dat_file(dat_path)
            self.add_message('Loaded file: {}\n'.format(dat_path.split('/')[-1]))
            self.dat_path = dat_path
        except:
            pass

        # Include a message saying what the first SOL ring is.
        irsep = int(self.op.nc.variables['IRSEP'][:])
        irwall = int(self.op.nc.variables['IRWALL'][:])
        self.add_message('First, last SOL ring: {}, {}\n'.format(irsep, irwall))

        # Print out how long the run took.
        try:
            time = self.op.cpu_time
            if time <= 60:
                time_str = '{} seconds'.format(time)
            elif time > 60 and time <= 3600:
                time_str = '{:.1f} minutes'.format(time / 60)
            elif time > 3600:
                time_str = '{:.1f} hours'.format(time / 3600)
            self.add_message('Run took {}.\n'.format(time_str))
        except:
            pass

        # Add a generic name for the Thomson output file.
        ts_out = netcdf_path.split('.nc')[0] + '_ts.pdf'
        self.ts_out_entry.delete(0, tk.END)
        self.ts_out_entry.insert(0, ts_out)

        # Move cursors in each entry to the end of the box so the user can see
        # the actual filename.
        self.netcdf_entry.xview("end")
        self.dat_entry.xview("end")
        self.ts_out_entry.xview("end")

    def help_command(self):
        """
        Show the help window.
        """

        self.help_window = tk.Toplevel(self.master)
        self.help_window.title('Help')

        tk.Label(self.help_window, text='OEDGE Plots GUI Help Window\n').grid(row=0, column=0, sticky='W', padx=padx, pady=pady)

        label = tk.Label(self.help_window, text='How do I use the GUI?')
        f = font.Font(label, label.cget('font'))
        f.configure(underline=True)
        label.configure(font=f)
        label.grid(row=1, column=0, sticky='W', padx=padx, pady=pady)

        tk.Label(self.help_window, text=
                  'This GUI is a means of performing some common OEDGE tasks. The first step you must always \n' +
                  'perform is supply the path to the NetCDF file that is output by an OEDGE run. The GUI will \n' +
                  'automatically search in the same directory for a .dat file, though \n' +
                  'this is not required for most tasks. The top plot command (above Plot Options...) creates\n' +
                  '2D plots, and is fairly self-explanatory. A few extra options are found in Plot Options...', justify=tk.LEFT).grid(row=2, column=0, sticky='W', padx=padx, pady=pady)

        label = tk.Label(self.help_window, text='Along a single ring plots')
        f = font.Font(label, label.cget('font'))
        f.configure(underline=True)
        label.configure(font=f)
        label.grid(row=3, column=0, sticky='W', padx=padx, pady=pady)

        tk.Label(self.help_window, text=
                 'Type in the ring number to get a target-to-target plot of whatever is selected in the above \n' +
                 'plot section (check the message window for the range of SOL rings). If you are not sure \n' +
                 'where a ring is, you can use the "Rings" plot option. Say you want ring 22. In Plot Options...,\n' +
                 'set the min and max to 21 and 23, respectively, and then click Plot in the 2D plot section.',
                 justify=tk.LEFT).grid(row=4, column=0, sticky='W', padx=padx, pady=pady)

        label = tk.Label(self.help_window, text='Thomson comparison plots')
        f = font.Font(label, label.cget('font'))
        f.configure(underline=True)
        label.configure(font=f)
        label.grid(row=5, column=0, sticky='W', padx=padx, pady=pady)

        tk.Label(self.help_window, text=
                 'The Thomson Input File is a file with Thomson scattering (TS) data created using the OMFITprofiles\n' +
                 'module via OMFIT. Instructions for creating this file can be found on the GitHub. Once the file\n' +
                 'has been created, use the Browse... button to locate it. Then clicking the Plot button just below\n' +
                 'it will generate a PDF file comparing the OSM solution to the provided experimental TS data.',
                 justify=tk.LEFT).grid(row=6, column=0, sticky='W', padx=padx, pady=pady)

        label = tk.Label(self.help_window, text='Simulate a Langmuir or Mach probe')
        f = font.Font(label, label.cget('font'))
        f.configure(underline=True)
        label.configure(font=f)
        label.grid(row=7, column=0, sticky='W', padx=padx, pady=pady)

        tk.Label(self.help_window, text=
                 'Under the message box are options to simulate a Mach probe. There are a couple other options as well. \n' +
                 'You can plot against R, Z or Psin. So if you choose along R, make sure Z start = Z end, etc.\n',
                 justify=tk.LEFT).grid(row=8, column=0, sticky='W', padx=padx, pady=pady)

    def quit_command(self):
        """
        Quit commands to avoid hanging up the terminal or anything.
        """
        self.master.quit()
        self.master.destroy()

    def extra_plot_opts(self):
        """
        Additional plot options to turn on and off via checkbuttons.
        """

        # Some constants.
        # R of tip of MiMES probe.
        rtip = 2.26
        ztip_dim = -1.10
        ztip_top = 0.93

        # Open a new window.
        self.opt_window = tk.Toplevel(self.master)
        self.opt_window.title("Plot Options")
        #self.opt_window.attributes('-topmost', 'true')

        row = 0

        # Create check boxes to turn on and off options.
        # Show MiMES probe.
        self.cp_cb_mid = tk.Checkbutton(self.opt_window, text='Show Collector Probe (MiMES) - Tip at R =', variable=self.cp_cb_mid_var)
        self.cp_cb_mid.grid(row=row, column=0, padx=padx, pady=pady, sticky='W')
        self.cp_mid_entry = tk.Entry(self.opt_window)
        self.cp_mid_entry.insert(0, rtip)
        self.cp_mid_entry.grid(row=row, column=1, padx=padx, pady=pady)

        row = row + 1

        # Show DiMES probe.
        self.cp_cb_dim = tk.Checkbutton(self.opt_window, text='Show Collector Probe (DiMES) - Tip at Z =', variable=self.cp_cb_dim_var)
        self.cp_cb_dim.grid(row=row, column=0, padx=padx, pady=pady, sticky='W')
        self.cp_dim_entry = tk.Entry(self.opt_window)
        self.cp_dim_entry.insert(0, ztip_dim)
        self.cp_dim_entry.grid(row=row, column=1, padx=padx, pady=pady)

        row = row + 1

        # Show top probe.
        self.cp_cb_top = tk.Checkbutton(self.opt_window, text='Show Collector Probe (top)      - Tip at Z =', variable=self.cp_cb_top_var)
        self.cp_cb_top.grid(row=row, column=0, padx=padx, pady=pady, sticky='W')
        self.cp_top_entry = tk.Entry(self.opt_window)
        self.cp_top_entry.insert(0, ztip_top)
        self.cp_top_entry.grid(row=row, column=1, padx=padx, pady=pady)

        row = row + 1

        # Show the metal rings.
        self.mr_cb_var = tk.IntVar()
        self.mr_cb = tk.Checkbutton(self.opt_window, text='Show Metal Rings', variable=self.mr_cb_var)
        self.mr_cb.grid(row=row, column=0, padx=padx, pady=pady, sticky='W')

        row = row + 1

        # Plot the core data or not.
        self.core_cb_var = tk.IntVar()
        self.core_cb = tk.Checkbutton(self.opt_window, text='Exclude Core Data', variable=self.core_cb_var)
        self.core_cb.grid(row=row, column=0, padx=padx, pady=pady, sticky='W')

        row = row + 1

        # Blank space
        tk.Label(self.opt_window, text=' ').grid(row=row, column=0)

        row = row + 1

        # Entry to choose what multiplier to use for FF calculation.
        tk.Label(self.opt_window, text='vz Multiplier (for FF):').grid(row=row, column=0, padx=padx, pady=pady, sticky='W')
        self.vzmult_entry = tk.Entry(self.opt_window)
        self.vzmult_entry.insert(0, 0)
        self.vzmult_entry.grid(row=row, column=1, padx=padx, pady=pady)

        row = row + 1

        # Option to choose which charge state to plot.
        tk.Label(self.opt_window, text='Charge State to Plot:').grid(row=row, column=0, padx=padx, pady=pady, sticky='W')
        self.charge_entry = tk.Entry(self.opt_window)
        self.charge_entry.insert(0, 15)
        self.charge_entry.grid(row=row, column=1, padx=padx, pady=pady)

        row = row + 1

        # Option to change vmin for plot colorbars.
        tk.Label(self.opt_window, text='Colorbar Scale Min:').grid(row=row, column=0, padx=padx, pady=pady, sticky='W')
        self.vmin_entry = tk.Entry(self.opt_window)
        self.vmin_entry.insert(0, 'auto')
        self.vmin_entry.grid(row=row, column=1, padx=padx, pady=pady)

        row = row + 1

        # Option to change vmax for plot colorbars.
        tk.Label(self.opt_window, text='Colorbar Scale Max:').grid(row=row, column=0, padx=padx, pady=pady, sticky='W')
        self.vmax_entry = tk.Entry(self.opt_window)
        self.vmax_entry.insert(0, 'auto')
        self.vmax_entry.grid(row=row, column=1, padx=padx, pady=pady)

    def browse_dat(self):
        """
        Function linked to 'Browse...' button to get location of dat file.
        """
        self.message_box.insert(tk.END, 'Loading...\n')
        root = tk.Tk(); root.withdraw()
        self.dat_path = tk.filedialog.askopenfilename(filetypes=(('Dat files', '*.dat'),))
        self.dat_entry.delete(0, tk.END)
        self.dat_entry.insert(0, self.dat_path)
        self.op.add_dat_file(self.dat_path)
        self.add_message('Loaded file: {}\n'.format(self.dat_path.split('/')[-1]))

    def browse_ts(self):
        """
        Function linked to 'Browse...' button to get location of Thomson file.
        """
        self.message_box.insert(tk.END, 'Loading...\n')
        root = tk.Tk(); root.withdraw()
        self.ts_path = tk.filedialog.askopenfilename(filetypes=(('Excel files', '*.xlsx'),))
        self.ts_entry.delete(0, tk.END)
        self.ts_entry.insert(0, self.ts_path)
        self.add_message('Loaded file: {}\n'.format(self.ts_path.split('/')[-1]))
        self.ts_entry.xview("end")

    def compare_ts_command(self):
        """
        Generate the pdf of Thomson comparisons.
        """

        self.add_message('Generating PDF...')
        ts_filename = self.ts_entry.get()
        #rings = np.append(np.arange(19, 49), np.arange(178, 190))
        #rings = np.arange(self.op.irsep, self.op.irsep + 30)  # Just do 30 rings out.
        rings = np.arange(self.op.irsep, self.op.irwall+1)
        self.op.compare_ts(ts_filename, rings, show_legend='short', output_file=self.ts_out_entry.get())
        self.add_message(' Done.\n')

    def plot_command_cp(self):
        """
        Plotting commands for the collector probe plots.
        """

        if self.current_option_cp.get() == 'R-Rsep OMP vs. Flux - Midplane':
            plot_args = {'cp_path':self.cp_path, 'xaxis':'ROMP', 'yaxis':'IMPFLUX', 'cp_num':1, 'log':True}

        elif self.current_option_cp.get() == 'R-Rsep OMP vs. Flux - Crown':
            plot_args = {'cp_path':self.cp_path, 'xaxis':'ROMP', 'yaxis':'IMPFLUX', 'cp_num':2, 'log':True}

        else:
            self.add_message('Plot option not found.')

        self.op.cp_plots(**plot_args)

    def plot_along_command(self):
        """
        Function to plot data along a specified ring.
        """

        # Just need the ring to plot along.
        ring = int(self.along_ring_entry.get())

        message = "Plotting " + self.current_option.get() + " along ring " + str(ring) + ".\n"
        self.add_message(message)

        # Can grab the relevant plot arguments from the main plot command function.
        try:
            plot_args = self.plot_command(plot_it=False)
        except AttributeError as e:
            err_msg = "Error: {}. Was DIVIMP ran, and not only OSM?\n".format(e)
            print(err_msg)
            self.add_message(err_msg)
            return None

        if 'charge' in plot_args.keys():
            charge = plot_args['charge']
        else:
            charge = None

        vz_mult = float(self.vzmult_entry.get())

        x, y = self.op.along_ring(ring, dataname=plot_args['dataname'],
                                  ylabel=plot_args['cbar_label'], charge=charge,
                                  vz_mult=vz_mult)

    def fake_plot(self):
        """
        Call the function that plots a hypothetical Mach/Langmuir probe (hence the
        name fake).
        """

        # Pick correct data to plot.
        if self.current_fake.get() == 'Te':
            plot_args = {'data':'Te'}
        elif self.current_fake.get() == 'ne':
            plot_args = {'data':'ne'}
        elif self.current_fake.get() == 'Mach':
            plot_args = {'data':'Mach'}
        elif self.current_fake.get() == 'Velocity':
            plot_args = {'data':'Velocity'}
        elif self.current_fake.get() == 'nz':
            plot_args = {'data':'nz'}
        elif self.current_fake.get() == 'L ITF':
            plot_args = {'data':'L ITF'}
        elif self.current_fake.get() == 'L OTF':
            plot_args = {'data':'L OTF'}
        else:
            self.message_box.insert(tk.END, "Error: Invalid plot argument for fake probe.\n")

        if self.current_fake_type.get() == 'Along R':
            plot_args['plot'] = 'R'
        elif self.current_fake_type.get() == 'Along Z':
            plot_args['plot'] = 'Z'
        elif self.current_fake_type.get() == 'Psin':
            plot_args['plot'] = 'psin'

        # Put in the coordinates for start and stop (R, Z).
        plot_args['r_start'] = float(self.rstart_entry.get())
        plot_args['r_end']   = float(self.rend_entry.get())
        plot_args['z_start'] = float(self.zstart_entry.get())
        plot_args['z_end']   = float(self.zend_entry.get())

        message = 'Plotting {} at constant {}.'.format(plot_args['data'], plot_args['plot'])
        self.add_message(message)
        self.op.fake_probe(**plot_args)

    def plot_command(self, plot_it=True):
        """
        For the selected DIVIMP plot, pass the correct parameters to the plotting
        function.
        """

        if self.current_option.get() == 'Electron Temperature':
            plot_args = {'dataname'  :'KTEBS',
                         'cmap'      :'inferno',
                         'cbar_label':'Te (eV)',
                         'normtype'  :'log'}

        elif self.current_option.get() == 'Ion Temperature':
            plot_args = {'dataname'  :'KTIBS',
                         'cmap'      :'inferno',
                         'cbar_label':'Ti (eV)',
                         'normtype'  :'log'}

        elif self.current_option.get() == 'Density':
            plot_args = {'dataname'  :'KNBS',
                         'cmap'      :'viridis',
                         'cbar_label':'ne (m-3)',
                         'normtype'  :'log'}

        elif self.current_option.get() == 'Temperature - Divertor':
            plot_args = {'dataname'  :'KTEBS',
                         'cmap'      :'inferno',
                         'cbar_label':'Te (eV)',
                         'normtype'  :'log',
                         'xlim'      :[1.3, 1.75],
                         'ylim'      :[-1.4, -0.9]}

        elif self.current_option.get() == 'Density - Divertor':
            plot_args = {'dataname'  :'KNBS',
                         'cmap'      :'viridis',
                         'cbar_label':'ne (m-3)',
                         'normtype'  :'log',
                         'xlim'      :[1.3, 1.75],
                         'ylim'      :[-1.4, -0.9]}

        elif self.current_option.get() == 'ExB Radial':
            plot_args = {'dataname'  :'EXB_R',
                         'cbar_label':'ExB Radial Drift (m/s?)',
                         'normtype'  :'symlog',
                         'scaling':1.0 / self.op.qtim}

        elif self.current_option.get() == 'ExB Poloidal':
            plot_args = {'dataname'  :'EXB_P',
                         'cbar_label':'ExB Poloidal Drift (m/s?)',
                         'normtype'  :'symlog',
                         'scaling':1.0 / self.op.qtim}

        elif self.current_option.get() == 'B Ratio':
            plot_args = {'dataname'  :'BRATIO',
                         'cbar_label':'B Ratio'}

        elif self.current_option.get() == 'Flow Velocity':

            # Flow velocity should be scaled by 1/QTIM.
            scaling = 1.0 / self.op.qtim
            plot_args = {'dataname'  :'KVHS',
                         'cbar_label':'Flow Velocity (m/s)',
                         'normtype'  :'symlog',
                         'scaling'   :scaling}

        elif self.current_option.get() == 'Flow Velocity - Mach':

            # Flow velocity should be scaled by 1/QTIM.
            scaling = 1.0 / self.op.qtim
            plot_args = {'dataname'  :'KVHS - Mach',
                         'cbar_label':'Mach Number',
                         'normtype'  :'symlin',
                         'scaling'   :scaling,
                         'vmin'      :-1,
                         'vmax'      :1}

        elif self.current_option.get() == 'Flow Velocity (with T13)':
            plot_args = {'dataname'  :'KVHSimp',
                         'cbar_label':'Flow Velocity (m/s)',
                         'normtype'  :'symlog'}

        elif self.current_option.get() == 'Flow Velocity (with T13) - Mach':
            plot_args = {'dataname'  :'KVHSimp - Mach',
                         'cbar_label':'Mach Number',
                         'normtype'  :'symlin',
                         'vmin'      :-1,
                         'vmax'      :1}

        elif self.current_option.get() == 'E Radial':
            plot_args = {'dataname'  :'E_RAD',
                         'cbar_label':'E Radial (V/m)',
                         'normtype'  :'symlog'}

        elif self.current_option.get() == 'E Poloidal':

            # MIGHT NEED TO BE SCALED BY e/fact.
            plot_args = {'dataname'  :'E_POL',
                         'cbar_label':'E Poloidal (V/m)',
                         'normtype'  :'symlog'}

        elif self.current_option.get() == 'Impurity Density':

            # Impurity density is scaled by ABSFAC.
            scaling = self.op.absfac
            plot_args = {'dataname'  :'DDLIMS',
                         'cbar_label':'Impurity Density (m-3)',
                         'normtype'  :'log',
                         'charge'    :'all',
                         'vmin'      :1e13,
                         'vmax'      :1e17,
                         'scaling'   :scaling,
                         'cmap'      :'nipy_spectral'}
                         #'cmap'      :'gist_ncar'}

        elif self.current_option.get() == 'Impurity Density - Charge':

            # Impurity density is scaled by ABSFAC.
            scaling = self.op.absfac
            charge = int(self.charge_entry.get())
            plot_args = {'dataname'  :'DDLIMS',
                         'cbar_label':'Impurity Density {}+ (m-3)'.format(charge),
                         'normtype'  :'log',
                         'charge'    :charge,
                         'vmin'      :1e13,
                         'vmax'      :1e17,
                         'scaling'   :scaling,
                         #'cmap'      :'nipy_spectral',
                         'cmap'      :'jet'}

        elif self.current_option.get() == 'Impurity Ionization':
            charge = int(self.charge_entry.get())
            plot_args = {'dataname'  :'KFIZS',
                         'cbar_label':'Impurity Ionization Rate {}+'.format(charge),
                         'charge'    :charge,
                         'normtype'  :'log',
                         'fix_fill'  :True}

        elif self.current_option.get() == 'S Coordinate':
            plot_args = {'dataname'  :'Snorm',
                         'cbar_label':'S Parallel (Normalized)',
                         'normtype'  :'linear',
                         'lut'       :10,
                         #'cmap'      :'nipy_spectral',
                         'cmap'      :'jet'}

        elif self.current_option.get() == 'Rings':
            plot_args = {'dataname'  :'Ring',
                         'cbar_label':'Ring'}

        elif self.current_option.get() == 'Force - FF':

            # FF does depend on charge state, so you need to pass it in here.
            charge = int(self.charge_entry.get())
            vz_mult = float(self.vzmult_entry.get())
            #plot_args = {'dataname'  :'FFI',
            #             'cbar_label':'Friction Force W{}+ (???)'.format(charge),
            #             'charge'    :charge,
            #             'normtype'  :'log'}

            plot_args = {'dataname'  : 'ff',
                         'cbar_label': 'Friction Force {}+ (N)'.format(charge),
                         'charge'    : charge,
                         'normtype'  : 'symlog',
                         'vz_mult'   : vz_mult}

        elif self.current_option.get() == 'Force - FiG':
            charge = int(self.charge_entry.get())
            plot_args = {'dataname'  : 'fig',
                         'charge'    : charge,
                         'cbar_label': 'Ion Temp. Gradient Force {}+ (N)'.format(charge),
                         'normtype'  :'symlog'}

        elif self.current_option.get() == 'Force - FE':
            charge = int(self.charge_entry.get())
            plot_args = {'dataname'  : 'fe',
                         'charge'    : charge,
                         'cbar_label': 'Electric Field Force (N)',
                         'normtype' : 'symlog'}

        elif self.current_option.get() == 'Force - FeG':
            charge = int(self.charge_entry.get())
            plot_args = {'dataname'  : 'feg',
                         'charge'    : charge,
                         'cbar_label': 'Electron Temp. Gradient Force (N)',
                         'normtype' : 'symlog'}

        elif self.current_option.get() == 'Force - Net':
            charge = int(self.charge_entry.get())
            vz_mult = float(self.vzmult_entry.get())
            plot_args = {'dataname'  : 'fnet',
                         'charge'    : charge,
                         'cbar_label': 'Net Force {}+ (N)'.format(charge),
                         'normtype'  : 'symlog',
                         'vz_mult'   : vz_mult}

        elif self.current_option.get() == 'Area of Cells':
            charge = int(self.charge_entry.get())
            plot_args = {'dataname'  : 'KAREAS',
                         'cbar_label': 'Area of Cells (m2?)'}

        elif self.current_option.get() == 'Psin':
            plot_args = {'dataname'  : 'PSIFL',
                         'cbar_label': 'Psin',
                         'vmin'      : 1.0,
                         'lut'       : 10}

        else:
            self.add_message('Plot option not found.\n')

        # Options for including collector probes or metal rings. Add into
        # dictionary if we want them.
        if self.cp_cb_mid_var.get() == 1:

            # Append the probe and tip value on if a list has already been made.
            if 'show_cp' in plot_args.keys():
                plot_args['show_cp'].append(1)
                plot_args['ptip'].append(float(self.cp_mid_entry.get()))
            else:
                plot_args['show_cp'] = [1]
                plot_args['ptip'] = [float(self.cp_mid_entry.get())]

        if self.cp_cb_top_var.get() == 1:

            # Append the probe and tip value on if a list has already been made.
            if 'show_cp' in plot_args.keys():
                plot_args['show_cp'].append(2)
                plot_args['ptip'].append(float(self.cp_top_entry.get()))
            else:
                plot_args['show_cp'] = [2]
                plot_args['ptip'] = [float(self.cp_top_entry.get())]

        if self.cp_cb_dim_var.get() == 1:

            # Append the probe and tip value on if a list has already been made.
            if 'show_cp' in plot_args.keys():
                plot_args['show_cp'].append(3)
                plot_args['ptip'].append(float(self.cp_dim_entry.get()))
            else:
                plot_args['show_cp'] = [3]
                plot_args['ptip'] = [float(self.cp_dim_entry.get())]

        # Option to show metal rings.
        if self.mr_cb_var.get() == 1:
            plot_args['show_mr'] = True

        # Option to show core data.
        if self.core_cb_var.get() == 1:
            plot_args['no_core'] = True

        # Option to supply a vmin or vmax. Put in a try since this may not even
        # be defined yet if the extra plot option window isn't open.
        try:
            if self.vmin_entry.get() != 'auto':
                plot_args['vmin'] = float(self.vmin_entry.get())
            if self.vmax_entry.get() != 'auto':
                plot_args['vmax'] = float(self.vmax_entry.get())
        except:
            pass

        if plot_it:
            message = "Plotting " + self.current_option.get() + ".\n"
            self.add_message(message)
            fig = self.op.plot_contour_polygon(**plot_args)

        return plot_args

def main():

    root = tk.Tk()
    window = Window(root)
    window.mainloop()

if __name__ == '__main__':
    main()
