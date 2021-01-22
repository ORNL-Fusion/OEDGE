'''

Author: Kenneth Bott
Email: kbott@vols.utk.edu

This file creates a GUI to make graphing of fusion data easier

'''

import tkinter as tk
from tkinter import filedialog
from tkinter import font
from PIL     import ImageTk, Image
import numpy as np
import dlim_plots as dlim
import lim_plots as limpt
import matplotlib.pyplot as plt

# Color name to indicate background color of section.
cname = 'gray75'
cname2 = 'gray90'

#spacing constants
padx = 3
pady = 4

#Width of columns
col1_width = 10
col2_width = 10
col3_width = 5

#plot options
plot_op = ['Please Select Option', 'Center Line', 'Contour', 'Poloidal Profiles',
           'Temperature Contour', 'Plasma Density', 'Impurity Velocity',
           'Impurity Contour']

class Window(tk.Frame):

    def __init__(self, master=None):

        # This line does something so you can like use the Tk.Frame methods or
        # something. All the examples have it at least.
        super().__init__(master)
        self.master = master
        self.master.title('3DLIM Plotting GUI')
        self.master.rowconfigure(0, weight=1)
        self.master.columnconfigure(0, weight=1)
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

    def create_widgets(self):
        """
        Create all the buttons, message box, etc. and lay them all out.
        """

        # Variable to keep track of row number so we don't have to!
        row = 0

        # Add a message box
        self.message_box = tk.Text(self.master, height = 7, width=50)
        self.message_box.grid(row=0, column=4, rowspan=5, padx=padx, sticky='NS', )
        self.message_box.insert(tk.END, "Click 'Browse...' to load path to netCDF file.\n")

        #Add scrollbar to message box
        self.scroll = tk.Scrollbar(self.master)
        self.scroll.grid(row=0, column=5, pady=pady, sticky='NS')
        self.scroll.config(command=self.message_box.yview)
        self.message_box.config(yscrollcommand=self.scroll.set)

        #Make new frame for file stuff
        self.netcdf_frame = tk.Frame(self.master, bg=cname)
        self.netcdf_frame.grid(row=row, column=0, columnspan=4, padx=padx, pady=pady, sticky='NSWE')
        tk.Label(self.netcdf_frame, text='Input Files:', bg=cname, width=col1_width).grid(row=row, column=1, sticky='E')

        #Browse Button to browse for netcdf files
        self.netcdf_button = tk.Button(self.netcdf_frame, text='Browse...',width=col3_width)
        self.netcdf_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.netcdf_button['command'] = self.browse_netcdf

        # Checkbox to allow combining repeat runs' NERODS3 together.
        self.repeat = tk.BooleanVar()
        self.repeat.set(False)
        self.repeat_selector = tk.Checkbutton(self.netcdf_frame, variable=self.repeat, bg=cname2, text='Repeat Runs?', onvalue=True, offvalue=False)
        self.repeat_selector.grid(row=row, column=3, padx=padx, pady=pady)
        self.repeat_selector['command'] = self.repeat_command

        row += 1

        #Place for netcdf file
        self.netcdf_entry = tk.Entry(self.netcdf_frame, width=col2_width)
        self.netcdf_entry.grid(row=row, column=0, columnspan=4, padx=padx, pady=pady, sticky='WE')

        tk.Label(self.netcdf_frame, text='.nc', bg=cname, width=col3_width, anchor='w').grid(row=row, column=4)

        row += 1

        #Place for dat file
        self.dat_entry = tk.Entry(self.netcdf_frame, width=col2_width)
        self.dat_entry.grid(row=row, column=0, columnspan=4, padx=padx, pady=pady, sticky='WE')

        tk.Label(self.netcdf_frame, text='.dat', bg=cname, width=col3_width, anchor='w').grid(row=row, column=4)

        row += 1

        #Place for lim file
        self.lim_entry = tk.Entry(self.netcdf_frame, width=col2_width)
        self.lim_entry.grid(row=row, column=0, columnspan=4, padx=padx, pady=pady, sticky='WE')

        tk.Label(self.netcdf_frame, text='.lim', bg=cname, width=col3_width, anchor='w').grid(row=row, column=4)

        row += 1

        #Plot overview button to plot simple graphs
        self.overview_button = tk.Button(self.netcdf_frame, text='Plot Overview')
        self.overview_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.overview_button['command'] = self.overview

        row += 1

        #Make a new frame for plot selection
        self.sel_frame = tk.Frame(self.master, bg=cname2)
        self.sel_frame.grid(row=row, column=0, padx=padx, pady=pady, sticky='WE')
        tk.Label(self.sel_frame, text='Plot Selection:', bg=cname2).grid(row=row, column=1)

        #Make a new frame for plot options
        self.opt_frame = tk.Frame(self.master, bg=cname)
        self.opt_frame.grid(row=row, column=4, padx=padx, pady=pady, sticky='NSWE')
        tk.Label(self.opt_frame, text='Additional Plot Options', bg=cname).grid(row=1, column=1)
        self.opt_frame.grid_rowconfigure(0, weight=1)
        self.opt_frame.grid_rowconfigure(2, weight=1)
        self.opt_frame.grid_columnconfigure(0, weight=1)
        self.opt_frame.grid_columnconfigure(2, weight=1)

        #row += 1

        #Option Menu for different plots
        self.current_option = tk.StringVar(self.sel_frame)
        self.current_option.set(plot_op[0])
        self.plot_options = tk.OptionMenu(self.sel_frame, self.current_option, *plot_op)
        self.plot_options.grid(row=row, column=2, padx=padx, pady=pady)
        self.current_option.trace('w', self.option_selection)

        #Power T image
        #self.image = Image.open("PowerTimage.png")
        #self.image = self.image.resize((100,100))
        #self.render = ImageTk.PhotoImage(image=self.image)
        #self.img = tk.Label(self.sel_frame, image=self.render)
        #self.img.image = self.render
        #self.img.grid(row=row-1, column=4, padx=padx, pady=pady, sticky='e', rowspan=3)

        row += 1

        #multiplot selector to plot multiple plots
        self.multiplot = tk.IntVar()
        self.multiplot.set(0)
        self.multiplot_selector = tk.Checkbutton(self.sel_frame, variable=self.multiplot, bg=cname2, text='Multiplot: ')
        self.multiplot_selector.grid(row=row, column=1, padx=padx, pady=pady)
        self.multiplot_selector['command'] = self.multiplot_tab

        row += 1

        #Plot button to plot the selected plot
        self.plot_button = tk.Button(self.sel_frame, text='Plot', width=col3_width)
        self.plot_button.grid(row=row, column=1, padx=padx, pady=pady)
        self.plot_button['command'] = self.plot_action

        row += 1

        self.quit_button = tk.Button(self.master, text='Quit', width=col1_width)
        self.quit_button.grid(row=row, column=0, padx=padx, pady=pady, rowspan=2)
        self.quit_button['command'] = self.quit_command

    def overview(self):
        """
        Function to plot the overview plot
        """

        # checks to see if there is a file selected
        if self.netcdf_entry.get() == '':
            self.message_box.insert(tk.END, 'No netCDF File Selected\n')

        else:
            self.dl.overviewplot()

    def plot_arguments(self):
        """
        Function to get the plot arguments for plotting the graphs with
        """

        if self.current_option.get() == 'Center Line':
            try:
                if self.centline_log == 1:
                    plot_args = {'log': True}
                else:
                    plot_args = {'log': False}

                if self.centline_exp ==1:
                    plot_args['fit_exp'] = True
                else:
                    plot_args['fit_exp'] = False
            except:
                plot_args = {}

        elif self.current_option.get() == 'Contour':

            if self.side_option.get() == 2:
                plot_args = {'side': 'ITF'}
            else:
                plot_args = {'side': 'OTF'}

            if self.width_option.get() == 2:
                plot_args['probe_width'] = 0.005
            elif self.width_option.get() == 3:
                plot_args['probe_width'] = 0.0025
            elif self.width_option.get() == 4:
                try:
                    plot_args['probe_width'] = float(self.cont_widEnt.get())
                except:
                    self.message_box.insert(tk.END, 'Custom width failed.\n')
                    plot_args['probe_width'] = 0.015
            else:
                plot_args['probe_width'] = 0.015

        elif self.current_option.get() == 'Poloidal Profiles':

            if self.width_option.get() == 2:
                plot_args = {'probe_width': 0.005}
            elif self.width_option.get() == 3:
                plot_args = {'probe_width': 0.0025}
            elif self.width_option.get() == 4:
                try:
                    plot_args = {'probe_width': float(self.cont_widEnt.get())}
                except:
                    self.message_box.insert(tk.END, 'Custom width failed.\n')
                    plot_args = {'probe_width': 0.015}
            else:
                plot_args = {'probe_width': 0.015}

        elif self.current_option.get() == 'Temperature Contour':

            plot_args = {}

        elif self.current_option.get() == 'Plasma Density':

            plot_args = {}

        elif self.current_option.get() == 'Impurity Velocity':

            plot_args = {}

        elif self.current_option.get() == 'Impurity Contour':

            plot_args = {}

        else:
            plot_args = {}

        return plot_args

    def plot_action(self):
        """
        Function to take current options and plot the selected graph.
        """

        #checks to see if there is a file selected
        if self.netcdf_entry.get() == '':
            self.message_box.insert(tk.END, 'No netCDF File Selected\n')
            return

        if self.multiplot.get() == 0:
            plot_args = self.plot_arguments()

            if self.current_option.get() == 'Center Line':

                self.dl.centerline(**plot_args)

            elif self.current_option.get() == 'Contour':

                self.dl.deposition_contour(**plot_args)

            elif self.current_option.get() == 'Poloidal Profiles':

                self.dl.avg_pol_profiles(**plot_args)

            elif self.current_option.get() == 'Temperature Contour':

                self.dl.te_contour(**plot_args)

            elif self.current_option.get() == 'Plasma Density':

                self.dl.ne_contour(**plot_args)

            elif self.current_option.get() == 'Impurity Velocity':

                self.dl.avg_imp_vely(**plot_args)

            elif self.current_option.get() == 'Impurity Contour':

                self.dl.imp_contour_plot(**plot_args)

            elif self.current_option.get() == 'Please Select Option':

                self.message_box.insert(tk.END, 'Please Select a Plot Option \n')

            else:
                self.message_box.insert(tk.END, 'Plotting Failed \n')

        elif self.multiplot == 1:

            self.dl.muliplot_end()

    def option_selection(self, var, ind, mode):
        """
        Function to make selection of plot and show additional plot options.
        """

        #The Center Line is selected this is what happens
        if self.current_option.get() == 'Center Line':

            #creates new frame for the center line selection
            self.delete_frames()
            self.opt_centline = tk.Frame(self.master, bg=cname)
            self.opt_centline.grid(row=5, column=4, columnspan=4, padx=padx, pady=pady, sticky='NSWE')

            #Adds log options checkbox
            tk.Label(self.opt_centline, text='Log plot', bg=cname, width=10, anchor='w').grid(row=1, column=2)
            self.log_option = tk.IntVar()
            self.centline_log = tk.Checkbutton(self.opt_centline, variable=self.log_option, onvalue=1, offvalue=0, bg=cname)
            self.centline_log.grid(row=1, column=1)

            #Adds exp options checkbox
            tk.Label(self.opt_centline, text='Exponential fit', bg=cname, anchor='w').grid(row=2, column=2)
            self.exp_option = tk.IntVar()
            self.centline_exp = tk.Checkbutton(self.opt_centline, variable=self.exp_option, onvalue=1, offvalue=0, bg=cname)
            self.centline_exp.grid(row=2, column=1)

            #Centers checkboxes in frame
            self.opt_centline.grid_rowconfigure(0, weight=1)
            self.opt_centline.grid_rowconfigure(2, weight=1)

        #The contour is selected this is what happens
        elif self.current_option.get() == 'Contour':

            #creates new frame for the contour selection
            self.delete_frames()
            self.opt_cont = tk.Frame(self.master, bg=cname)
            self.opt_cont.grid(row=5, column=4, columnspan=4, padx=padx, pady=pady, sticky='NSWE')

            tk.Label(self.opt_cont, text='Side:', bg=cname).grid(row=5, column=4, columnspan=2)

            tk.Label(self.opt_cont, text='Probe Width:', bg=cname).grid(row=5, column=6, columnspan=2)

            #next row

            #Adds side selection OTF
            tk.Label(self.opt_cont, text='OTF: ', width=10, bg=cname, anchor='e').grid(row=6, column=4)
            self.side_option = tk.IntVar()
            self.side_option.set(1)
            self.cont_OTF = tk.Radiobutton(self.opt_cont, variable=self.side_option, value=1, bg=cname)
            self.cont_OTF.grid(row=6, column=5)

            tk.Label(self.opt_cont, text='0.015: ', width=10, bg=cname, anchor='e').grid(row=6, column=6)
            self.width_option = tk.IntVar()
            self.width_option.set(1)
            self.cont_widA = tk.Radiobutton(self.opt_cont, variable=self.width_option, value=1, bg=cname)
            self.cont_widA.grid(row=6, column=7)

            #next row

            #Adds side selection ITF
            tk.Label(self.opt_cont, text='ITF: ', width=10, bg=cname, anchor='e').grid(row=7, column=4)
            self.cont_ITF = tk.Radiobutton(self.opt_cont, variable=self.side_option, value=2, bg=cname)
            self.cont_ITF.grid(row=7, column=5)

            tk.Label(self.opt_cont, text='0.005: ', width=10, bg=cname, anchor='e').grid(row=7, column=6)
            self.cont_widB = tk.Radiobutton(self.opt_cont, variable=self.width_option, value=2, bg=cname)
            self.cont_widB.grid(row=7, column=7)

            #next row

            tk.Label(self.opt_cont, text='0.0025: ', width=10, bg=cname, anchor='e').grid(row=8, column=6)
            self.cont_widC = tk.Radiobutton(self.opt_cont, variable=self.width_option, value=3, bg=cname)
            self.cont_widC.grid(row=8, column=7)

            #next row

            tk.Label(self.opt_cont, text='Other: ', width=10, bg=cname, anchor='e').grid(row=9, column=6)
            self.cont_widO = tk.Radiobutton(self.opt_cont, variable=self.width_option, value=4, bg=cname)
            self.cont_widO.grid(row=9, column=7)

            self.cont_widEnt = tk.Entry(self.opt_cont)
            self.cont_widEnt.grid(row=9, column=8)

        elif self.current_option.get() == 'Poloidal Profiles':

            self.delete_frames()
            self.opt_polo = tk.Frame(self.master, bg=cname)
            self.opt_polo.grid(row=5, column=4, columnspan=4, padx=padx, pady=pady, sticky='NSWE')

            tk.Label(self.opt_polo, text='',bg=cname, width=10).grid(row=5, column=4, columnspan=2)

            tk.Label(self.opt_polo, text='Probe Width:', bg=cname).grid(row=5, column=6, columnspan=2)

            tk.Label(self.opt_polo, text='0.015: ', width=10, bg=cname, anchor='e').grid(row=6, column=6)
            self.width_option = tk.IntVar()
            self.width_option.set(1)
            self.cont_widA = tk.Radiobutton(self.opt_polo, variable=self.width_option, value=1, bg=cname)
            self.cont_widA.grid(row=6, column=7)

            tk.Label(self.opt_polo, text='0.005: ', width=10, bg=cname, anchor='e').grid(row=7, column=6)
            self.cont_widB = tk.Radiobutton(self.opt_polo, variable=self.width_option, value=2, bg=cname)
            self.cont_widB.grid(row=7, column=7)

            tk.Label(self.opt_polo, text='0.0025: ', width=10, bg=cname, anchor='e').grid(row=8, column=6)
            self.cont_widC = tk.Radiobutton(self.opt_polo, variable=self.width_option, value=3, bg=cname)
            self.cont_widC.grid(row=8, column=7)

            tk.Label(self.opt_polo, text='Other: ', width=10, bg=cname, anchor='e').grid(row=9, column=6)
            self.cont_widO = tk.Radiobutton(self.opt_polo, variable=self.width_option, value=4, bg=cname)
            self.cont_widO.grid(row=9, column=7)

            self.cont_widEnt = tk.Entry(self.opt_polo)
            self.cont_widEnt.grid(row=9, column=8)

        #The temperature is selected this is what happens
        elif self.current_option.get() == 'Temperature Contour':

            #creates new frame for the temperature selection
            self.delete_frames()
            self.opt_Temp = tk.Frame(self.master, bg=cname)
            self.opt_Temp.grid(row=5, column=4, padx=padx, pady=pady, sticky='WENS')

            #Adds centered text in frame
            tk.Label(self.opt_Temp, text='No Extra Options for Temperature Contour', bg=cname).grid(row=1, column=1)
            self.opt_Temp.grid_rowconfigure(0, weight=1)
            self.opt_Temp.grid_rowconfigure(2, weight=1)
            self.opt_Temp.grid_columnconfigure(0, weight=1)
            self.opt_Temp.grid_columnconfigure(2, weight=1)

        elif self.current_option.get() == 'Plasma Density':

            #creates new frame for the Plasma Density selection
            self.delete_frames()
            self.opt_Plasma = tk.Frame(self.master, bg=cname)
            self.opt_Plasma.grid(row=5, column=4, padx=padx, pady=pady, sticky='NSEW')

            #Adds centered text in frame
            tk.Label(self.opt_Plasma, text='No Extra Options for Plasma Density', bg=cname).grid(row=1, column=1)
            self.opt_Plasma.grid_rowconfigure(0, weight=1)
            self.opt_Plasma.grid_rowconfigure(2, weight=1)
            self.opt_Plasma.grid_columnconfigure(0, weight=1)
            self.opt_Plasma.grid_columnconfigure(2, weight=1)

        elif self.current_option.get() == 'Impurity Contour':

            #creates new frame for the Impurtiy Contour selection
            self.delete_frames()
            self.opt_ImpCon = tk.Frame(self.master, bg=cname)
            self.opt_ImpCon.grid(row=5, column=4, columnspan=4, padx=padx, pady=pady, sticky='NSEW')

            #Adds rmin entry box
            tk.Label(self.opt_ImpCon, text='rmin: ', bg=cname).grid(row=1, column=1, padx=padx, pady=pady)
            tk.Entry(self.opt_ImpCon).grid(row=1, column=2, padx=padx, pady=pady)

    def delete_frames(self):
        try:
            self.opt_frame.grid_forget()
        except:
            pass

        try:
            self.opt_Temp.grid_forget()
        except:
            pass

        try:
            self.opt_polo.grid_forget()
        except:
            pass

        try:
            self.opt_cont.grid_forget()
        except:
            pass

        try:
            self.opt_centline.grid_forget()
        except:
            pass

        try:
            self.opt_Plasma.grid_forget()
        except:
            pass

        try:
            self.opt_ImpCon.grid_forget()
        except:
            pass

    def multiplot_tab(self):
        """
        Function to either open or close the multiplot tab
        """

        if self.multiplot.get() == 0:
            self.mul_frame.grid_forget()

        elif self.multiplot.get() == 1:

            #self.dl.multiplot_start()

            row = 0

            self.mul_frame = tk.Frame(self.master, bg=cname)
            self.mul_frame.grid(row=0, rowspan=8, column=6, sticky='NS')

            tk.Label(self.mul_frame, text='Multiplot Selector:', bg=cname).grid(row=row, column=1)

            row += 1

            self.plot1_button = tk.Button(self.mul_frame, text='Plot 1')
            self.plot1_button.grid(row=row, column=0, padx=padx, pady=pady)
            self.plot1_button['command'] = self.plot1_command

            self.plot2_button = tk.Button(self.mul_frame, text='Plot 2')
            self.plot2_button.grid(row=row, column=1, padx=padx, pady=pady)

            self.plot3_button = tk.Button(self.mul_frame, text='Plot 3')
            self.plot3_button.grid(row=row, column=2, padx=padx, pady=pady)

            row +=1

            self.plot4_button = tk.Button(self.mul_frame, text='Plot 4')
            self.plot4_button.grid(row=row, column=0, padx=padx, pady=pady)

            self.plot5_button = tk.Button(self.mul_frame, text='Plot 5')
            self.plot5_button.grid(row=row, column=1, padx=padx, pady=pady)

            self.plot6_button = tk.Button(self.mul_frame, text='Plot 6')
            self.plot6_button.grid(row=row, column=2, padx=padx, pady=pady)

            row += 1

            self.plot7_button = tk.Button(self.mul_frame, text='Plot 7')
            self.plot7_button.grid(row=row, column=0, padx=padx, pady=pady)

            self.plot8_button = tk.Button(self.mul_frame, text='Plot 8')
            self.plot8_button.grid(row=row, column=1, padx=padx, pady=pady)

            self.plot9_button = tk.Button(self.mul_frame, text='Plot 9')
            self.plot9_button.grid(row=row, column=2, padx=padx, pady=pady)

            row += 1

            tk.Label(self.mul_frame, text='', bg=cname).grid(row=row, column=1)

            row += 1

            tk.Label(self.mul_frame, text='Selected Plots:', bg=cname).grid(row=row, column=1, padx=padx, pady=pady)

            row += 1

            self.plot1_label = tk.Label(self.mul_frame, text='N/A').grid(row=row, column=0, padx=padx, pady=pady)
            self.plot2_label = tk.Label(self.mul_frame, text='N/A').grid(row=row, column=1, padx=padx, pady=pady)
            self.plot3_label = tk.Label(self.mul_frame, text='N/A').grid(row=row, column=2, padx=padx, pady=pady)

            row += 1

            self.plot4_label = tk.Label(self.mul_frame, text='N/A').grid(row=row, column=0, padx=padx, pady=pady)
            self.plot5_label = tk.Label(self.mul_frame, text='N/A').grid(row=row, column=1, padx=padx, pady=pady)
            self.plot6_label = tk.Label(self.mul_frame, text='N/A').grid(row=row, column=2, padx=padx, pady=pady)

            row += 1

            self.plot7_label = tk.Label(self.mul_frame, text='N/A').grid(row=row, column=0, padx=padx, pady=pady)
            self.plot8_label = tk.Label(self.mul_frame, text='N/A').grid(row=row, column=1, padx=padx, pady=pady)
            self.plot9_label = tk.Label(self.mul_frame, text='N/A').grid(row=row, column=2, padx=padx, pady=pady)

    def repeat_command(self):

        if self.repeat.get() == 0:
            combine_repeat_runs = False
        elif self.repeat.get() == 1:
            combine_repeat_funs = True

    def plot1_command(self):

        plot_args = self.plot_arguments()
        plot_args = {'plotnum': 1}

        if self.current_option.get() == 'Center Line':

            self.dl.centerline(**plot_args)

        elif self.current_option.get() == 'Contour':

            self.dl.deposition_contour(**plot_args)

        elif self.current_option.get() == 'Poloidal Profiles':

            self.dl.avg_pol_profiles(**plot_args)

    def browse_netcdf(self):
        """
        Function for Browse button to get netcdf file.
        """

        #adds the netcdf file
        self.message_box.insert(tk.END, 'Loading...\n')
        root = tk.Tk(); root.withdraw()
        netcdf_path = tk.filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),))
        self.netcdf_entry.delete(0, tk.END)
        self.netcdf_entry.insert(0, netcdf_path)
        self.dl = limpt.LimPlots(self.netcdf_entry.get(), combine_repeat_runs=self.repeat.get())
        self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(netcdf_path.split('/')[-1]))

        cp_path = netcdf_path.split('.nc')[0] + '.collector_probe'
        try:
            f = open(cp_path, 'r')
            self.cp_entry.delete(0, tk.END)
            self.cp_entry.insert(0, cp_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(cp_path.split('/')[-1]))
            self.cp_path = cp_path
        except:
            pass

        #adds the dat file
        dat_path = netcdf_path.split('.nc')[0] + '.dat'
        try:
            f = open(dat_path, 'r')
            self.dat_entry.delete(0, tk.END)
            self.dat_entry.insert(0, dat_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(dat_path.split('/')[-1]))
            self.dat_path = dat_path
        except:
            pass

        #adds the lim file
        lim_path = netcdf_path.split('.nc')[0] + '.lim'
        try:
            f = open(lim_path, 'r')
            self.lim_entry.delete(0,tk.END)
            self.lim_entry.insert(0, lim_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(lim_path.split('/')[-1]))
            self.lim_path = lim_path
        except:
            pass

    def quit_command(self):

        self.master.quit()
        self.master.destroy()


def main():

    plt.ion()
    root = tk.Tk()
    window = Window(root)
    window.mainloop()

if __name__ == '__main__':
    main()
