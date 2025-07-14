
def segy_top_grid(self):
    ''' create a toplevel window (pop-out) to get plot parameters
    # self.master is the root window in a form which can be passed to functions
    may need to get tt byte, ttscalar, DF_ASL, SrcElev, SRD_ASL
    '''
    import tkinter as tk
    from tkinter import ttk
    import procvsp.utils_gui as utilg
    #import timepick.gui_funcs as guifuncs
    from iovsp import segyin_gui


    # if segy in frame (or any other)already exists, delete it and 
    # create a new one in it's place
    
    utilg.grid_clear(self)

    # create a fancy frame for segy load parameter entry in column 0, row 1 of master window
    self.segyin_paramframe = ttk.LabelFrame(self.mframe, text = 'Header Adjustment Options')#, height = 5, width = 800)
    self.segyin_paramframe.grid(row = 0, column = 0, padx = 2, pady=2, sticky='nw')

    # geometry updating variables frame - insert text for user input boxes
    labelttst = tk.Label(master = self.segyin_paramframe,text="TT start byte" )
    labelttsc = tk.Label(master = self.segyin_paramframe,text="Scalar to adjust TT" )
    labeldfel= tk.Label(master = self.segyin_paramframe,text="Drill Floor Elevation" )
    labelsrcel = tk.Label(master = self.segyin_paramframe,text="Source Elevation" )
    labelsrdel = tk.Label(master = self.segyin_paramframe,text="SRD Elevation" )

    labelttst.grid(row=0, column=0,padx=5, pady=5, sticky='w')
    labelttsc.grid(row=1, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labeldfel.grid(row=0, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelsrcel.grid(row=1, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelsrdel.grid(row=0, column=4,ipadx=0,padx=5, pady=5, sticky='w')
    
   # geometry updating frame- get user inputs, insert defaults
    self.entryttst = tk.Entry(master=self.segyin_paramframe, width=8)
    self.entryttsc = tk.Entry(master=self.segyin_paramframe, width=8)
    self.entrydfel = tk.Entry(master=self.segyin_paramframe, width=8)
    self.entrysrcel = tk.Entry(master=self.segyin_paramframe, width=8)
    self.entrysrdel = tk.Entry(master=self.segyin_paramframe, width=8)
   
    self.entryttst.grid(row=0, column=1, padx=5, pady=5, sticky='w')
    self.entryttsc.grid(row=1, column=1, padx=5, pady=5, sticky='w')
    self.entrydfel.grid(row=0, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entrysrcel.grid(row=1, column=3, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrysrdel.grid(row=0, column=5, ipadx=0, padx=5, pady=5, sticky='w')

   # Set some default values which are updated by SEG-Y header values
    self.entryttst.insert(0, '197')
    self.entryttsc.insert(0, '100')
    self.entrydfel.insert(0, '0')
    self.entrysrcel.insert(0, '0')
    self.entrysrdel.insert(0, '0')

    # create checkbuttons for writing tables
    self.hdrschk = tk.StringVar()
    self.tblchk = tk.StringVar()
    
    self.hdrschk.set('y')
    self.tblchk.set('n')

    buttontrplt = ttk.Button(self.segyin_paramframe, text = "Read SEG-Y", command = segyin_gui.pick_segy(self))#self.pick_segy)        
    buttontrplt.grid(row=0, column=6, columnspan=1, sticky='w',padx=10, pady=5 )   
 
    self.segyin_menu_exists = self.segyin_paramframe.winfo_exists()

def get_segyin_params(self):
    ''' Get SEG-Y input variables from gui and write to dictionary which
    can be passed around.
    '''    
    segyin_params=dict()
    segyin_params['SrcElev'] = int(self.entrysrcel.get())
    segyin_params['DF_ASL'] = float(self.entrydfel.get())
    segyin_params['SRD_ASL'] = float(self.entrysrdel.get())
    segyin_params['ttbyte'] = int(self.entryttst.get())
    segyin_params['ttscalar'] = float(self.entryttsc.get())
    segyin_params['file_headers'] = 'n'
    segyin_params['PTS'] = 'n'

    # variable to indicate tkinter useage
    segyin_params['tkinter_use']='y'
    
    return segyin_params

def segyread(self):
    ''' set up and plot parameter entry grid.
    Action button invokes pick_segy'''

    from iovsp import segyin_gui

    segyin_gui.segy_top_grid(self)

def pick_segy(self):
    ''' started from segyin action button. Read the segy file
    and plot all the traces and samples'''    

    import numpy as np
    from iovsp import segyin, segyin_gui
    from tkinter import filedialog as fd
    from plotvsp import seisplots_gui
    
    # place parameter grid on main window
    segyin_params=get_segyin_params(self)

    # use a windows dialog box to select a seg-y file
    segy_path = fd.askopenfilename(initialdir = "c:\\Users\\acampbell45\\Documents",
                title = "Select file",filetypes = (("seg-y files","*.sgy"),("all files","*.*")))
    # self.segy is referenced in some older functions
    self.segy = segy_path

    #self.data, self.numsamp, self.samprate, self.fs, self.headers = segyin_gui.readsegy_gui(self)
    self.data, self.numsamp, self.samprate, self.fs, self.headers = segyin.readsegyio3(segy_path,**segyin_params)

    self.tindex = np.arange(0, (self.numsamp+1)*(1000/self.fs),(1000/self.fs) )  # convert fs to msec.
    
    self.plotdata = np.copy(self.data)
    self.plottindex = np.copy(self.tindex)
    self.plotthead = np.copy(self.headers)
    
    # plot everything with defaults from data shape and SEG-Y headers
    # no plotting parameter grid is generated
    self.plot_all='y'
    self.sta_pick = 'n'
    # make a plot with all data - no options
    seisplots_gui.stack_plot(self)

    print ('\n Number of traces, Number of samples :',self.data.shape[0],self.data.shape[1])
    print (' Number of traces, Number of headers :',self.headers.shape[0],self.headers.shape[1])
    
    # call the plot options with parameter grid for updating plot
    self.plot_all='n'
    seisplots_gui.plot_param_top_grid(self) 
        