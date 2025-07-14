r'''plot routines for gui presentations
'''
def plot_param_top_grid(self):
    r''' 
    Insert a label frame into the master frame to get plot parameters
    We put the label frame into row 0 of mframe
    we put plot_frame into  row 1 of mframe
    '''    
    # create a toplevel window (pop-out) to get plot parameters
    # self.master is the root window in a form which can be passed to functions
    
    import tkinter as tk
    from tkinter import ttk
    import procvsp.utils_gui as utilg
    from plotvsp import seisplots_gui

    # remove the other parameter grids or they will get overlayed by this new one
    # need to find a more comprehensive way
    utilg.grid_clear(self)
  
    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.plot_paramframe = ttk.LabelFrame(self.mframe, text = 'Trace Display Options')#, height = 5, width = 800)
    self.plot_paramframe.grid(row = 0, column = 0, padx = 2, pady=2, sticky='nw')

    # data plotting frame - insert text for user input boxes
    labelscal = tk.Label(master = self.plot_paramframe,text="Amp. Multiplier" )
    labelfsplt = tk.Label(master = self.plot_paramframe,text="Samp. Rate (Hz)" )
    labelsplt= tk.Label(master = self.plot_paramframe,text="Start Time (ms)" )
    labeleplt = tk.Label(master = self.plot_paramframe,text="End Time (ms)" )
    labelstr = tk.Label(master = self.plot_paramframe,text="Start Trace" )
    labeletr = tk.Label(master = self.plot_paramframe,text="End Trace" )
    labeltpi = tk.Label(master = self.plot_paramframe,text="Trace Scale" )
    labelskip = tk.Label(master = self.plot_paramframe,text="Header Deci." ) 

    labelscal.grid(row=0, column=0,padx=5, pady=5, sticky='w')
    labelfsplt.grid(row=1, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labelsplt.grid(row=0, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labeleplt.grid(row=1, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelstr.grid(row=0, column=4,ipadx=0,padx=5, pady=5, sticky='w')
    labeletr.grid(row=1, column=4, ipadx=0,padx=5, pady=5, sticky='w')
    labeltpi.grid(row=1, column=6, ipadx=0,padx=5, pady=5, sticky='w')
    labelskip.grid(row=0, column=6, ipadx=0,padx=5, pady=5, sticky='w')

   # data plotting frame- get user inputs, insert defaults
    self.entryscal = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryfsplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrysplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryeplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrystr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryetr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrytpi = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryskip = tk.Entry(master=self.plot_paramframe, width=8)
    
    self.entryscal.grid(row=0, column=1, padx=5, pady=5, sticky='w')
    self.entryfsplt.grid(row=1, column=1, padx=5, pady=5, sticky='w')
    self.entrysplt.grid(row=0, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryeplt.grid(row=1, column=3, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrystr.grid(row=0, column=5, ipadx=0,padx=5, pady=5, sticky='w')
    self.entryetr.grid(row=1, column=5, ipadx=0, padx=5, pady=5, sticky='w')
    self.entrytpi.grid(row=1, column=7, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryskip.grid(row=0, column=7, ipadx=0, padx=5, pady=5, sticky='w')

    # data plotting frame - create action  buttons
    #if self.load_type == 1:
    #    buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.raw_plot)
    #elif self.load_type == 2:
    #    buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.xcorr_plot)
    #elif self.load_type == 3:
    buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = lambda: seisplots_gui.stack_plot(self))        
    buttontrplt.grid(row=0, column=12, columnspan=1, sticky='w',padx=10, pady=5 )
    # data plotting frame - create checkbuttons for plotting controls
    self.normchk = tk.StringVar()
    self.polchk = tk.StringVar()
    self.vachk = tk.StringVar()
    self.decichk = tk.StringVar()
    self.grdchk = tk.StringVar()
    
    self.normchk.set('y')
    self.polchk.set('n')
    self.vachk.set('n')
    self.decichk.set('n')
    self.grdchk.set('n')

    cknorm = ttk.Checkbutton(self.plot_paramframe, text='Normalize', variable=self.normchk, onvalue='y', offvalue = 'n')
    ckpol = ttk.Checkbutton(self.plot_paramframe, text='Flip Pol.', variable=self.polchk, onvalue='r', offvalue = 'n')
    ckva = ttk.Checkbutton(self.plot_paramframe, text='VA plot', variable=self.vachk, onvalue='y', offvalue = 'n')
    ckdeci = ttk.Checkbutton(self.plot_paramframe, text='Decimate', variable=self.decichk, onvalue='y', offvalue = 'n')
    ckgrd = ttk.Checkbutton(self.plot_paramframe, text='Grid Lines', variable=self.grdchk, onvalue='y', offvalue = 'n')

    cknorm.grid(row=0, column=9, ipadx=0, padx=5, pady=5, sticky='w')
    ckpol.grid(row=1, column=9, ipadx=0, padx=5, pady=5, sticky='w')
    ckva.grid(row=0, column=10, ipadx=0, padx=5, pady=5, sticky='w')
    ckdeci.grid(row=1, column=10, ipadx=0, padx=5, pady=5, sticky='w')
    ckgrd.grid(row=0, column=11, ipadx=0, padx=5, pady=5, sticky='w')

   # Set some default values which are updated by SEG-D header values
    self.entryfsplt.insert(0, '1000')
    self.entrysplt.insert(0, '0')
    self.entryeplt.insert(0, '2000')
    self.entrystr.insert(0, '1')
    self.entryetr.insert(0, '10')
    # Set some default values which can be manually updated
    self.entryscal.insert(0, '1')
    #self.entryips.insert(0,'5')
    self.entrytpi.insert(0,'10')
    self.entryskip.insert(0,'1')    
    # Insert default values that are read from SEG-Y data file
    self.entryfsplt.delete(0, tk.END)
    self.entryfsplt.insert(0, '%s'%(self.fs))
    self.entrysplt.delete(0, tk.END)
    self.entrysplt.insert(0, '0')
    self.entryeplt.delete(0, tk.END)
    self.entryeplt.insert(0, '%s'%(self.plottindex.max()))
    self.entrystr.delete(0, tk.END)
    self.entrystr.insert(0, '1')
    self.entryetr.delete(0, tk.END)
    self.entryetr.insert(0, '%s'%(self.plotdata.shape[0]))

    self.plot_menu_exists = self.plot_paramframe.winfo_exists()

def scrlplot_get_params(self):
    ''' Get prameters from gui for scrolling plot of traces
    '''
    scrplt=dict()

    if self.plot_all=='y':
        scrplt['first'] = 0
        scrplt['last'] = int(self.plotdata.shape[0])
        # get user input plot parameters, .get returns strings!
        scrplt['scal'] = 1
        scrplt['stime'] = 0
        scrplt['etime'] = int(self.plotdata.shape[1]*1000/self.fs)
        scrplt['sta_pick'] = 'n'
        # use start time and end time from STA-LTA as plot limit
        if self.sta_pick == 'y':
            scrplt['sta_pick'] = 'y'
            scrplt['stime'] = float(self.entrytims.get())
            scrplt['etime'] = float(self.entrytime.get())
        scrplt['pol'] = 'n'
        scrplt['norm'] = 'y'
        scrplt['va'] = 'n'
        scrplt['skiplabel'] = 1
        scrplt['tpi'] = 10
        scrplt['dec'] = 'y'
        scrplt['ygrids'] = 'n' 
        scrplt['plot_all'] = 'y'    
    else:
        scrplt['first'] = int(self.entrystr.get())-1
        scrplt['last'] = int(self.entryetr.get())
        scrplt['scal'] = float(self.entryscal.get())
        scrplt['stime'] = int(self.entrysplt.get())
        scrplt['etime'] = int(float(self.entryeplt.get()))
        scrplt['pol'] = self.polchk.get()
        scrplt['norm'] = self.normchk.get()
        scrplt['va'] = self.vachk.get()
        scrplt['dec'] = self.decichk.get()
        scrplt['skiplabel'] = int(self.entryskip.get())
        scrplt['tpi'] = int(self.entrytpi.get())
        scrplt['ygrids'] = self.grdchk.get()
        scrplt['plot_all'] = self.plot_all
        scrplt['sta_pick'] = self.sta_pick

    return scrplt

def plot_stack(self):
    ''' set up the grid for plotting parameter inputs
    '''
    from plotvsp import seisplots_gui

    # to use widgets to get plot parameters, set to 'n'
    self.plot_all='n'
    # use start time and end time from STA-LTA as plot limit
    self.sta_pick = 'n'
    # parameter grid setup
    seisplots_gui.plot_param_top_grid(self)

def stack_plot(self):
    '''Plot data traces
    
    Separate canvas from figure for re-usability of plotting code.

    1. Create a tkinter top level window 
    2. Crete a tkinter canvas with scrollbars place on toplevel 
    3. Create a matplotlib figure
    4. Place figure on canvas
    '''
    import procvsp.utils_gui as utilg
    import tkinter as tk
    from plotvsp import seisplots_gui
    from plotvsp import seisplots
    
    # get the plotting parameters from gui grid or
    # from data extents if plot_all='y'
    scroll_params=seisplots_gui.scrlplot_get_params(self)
    # create a tkinter canvas and scroll bars
    canvas_scroll=utilg.add_scrollbars(self.plot_frame)
    # create the matplotlib figure
    figs = seisplots.scroll_plot(self.headers,self.data,self.fs,**scroll_params)
    figs.canvas.draw_idle()  # should be faster then plt.show but in reality is not noticeable
    # copy the figure to the canvas
    canvas_scroll = utilg.fig_canvas(figs,canvas_scroll)
    # Bring window plotting parameters to front
    if self.plot_all=='n':
        tk.Misc.lift(self.plot_paramframe)