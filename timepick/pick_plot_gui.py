
def manual_pick_grid(self,root):
    r''' 
    Insert a label frame into the master frame to get manual
    picking and plot parameters
    We put the label frame into row 0 of mframe
    we put plot_frame into  row 1 of mframe
    '''    

    import tkinter as tk
    from tkinter import ttk
    import procvsp.utils_gui as utilg
    from timepick import pick_plot_gui

    # remove the other parameter grids or they will get overlayed by this new one
    # need to find a more comprehensive way
    utilg.grid_clear(self)
  
    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.mpick_paramframe = ttk.LabelFrame(self.mframe, text = 'Manual Picker Options')#, height = 5, width = 800)
    self.mpick_paramframe.grid(row = 0, column = 0, padx = 2, pady=2, sticky='nw')

    # data plotting frame - insert text for user input boxes
    labelscal = tk.Label(master = self.mpick_paramframe,text="Amp. Multiplier" )
    labelfsplt = tk.Label(master = self.mpick_paramframe,text="Samp. Rate (Hz)" )
    labelsplt= tk.Label(master = self.mpick_paramframe,text="Start Time (ms)" )
    labeleplt = tk.Label(master = self.mpick_paramframe,text="End Time (ms)" )
    labelstr = tk.Label(master = self.mpick_paramframe,text="Start Trace" )
    labeletr = tk.Label(master = self.mpick_paramframe,text="End Trace" )
    labeltpi = tk.Label(master = self.mpick_paramframe,text="Trace Scale" )
    labelskip = tk.Label(master = self.mpick_paramframe,text="Header Deci." )
    labelsubsmp = tk.Label(master = self.mpick_paramframe,text="Subsample." )

    labelscal.grid(row=0, column=0,padx=5, pady=5, sticky='w')
    labelfsplt.grid(row=1, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labelsplt.grid(row=0, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labeleplt.grid(row=1, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelstr.grid(row=0, column=4,ipadx=0,padx=5, pady=5, sticky='w')
    labeletr.grid(row=1, column=4, ipadx=0,padx=5, pady=5, sticky='w')
    labeltpi.grid(row=1, column=6, ipadx=0,padx=5, pady=5, sticky='w')
    labelskip.grid(row=0, column=6, ipadx=0,padx=5, pady=5, sticky='w')
    labelsubsmp.grid(row=0, column=8, ipadx=0,padx=5, pady=5, sticky='w')

   # data plotting frame- get user inputs, insert defaults
    self.entryscal = tk.Entry(master=self.mpick_paramframe, width=8)
    self.entryfsplt = tk.Entry(master=self.mpick_paramframe, width=8)
    self.entrysplt = tk.Entry(master=self.mpick_paramframe, width=8)
    self.entryeplt = tk.Entry(master=self.mpick_paramframe, width=8)
    self.entrystr = tk.Entry(master=self.mpick_paramframe, width=8)
    self.entryetr = tk.Entry(master=self.mpick_paramframe, width=8)
    self.entrytpi = tk.Entry(master=self.mpick_paramframe, width=8)
    self.entryskip = tk.Entry(master=self.mpick_paramframe, width=8)
    self.entrysubsmp = tk.Entry(master=self.mpick_paramframe, width=8)
    
    self.entryscal.grid(row=0, column=1, padx=5, pady=5, sticky='w')
    self.entryfsplt.grid(row=1, column=1, padx=5, pady=5, sticky='w')
    self.entrysplt.grid(row=0, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryeplt.grid(row=1, column=3, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrystr.grid(row=0, column=5, ipadx=0,padx=5, pady=5, sticky='w')
    self.entryetr.grid(row=1, column=5, ipadx=0, padx=5, pady=5, sticky='w')
    self.entrytpi.grid(row=1, column=7, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryskip.grid(row=0, column=7, ipadx=0, padx=5, pady=5, sticky='w')
    self.entrysubsmp.grid(row=0, column=9, ipadx=0, padx=5, pady=5, sticky='w')    

    #buttontrplt = ttk.Button(self.mpick_paramframe, text = "Picking Window", 
    #                         command = lambda: pick_plot_gui.mpick_plot(self,root))
    #
    buttontrplt = ttk.Button(self.mpick_paramframe, text = "Picking Window", 
                             command = lambda: pick_plot_gui.mpick_plot(self,root))         
    buttontrplt.grid(row=0, column=13, columnspan=1, sticky='w',padx=10, pady=5 )

    # data plotting frame - create checkbuttons for plotting controls
    self.normchk = tk.StringVar()
    self.polchk = tk.StringVar()
    self.vachk = tk.StringVar()
    self.decichk = tk.StringVar()
    self.grdchk = tk.StringVar()
    self.update = tk.StringVar()
    
    self.normchk.set('y')
    self.polchk.set('n')
    self.vachk.set('n')
    self.decichk.set('n')
    self.grdchk.set('n')
    self.update.set('n') 

    cknorm = ttk.Checkbutton(self.mpick_paramframe, text='Normalize', variable=self.normchk, onvalue='y', offvalue = 'n')
    ckpol = ttk.Checkbutton(self.mpick_paramframe, text='Flip Pol.', variable=self.polchk, onvalue='r', offvalue = 'n')
    ckva = ttk.Checkbutton(self.mpick_paramframe, text='VA plot', variable=self.vachk, onvalue='y', offvalue = 'n')
    ckdeci = ttk.Checkbutton(self.mpick_paramframe, text='Decimate', variable=self.decichk, onvalue='y', offvalue = 'n')
    ckgrd = ttk.Checkbutton(self.mpick_paramframe, text='Grid Lines', variable=self.grdchk, onvalue='y', offvalue = 'n')
    ckupdate = ttk.Checkbutton(self.mpick_paramframe, text='Update Field Headers', variable=self.update, onvalue='y', offvalue = 'n')

    cknorm.grid(row=0, column=10, ipadx=0, padx=5, pady=5, sticky='w')
    ckpol.grid(row=1, column=10, ipadx=0, padx=5, pady=5, sticky='w')
    ckva.grid(row=0, column=11, ipadx=0, padx=5, pady=5, sticky='w')
    ckdeci.grid(row=1, column=11, ipadx=0, padx=5, pady=5, sticky='w')
    ckgrd.grid(row=0, column=12, ipadx=0, padx=5, pady=5, sticky='w')
    ckupdate.grid(row=1, column=12, ipadx=0, padx=5, pady=5, sticky='w')

   # Set some default values which are updated by SEG-D header values
    self.entryfsplt.insert(0, '1000')
    self.entrysplt.insert(0, '0')
    self.entryeplt.insert(0, '2000')
    self.entrystr.insert(0, '1')
    self.entryetr.insert(0, '10')
    self.entryscal.insert(0, '1')
    self.entrytpi.insert(0,'10')
    self.entryskip.insert(0,'1')
    self.entrysubsmp.insert(0,'1')    
    # Insert default values that are read from SEG-Y data file
    self.entryfsplt.delete(0, tk.END)
    self.entryfsplt.insert(0, '%s'%(self.fs))
    self.entrysplt.delete(0, tk.END)
    self.entrysplt.insert(0, '0')
    self.entrystr.delete(0, tk.END)
    self.entrystr.insert(0, '1')
    self.entryetr.delete(0, tk.END)
    self.entryetr.insert(0, '%s'%(self.plotdata.shape[0]))

    self.mpick_menu_exists = self.mpick_paramframe.winfo_exists()

def mpick_get_params(self):
    ''' create a dictionary of parameters read from the manual trace 
    picking option.

    THese are passed to picker_multitrace in pick_plot.py
    '''

    mpik_params=dict()

    mpik_params['Scalar'] = float(self.entryscal.get())
    trange = [float(self.entrysplt.get()),float(self.entryeplt.get())]
    mpik_params['time_range'] = trange
    tracerange = [int(self.entrystr.get())-1,int(self.entryetr.get())]
    mpik_params['trace_range'] = tracerange
    mpik_params['polarity'] = self.polchk.get()
    mpik_params['header_spacing'] = int(self.entryskip.get())
    mpik_params['plot_norm'] = self.normchk.get()
    mpik_params['var_area'] = self.vachk.get()
    mpik_params['tr_per_in'] = int(self.entrytpi.get())
    mpik_params['subsamp_ratio'] = int(self.entrysubsmp.get())
    mpik_params['update_h'] = self.update.get()
    mpik_params['tkinter_use'] = 'y'

    return mpik_params

def man_pick(self,root):
    ''' Get manual picking paramaters
    1. set up tpick_param menu
    '''
    from timepick import pick_plot_gui

    pick_plot_gui.manual_pick_grid(self,root)

def header_update(self,mpikframe):
    ''' Copy the manually picked headers in column 16 to
    the input trace header file column 8
    Overwrites original times with STA-LTA with manual picks
    '''
    import numpy as np
    if self.headers.shape[1]!=15:
        print ('\n Header Save button results ')       
        print (' self.headers.shape :',self.headers.shape)
        diff =np.max( np.abs(self.headers[:,8]-self.picks[:,16]))
        print ('\n header update worked STA LTA, diff : ',diff)    
        self.headers[:,8] = self.picks[:,16]
    else:
        self.headers[:,8] = self.picks[:,8]
        diff =np.max( np.abs(self.headers[:,8]-self.picks[:,8]))
        print ('\n header update worked , diff : ',diff)

    mpikframe.destroy()

def test_plot_canvas(self,root,pickparams,*qcdata,**kwargs2):
    ''' STA LTA test plots on one VSP trace
    
    Separate canvas from figure for re-usability of plotting code.

    1. Create a tkinter top level window 
    2. Crete a tkinter canvas with scrollbars place on toplevel 
    3. Create a matplotlib figure
    4. Place figure on canvas
    '''
    import procvsp.utils_gui as utilg
    import timepick.pick_plot as pickplot
    import tkinter as tk    

    # create a toplevel window (pop-out) to plot response
    fname='STA LTA QC'
    staqcframe=utilg.create_toplevel(root,fname)
    # create a tkinter canvas and scroll bars
    canvas_qcsta=utilg.add_scrollbars(staqcframe)
    # create the matplotlib figure
    figt, ax= pickplot.test_plots_guts(pickparams,*qcdata, **kwargs2) 
    figt.canvas.draw_idle()  # should be faster then plt.show but in reality is not noticeable
    # copy the figure to the canvas
    canvas_qcsta = utilg.fig_canvas(figt,canvas_qcsta)
    # Bring window with figure to front
    tk.Misc.lift(staqcframe)

def mpick_plot(self,root):
    ''' Manual time picking window
    
    Separate canvas from figure for re-usability of plotting code.

    1. Create a tkinter top level window 
    2. Crete a tkinter canvas with scrollbars place on toplevel 
    3. Create a matplotlib figure
    4. Place figure on canvas
    '''
    import procvsp.utils_gui as utilg
    import timepick.pick_plot as pickplot
    import tkinter as tk
    from timepick import pick_plot_gui
    
    # get the plotting parameters from gui grid
    manual_params=pick_plot_gui.mpick_get_params(self)
    # create a toplevel window (pop-out) to plot response
    fname='Manual Time Picking'
    mpikframe=utilg.create_toplevel(root,fname)
    # Bring window with figure to front
    tk.Misc.lift(mpikframe)
    # create a tkinter canvas and scroll bars
    canvas_mpik=utilg.add_scrollbars(mpikframe)
    # add a save button to plot frame    
    save_button = tk.Button(mpikframe, text="Save Headers", command=lambda: header_update(self,mpikframe))
    save_button.grid(row = 2, column = 0, )     
    # create the matplotlib figure
    figm, self.picks= pickplot.picker_multitrace(self.data,self.headers,self.fs,**manual_params)

    figm.canvas.draw_idle()  # should be faster then plt.show but in reality is not noticeable
    # copy the figure to the canvas
    canvas_mpik = utilg.fig_canvas(figm,canvas_mpik)

    #mpikframe.destroy()


