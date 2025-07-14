r'''1d and 2d Spectral plots.
Tkinter widgets to get parameters
'''
def nextpow2(x):
    import numpy as np
 
    ''' find a number that is a power of 2
    from https://stackoverflow.com/questions/14267555/find-the-
    smallest-power-of-2-greater-than-n-in-python
    ''' 
    if x == 0:            
        y = 1
       
    else: y = np.ceil(np.log2(x))
        
    print("\u0332".join('\nnextpow2 info and parameters') )
    print(' x :', x)
    print(' nextpow2 :', y)    

    return int(y)
    
def spectra(idata, timerange, frange, fs,twin):

    import numpy as np
    import scipy.signal as sig
    import scipy.fft    

    ####### window the trace to prevent edge effect in fft  #############
    
    dt =1/fs *1000             # sample rate in ms
    samprate = 1/fs            #sample rate in seconds

    ####### optimize number of samples for transform  #############

    # for a segment of trace
    if (twin =='y')or(twin=='Y'):
        start = int(timerange[0]*(fs/1000))
        stop = int(timerange[1]*(fs/1000))
        idata_trimd = idata[:,start:stop]
        
        # Apply window to segment
        w=sig.tukey(idata_trimd.shape[1], alpha=.1)    
        idata_win = idata_trimd[:,:]*w.reshape(-1, w.shape[0])        # Multiply trace by window
        N = scipy.fft.next_fast_len(int((timerange[1]-timerange[0])*(fs/1000)))
    
    # for a full trace
    else:
        # Apply window to whole trace
        w=sig.tukey(idata.shape[1], alpha=.1)    
        idata_win = idata[:,:]*w.reshape(-1, w.shape[0])        # Multiply trace by window    
        N = scipy.fft.next_fast_len(idata_win.shape[1]) # in samples, best for scipy fft
    
    # pad with zeros if optimal N greater than trace or segment
    #if (N > idata_win.shape[1]):                    
    #    pad = N - idata_win.shape[1]
    #    idata_win = np.pad(idata_win, ((0,0),(0,int(pad))), 'constant')        

    ####################### do the fft  #####################################
    X = np.zeros(shape = (idata_win.shape[0], N),dtype=np.complex_)
    X_db = np.zeros(shape = (idata_win.shape[0],N))
    X_pow = np.zeros(shape = (idata_win.shape[0],N))

    for k in range(0,(idata_win.shape[0])):
        
        X[k,:] = scipy.fft.fft(idata_win[k,:], n = N)      # from 0, to TT plus window
        X_db[k,:] = 20*np.log10(np.abs(X[k,:])/(np.abs(np.max(X[0,:])))) # db=20*np.log10(S/np.max(S)
        #X_pow[k,:] = np.abs(X[k,:])**2 # power spectrum
    # Only keep positive frequencies #########    
    freq = scipy.fft.fftfreq(N, d=samprate)    # Generate plot frequency axis

    keep = freq>=0
    freq = freq[keep]
    
    X_posfreq = np.zeros(shape = (idata_win.shape[0], freq.shape[0]))
    X_db_posfreq = np.zeros(shape = (idata_win.shape[0], freq.shape[0]))
    #X_pow_posfreq = np.zeros(shape = (idata_win.shape[0], freq.shape[0]))
    
    X_posfreq = X[:, keep]
    X_db_posfreq = X_db[:, keep]
    #X_pow_posfreq = X_pow[:, keep]
    
    return X_posfreq, X_db_posfreq,freq

def spec_param_grid(self,root):
    r''' 
    Insert a label frame into the master frame to get spectral display
    parameters. 
    
    We put the label frame into row 0 of mframe
    we put plot_frame into  row 1 of mframe
    '''    
    
    import tkinter as tk
    from tkinter import ttk
    import procvsp.utils_gui as utilg
    from procvsp import spec_gui


    # clear any existing parameter grids
    utilg.grid_clear(self)
    
    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.spec_paramframe = ttk.LabelFrame(self.mframe, text = 'Spectral Display Options')#, height = 5, width = 800)
    self.spec_paramframe.grid(row = 0, column = 0, padx = 2, pady=2, sticky='nw')

    # Bandpass filter frame - insert text for user input boxes
    labelsfreq = tk.Label(master = self.spec_paramframe,text="Start Freq. (Hz)" )
    labelefreq = tk.Label(master = self.spec_paramframe,text="End Freq. (Hz)" )
    labelstim = tk.Label(master = self.spec_paramframe,text="Start Time (ms)" )
    labeletim = tk.Label(master = self.spec_paramframe,text="End Time (ms)" )
    labelspectrc = tk.Label(master = self.spec_paramframe,text="Analysis Trace" )
    labelsfreq.grid(row=0, column=0, ipadx=0,padx=5, pady=5, sticky='w')
    labelefreq.grid(row=1, column=0, ipadx=0,padx=5, pady=5, sticky='w')
    labelstim.grid(row=0, column=2, ipadx=0,padx=5, pady=5, sticky='w')
    labeletim.grid(row=1, column=2, ipadx=0,padx=5, pady=5, sticky='w')
    labelspectrc.grid(row=1, column=4, ipadx=0,padx=5, pady=5, sticky='w')

    # Bandpass filter frame- get user inputs, insert defaults
    self.entrysfreq = tk.Entry(master=self.spec_paramframe, width=8)
    self.entryefreq = tk.Entry(master=self.spec_paramframe, width=8)
    self.entrystim = tk.Entry(master=self.spec_paramframe, width=8)
    self.entryetim = tk.Entry(master=self.spec_paramframe, width=8)
    self.entryspectrc = tk.Entry(master=self.spec_paramframe, width=8)

    self.entrysfreq.insert(0, '5')
    self.entryefreq.insert(0, '15')
    self.entrystim.insert(0, '0')
    self.entryetim.insert(0, '2000')
    self.entryspectrc.insert(0, '1')

    self.entrysfreq.grid(row=0, column=1, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryefreq.grid(row=1, column=1, ipadx=0, padx=5,pady=5, sticky='w')
    self.entrystim.grid(row=0, column=3, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryetim.grid(row=1, column=3, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryspectrc.grid(row=1, column=5, ipadx=0, padx=5,pady=5, sticky='w')
    # Insert default values that are read from SEG-Y data file
    self.entryetim.delete(0, tk.END)
    self.entryetim.insert(0, '%s'%(self.plotdata.shape[1]))
    self.entryefreq.delete(0, tk.END)
    self.entryefreq.insert(0, '%s'%(int(self.fs//2))) # default to nyquist frequency
    
    # data plotting frame - create checkbuttons for plotting controls
    self.window = tk.StringVar()       
    self.window.set('y')    
    cktwin = ttk.Checkbutton(self.spec_paramframe, text='Apply Window?', variable=self.window, onvalue='y', offvalue = 'n')    
    cktwin.grid(row=0, column=7, ipadx=0, padx=5, pady=5, sticky='w')
    
    #1d spectrum frame - create action  buttons
    buttonwindow = ttk.Button(self.spec_paramframe, text = "Generate spectra ", command = lambda: spec_gui.spectra1d(self,root))
    buttonwindow.grid(row=1, column=7, columnspan=1,sticky='w',padx=10, pady=5 )
    
    self.spec_menu_exists = self.spec_paramframe.winfo_exists()

def spec1d_get_params(self):
    ''' create a dictionary of parameters read from the spectrum 
    options

    THese are passed to spectra1d 
    '''
    spc1_params=dict()

    spc1_params['time_range'] = [float(self.entrystim.get()),float(self.entryetim.get())]
    spc1_params['freq_range'] = [float(self.entrysfreq.get()),float(self.entryefreq.get())]
    spc1_params['trace'] = int(self.entryspectrc.get())
    spc1_params['time_win'] = self.window.get()
    spc1_params['title_fran'] = ' VSP trace'
    spc1_params['tkinter_use'] = 'y'

    return spc1_params

def spectrum(self,root):
        ''' display a 1D T-F spectrum
        1. set up spec_param menu
        The menu allows starting the analysis
        '''
        from procvsp import spec_gui

        spec_gui.spec_param_grid(self,root)
    
def spectra1d(self,root):
    ''' Spectra from one VSP trace
    
    Separate canvas from figure for re-usability of plotting code.

    1. Create a tkinter top level window 
    2. Crete a tkinter canvas with scrollbars place on toplevel 
    3. Create a matplotlib figure
    4. Place figure on canvas
    '''
    import procvsp.utils_gui as utilg
    import procvsp.spec as spec
    import tkinter as tk 
    # get parameters from gui
    spec_params = spec1d_get_params(self)   
    # create a toplevel window (pop-out) to plot filter QC
    fname='Trace Spectra'
    spec1dframe=utilg.create_toplevel(root,fname)
    # create a tkinter canvas and scroll bars
    canvas_spec1d=utilg.add_scrollbars(spec1dframe)
    # create the matplotlib figure
    figt = spec.spec_1d(self.data,self.headers,self.fs, **spec_params) 
    figt.canvas.draw_idle()  # should be faster then plt.show but in reality is not noticeable
    # copy the figure to the canvas
    canvas_spec1d = utilg.fig_canvas(figt,canvas_spec1d)
    # Bring window with figure to front
    tk.Misc.lift(spec1dframe)    

def fzspec_param_grid(self,root):
    r''' 
    Insert a label frame into the master frame to get 2d (FZ) spectral 
    display parameters. 
    
    We put the label frame into row 0 of mframe
    we put plot_frame into  row 1 of mframe
    '''    
    
    import tkinter as tk
    from tkinter import ttk
    import procvsp.utils_gui as utilg
    from procvsp import spec_gui


    # clear any existing parameter grids
    utilg.grid_clear(self)
    
    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.fzspec_paramframe = ttk.LabelFrame(self.mframe, text = ' VSP Spectral Display Options')#, height = 5, width = 800)
    self.fzspec_paramframe.grid(row = 0, column = 0, padx = 2, pady=2, sticky='nw')

    # Bandpass filter frame - insert text for user input boxes
    labelsfreq = tk.Label(master = self.fzspec_paramframe,text="Start Freq. (Hz)" )
    labelefreq = tk.Label(master = self.fzspec_paramframe,text="End Freq. (Hz)" )
    labelstim = tk.Label(master = self.fzspec_paramframe,text="Start Time (ms)" )
    labeletim = tk.Label(master = self.fzspec_paramframe,text="End Time (ms)" )
    labelmindb = tk.Label(master = self.fzspec_paramframe,text="Min. db)" )
    labelmaxdb = tk.Label(master = self.fzspec_paramframe,text="Max. db" )
    labelftrc = tk.Label(master = self.fzspec_paramframe,text="First trc.)" )
    labeletrc = tk.Label(master = self.fzspec_paramframe,text="Last trc." )
    labelscal = tk.Label(master = self.fzspec_paramframe,text="Multiplier" )
    labelskip = tk.Label(master = self.fzspec_paramframe,text="Skip nth." )
    
    labelsfreq.grid(row=0, column=0, ipadx=0,padx=5, pady=5, sticky='w')
    labelefreq.grid(row=1, column=0, ipadx=0,padx=5, pady=5, sticky='w')
    labelstim.grid(row=0, column=2, ipadx=0,padx=5, pady=5, sticky='w')
    labeletim.grid(row=1, column=2, ipadx=0,padx=5, pady=5, sticky='w')
    labelmindb.grid(row=0, column=4, ipadx=0,padx=5, pady=5, sticky='w')
    labelmaxdb.grid(row=1, column=4, ipadx=0,padx=5, pady=5, sticky='w')
    labelftrc.grid(row=0, column=6, ipadx=0,padx=5, pady=5, sticky='w')
    labeletrc.grid(row=1, column=6, ipadx=0,padx=5, pady=5, sticky='w')
    labelscal.grid(row=0, column=8, ipadx=0,padx=5, pady=5, sticky='w')
    labelskip.grid(row=1, column=8, ipadx=0,padx=5, pady=5, sticky='w')

    # Bandpass filter frame- get user inputs, insert defaults
    self.entrysfreq = tk.Entry(master=self.fzspec_paramframe, width=8)
    self.entryefreq = tk.Entry(master=self.fzspec_paramframe, width=8)
    self.entrystim = tk.Entry(master=self.fzspec_paramframe, width=8)
    self.entryetim = tk.Entry(master=self.fzspec_paramframe, width=8)
    self.entrymindb = tk.Entry(master=self.fzspec_paramframe, width=8)
    self.entrymaxdb = tk.Entry(master=self.fzspec_paramframe, width=8)
    self.entryftrc = tk.Entry(master=self.fzspec_paramframe, width=8)
    self.entryetrc = tk.Entry(master=self.fzspec_paramframe, width=8)
    self.entryscal = tk.Entry(master=self.fzspec_paramframe, width=8)
    self.entryskip = tk.Entry(master=self.fzspec_paramframe, width=8)

    self.entrysfreq.insert(0, '5')
    self.entryefreq.insert(0, '15')
    self.entrystim.insert(0, '0')
    self.entryetim.insert(0, '2000')
    self.entrymindb.insert(0, '-90')
    self.entrymaxdb.insert(0, '0')
    self.entryftrc.insert(0, '1')
    self.entryetrc.insert(0, '1')    
    self.entryscal.insert(0, '1')    
    self.entryskip.insert(0, '1')    

    self.entrysfreq.grid(row=0, column=1, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryefreq.grid(row=1, column=1, ipadx=0, padx=5,pady=5, sticky='w')
    self.entrystim.grid(row=0, column=3, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryetim.grid(row=1, column=3, ipadx=0, padx=5,pady=5, sticky='w')
    self.entrymindb.grid(row=0, column=5, ipadx=0, padx=5,pady=5, sticky='w')
    self.entrymaxdb.grid(row=1, column=5, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryftrc.grid(row=0, column=7, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryetrc.grid(row=1, column=7, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryscal.grid(row=0, column=9, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryskip.grid(row=1, column=9, ipadx=0, padx=5,pady=5, sticky='w')

    # Insert default values that are read from SEG-Y data file
    self.entryetim.delete(0, tk.END)
    self.entryetim.insert(0, '%s'%(self.plotdata.shape[1]))
    self.entryefreq.delete(0, tk.END)
    self.entryefreq.insert(0, '%s'%(int(self.fs//2))) # default to nyquist frequency
    self.entryetrc.delete(0, tk.END)
    self.entryetrc.insert(0, '%s'%(self.plotdata.shape[0]))
    
    # data plotting frame - create checkbuttons for plotting controls
    self.interpchk = tk.StringVar()       
    self.interpchk.set('n')    
    ckinterp = ttk.Checkbutton(self.fzspec_paramframe, text='Gaussian Smoothing?', variable=self.interpchk, onvalue='y', offvalue = 'n')    
    ckinterp.grid(row=0, column=10, ipadx=0, padx=5, pady=5, sticky='w')

    self.dbchk = tk.StringVar()       
    self.dbchk.set('y')    
    ckdb = ttk.Checkbutton(self.fzspec_paramframe, text='db plot?', variable=self.dbchk, onvalue='y', offvalue = 'n')    
    ckdb.grid(row=1, column=10, ipadx=0, padx=5, pady=5, sticky='w')
    
    #1d spectrum frame - create action  buttons
    buttonwindow = ttk.Button(self.fzspec_paramframe, text = "Generate spectra ", command = lambda: spec_gui.spectra2d(self,root))
    buttonwindow.grid(row=1, column=11, columnspan=1,sticky='w',padx=10, pady=5 )
    
    self.fzspec_menu_exists = self.fzspec_paramframe.winfo_exists()

def specfz_get_params(self):

    specfz_params = dict()
    specfz_params['spec_type'] = self.dbchk.get() # 'y' for dB, 'n' for amplitude
    
    specfz_params['scale'] = float(self.entryscal.get()) # scale up image apmlitude plot
    specfz_params['title_fran'] = ' VSP Frequency-Depth %sms to %sms'%(float(self.entrystim.get()),
                                                                       float(self.entryetim.get()))
    specfz_params['time_win'] = 'y' # window trace prior to fft, default to yes
    specfz_params['time_range'] = [float(self.entrystim.get()),float(self.entryetim.get())]
    specfz_params['freq_range'] = [float(self.entrysfreq.get()),float(self.entryefreq.get())]    
    specfz_params['trace_range'] = [int(self.entryftrc.get()),int(self.entryetrc.get())]
    specfz_params['db_range'] = [float(self.entrymindb.get()),float(self.entrymaxdb.get())] 
    specfz_params['header_skip'] = int(float(self.entryskip.get())) # plot every nth header
    specfz_params['plt_rotate'] = 'n' # rotate the plot 90 degrees if rotate =='y'
    specfz_params['savepng'] = 'n'
    specfz_params['tkinter_use'] = 'n'
    specfz_params['interp'] = self.interpchk.get()

    return specfz_params

def spectrumFZ(self,root):
    ''' display a 2D T-F spectrum one
    for each trace, color coded.
    1. set up fzspec_param menu
    The menu allows starting the analysis
    '''
    from procvsp import spec_gui

    spec_gui.fzspec_param_grid(self,root)

def spectra2d_dummy(self,root):
    ''' do the forward transforms for each trace
    and display the results
    '''
    from procvsp import spec_gui
    spec_gui.spec_FZ_SLB(self,root)

def spectra2d(self,root):
    ''' Spectra from one VSP trace
    
    Separate canvas from figure for re-usability of plotting code.

    1. Create a tkinter top level window 
    2. Crete a tkinter canvas with scrollbars place on toplevel 
    3. Create a matplotlib figure
    4. Place figure on canvas
    '''
    import procvsp.utils_gui as utilg
    import procvsp.spec as spec
    import tkinter as tk 
    # get parameters from gui
    spec2d_params = specfz_get_params(self)   
    # create a toplevel window (pop-out) to plot filter QC
    fname='Trace Spectra'
    spec2dframe=utilg.create_toplevel(root,fname)
    # create a tkinter canvas and scroll bars
    canvas_spec2d=utilg.add_scrollbars(spec2dframe)
    # create the matplotlib figure
    fig2d = spec.spec_FZ(self.data, self.headers, self.fs,**spec2d_params) 
    fig2d.canvas.draw_idle()  # should be faster then plt.show but in reality is not noticeable
    # copy the figure to the canvas
    canvas_spec2d = utilg.fig_canvas(fig2d,canvas_spec2d)
    # Bring window with figure to front
    tk.Misc.lift(spec2dframe)    

def bpf_param_grid(self,root):
    r''' 
    Insert a label frame into the master frame to get BPF parameters and 
    plot parameters

    There is an option to re-plot the filtered seismic.

    We put the label frame into row 0 of mframe
    we put plot_frame into  row 1 of mframe

    QC plots are on a toplevel window (pop-out)
    '''    
    import tkinter as tk
    from tkinter import ttk
    import procvsp.utils_gui as utilg
    from procvsp import spec_gui
    from plotvsp import seisplots_gui
    
    # clear any existing parameter grids
    utilg.grid_clear(self)

    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.bpf_paramframe = ttk.LabelFrame(self.mframe, text = 'Butterworth Filter Options')#, height = 5, width = 800)
    self.bpf_paramframe.grid(row = 0, column = 0, padx = 2, pady=2, sticky='nw')

    # Bandpass filter frame - insert text for user input boxes
    labelsbpf = tk.Label(master = self.bpf_paramframe,text="Low Cut Freq. (Hz)" )
    labelebpf = tk.Label(master = self.bpf_paramframe,text="High Cut Freq. (Hz)" )
    labelordbpf = tk.Label(master = self.bpf_paramframe,text="Butterworth Order" )
    labelsbpf.grid(row=0, column=0, ipadx=0,padx=5, pady=5, sticky='w')
    labelebpf.grid(row=0, column=2, ipadx=0,padx=5, pady=5, sticky='w')
    labelordbpf.grid(row=1, column=0, ipadx=0,padx=5, pady=5, sticky='w')
    # Bandpass filter frame- get user inputs, insert defaults
    self.entrysbpf = tk.Entry(master=self.bpf_paramframe, width=8)
    self.entryebpf = tk.Entry(master=self.bpf_paramframe, width=8)
    self.entryordbpf = tk.Entry(master=self.bpf_paramframe, width=8)
    self.entrysbpf.insert(0, '5')
    self.entryebpf.insert(0, '15')
    self.entryordbpf.insert(0, '2')
    self.entrysbpf.grid(row=0, column=1, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryebpf.grid(row=0, column=3, ipadx=0, padx=5,pady=5, sticky='w')
    self.entryordbpf.grid(row=1, column=1, ipadx=0, padx=5,pady=5, sticky='w')
    # set default end frequency to nyquist
    self.entryebpf.delete(0, tk.END)
    self.entryebpf.insert(0, '%s'%(int(self.fs//2)-1)) # default to nyquist frequency

    #Bandpass filter frame - create action  buttons
    buttonbpfapp = ttk.Button(self.bpf_paramframe, text = "Run Filter", 
                              command = lambda: spec_gui.filter(self,root))
    buttonbpfapp.grid(row=0, column=4, columnspan=1,sticky='w',padx=10, pady=5 )
    buttonbpfplt = ttk.Button(self.bpf_paramframe, text = "Filtered Plot Parameters", 
                              command = lambda: seisplots_gui.plot_stack(self))
    buttonbpfplt.grid(row=1, column=4, columnspan=1,sticky='w',padx=10, pady=5 )
    
    self.bpf_menu_exists = self.bpf_paramframe.winfo_exists()

def bpf_get_params(self):
    ''' create a dictionary of parameters read from the manual trace 
    picking option.

    THese are passed to picker_multitrace in pick_plot.py
    '''

    bpf_params=dict()

    bpf_params['lowcut'] = int(self.entrysbpf.get())
    bpf_params['highcut'] = int(self.entryebpf.get())
    bpf_params['order'] = int(self.entryordbpf.get())
    bpf_params['polarity'] = self.polchk.get()
    bpf_params['numfsamp'] = 1024
    bpf_params['qcplot'] = 'y'
    bpf_params['tkinter_use'] = 'y'

    return bpf_params

def bpf_filt(self,root):
    ''' apply a Butterworth bandpass filter
    1. set up BPF_param menu
    The menu allows applying the filter, then plotting the data
    '''
    from procvsp import spec_gui
    # do not try to delete as plot has deleted it
    #self.spec_menu_exists = 0

    spec_gui.bpf_param_grid(self,root)
    
def filter(self,root):
    ''' BPF with test plots on one VSP trace
    
    Separate canvas from figure for re-usability of plotting code.

    1. Create a tkinter top level window 
    2. Crete a tkinter canvas with scrollbars place on toplevel 
    3. Create a matplotlib figure
    4. Place figure on canvas
    '''
    from plotvsp import seisplots_gui
    import procvsp.utils_gui as utilg
    import procvsp.spec as spec
    import tkinter as tk 
    # get parameters from gui
    bpf_params = bpf_get_params(self)  
    # apply the filter
    self.bpf, sos_out=spec.bandpass_filter(self.data, self.fs, **bpf_params)
    # plot everything with defaults from data shape and SEG-Y headers
    # no plotting parameter grid is generated
    self.plotdata = self.bpf
    self.data = self.bpf
    self.plot_all='y'
    self.sta_pick = 'n'
    # make a plot with all data - no options
    seisplots_gui.stack_plot(self) 
    # create a toplevel window (pop-out) to plot filter QC
    fname='Band Pass QC'
    specqcframe=utilg.create_toplevel(root,fname)
    # create a tkinter canvas and scroll bars
    canvas_qcspec=utilg.add_scrollbars(specqcframe)
    # create the matplotlib figure
    figt = spec.BPF_QCplot(sos_out,self.fs, **bpf_params) 
    figt.canvas.draw_idle()  # should be faster then plt.show but in reality is not noticeable
    # copy the figure to the canvas
    canvas_qcspec = utilg.fig_canvas(figt,canvas_qcspec)
    # Bring window with figure to front
    tk.Misc.lift(specqcframe)
    
      