''' Routines to integrate pickers into gui.

I avoided using objects to make the code more easily useable in notebooks etc.   
ie. did not want multiple copies of functions, too hard to update, find errors etc.

'''

def sta_params_grid(self,root):
    r''' 
    Insert a label frame into the master frame to get STA LTA parameters and 
    plot parameters
    We put the label frame into row 0 of mframe
    we put plot_frame into  row 1 of mframe
    '''    
    import tkinter as tk
    from tkinter import ttk
    import procvsp.utils_gui as utilg
    import timepick.pickers_gui as tpick
    from plotvsp import seisplots_gui

    # remove the other parameter grids or they will get overlayed by this new one
    # need to find a more comprehensive way
    utilg.grid_clear(self)
  
    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.tpick_paramframe = ttk.LabelFrame(self.mframe, text = 'Time Picking Options')#, height = 5, width = 800)
    self.tpick_paramframe.grid(row = 0, column = 0, padx = 2, pady=2, sticky='nw')

    # data plotting frame - insert text for user input boxes
    labeltims = tk.Label(master = self.tpick_paramframe,text="Start Time (ms)" )
    labeltime = tk.Label(master = self.tpick_paramframe,text="End Time (ms)" )
    labelsta= tk.Label(master = self.tpick_paramframe,text="STA Window (ms)" )
    labellta = tk.Label(master = self.tpick_paramframe,text="LTA Window (ms)" )
    labelttrc = tk.Label(master = self.tpick_paramframe,text="Test Trace" )
    labelsubrat = tk.Label(master = self.tpick_paramframe,text="Sub Sample Ratio" )
    labeltune = tk.Label(master = self.tpick_paramframe,text="Tune Window" )
    labelstab = tk.Label(master = self.tpick_paramframe,text="Stab Factor." )
    labelfs =  tk.Label(master = self.tpick_paramframe,text="Sample Rate %s hertz - info only"%(self.fs) )
    
    labeltims.grid(row=0, column=0,padx=5, pady=5, sticky='w')
    labeltime.grid(row=1, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labelsta.grid(row=0, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labellta.grid(row=1, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelstab.grid(row=0, column=4, ipadx=0,padx=5, pady=5, sticky='w')
    labelsubrat.grid(row=1, column=4, ipadx=0,padx=5, pady=5, sticky='w')
    labeltune.grid(row=0, column=6, ipadx=0,padx=5, pady=5, sticky='w')
    labelttrc.grid(row=1, column=6,ipadx=0,padx=5, pady=5, sticky='w')
    labelfs.grid(row=0, column=16, ipadx=0,padx=5, pady=5, sticky='w')

   # data plotting frame- get user inputs, insert defaults
    self.entrytims = tk.Entry(master=self.tpick_paramframe, width=8)
    self.entrytime = tk.Entry(master=self.tpick_paramframe, width=8)
    self.entrysta = tk.Entry(master=self.tpick_paramframe, width=8)
    self.entrylta = tk.Entry(master=self.tpick_paramframe, width=8)
    self.entryttrc = tk.Entry(master=self.tpick_paramframe, width=8)
    self.entrysubrat = tk.Entry(master=self.tpick_paramframe, width=8)
    self.entrytune = tk.Entry(master=self.tpick_paramframe, width=8)
    self.entrystab = tk.Entry(master=self.tpick_paramframe, width=8)
    
    self.entrytims.grid(row=0, column=1, padx=5, pady=5, sticky='w')
    self.entrytime.grid(row=1, column=1, padx=5, pady=5, sticky='w')
    self.entrysta.grid(row=0, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entrylta.grid(row=1, column=3, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrystab.grid(row=0, column=5, ipadx=0, padx=5, pady=5, sticky='w')
    self.entrysubrat.grid(row=1, column=5, ipadx=0, padx=5, pady=5, sticky='w')
    self.entrytune.grid(row=0, column=7, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryttrc.grid(row=1, column=7, ipadx=0,padx=5, pady=5, sticky='w')

    buttonttpick = ttk.Button(self.tpick_paramframe, text = "Run STA LTA", 
                              command = lambda: tpick.tpicks_sta(self,root))        
    buttonttpick.grid(row=0, column=15, columnspan=1, sticky='w',padx=10, pady=5 )

    buttontrplt = ttk.Button(self.tpick_paramframe, text = "Make Plot", 
                             command = lambda: seisplots_gui.plot_stack(self))        
    buttontrplt.grid(row=1, column=15, columnspan=1, sticky='w',padx=10, pady=5 )
    
    # create checkbuttons for plotting controls 
    
    self.hedup = tk.StringVar()
    self.qcplt = tk.StringVar()
    self.hilenv = tk.StringVar()
    self.peaks = tk.StringVar()
    self.verbose = tk.StringVar()
       
    self.hedup.set('n')
    self.qcplt.set('n')
    self.hilenv.set('n')
    self.peaks.set('y')
    self.verbose.set('n')

    ckhedupdate = ttk.Checkbutton(self.tpick_paramframe, text='Update Headers ?', variable=self.hedup, onvalue='y', offvalue = 'n')
    ckhedqcplot = ttk.Checkbutton(self.tpick_paramframe, text='QC plots ?', variable=self.qcplt, onvalue='y', offvalue = 'n')
    ckhilbert = ttk.Checkbutton(self.tpick_paramframe, text='Use Hilbert Envelope ?', variable=self.hilenv, onvalue='y', offvalue = 'n')
    ckpeaks = ttk.Checkbutton(self.tpick_paramframe, text='Tune on peak ?', variable=self.peaks, onvalue='y', offvalue = 'n')
    ckverbose = ttk.Checkbutton(self.tpick_paramframe, text='Verbose ?', variable=self.verbose, onvalue='y', offvalue = 'n')

    ckhedupdate.grid(row=0, column=12, ipadx=0, padx=5, pady=5, sticky='w')
    ckhedqcplot.grid(row=1, column=12, ipadx=0, padx=5, pady=5, sticky='w')
    ckhilbert.grid(row=0, column=13, ipadx=0, padx=5, pady=5, sticky='w')    
    ckpeaks.grid(row=1, column=13, ipadx=0, padx=5, pady=5, sticky='w')
    ckverbose.grid(row=0, column=14, ipadx=0, padx=5, pady=5, sticky='w')

    # Set some default values which can be manually updated
    self.entrytims.insert(0, '0')
    self.entrytime.insert(0, '1000')
    self.entrysta.insert(0, '25')
    self.entrylta.insert(0, '125')
    self.entryttrc.insert(0, '1')
    self.entrysubrat.insert(0, '1')
    self.entrytune.insert(0, '25')
    self.entrystab.insert(0,'.1')   
   
    # Insert default values that are read from data file
    self.entrytime.delete(0, tk.END)
    self.entrytime.insert(0, '%s'%(self.plotdata.shape[1]*(1000/self.fs)))
    
    # use to remove unwanted parameter grid when starting other grids
    self.tpick_menu_exists = self.tpick_paramframe.winfo_exists()

def stalta_params(self):
    ''' Get STA-LTA variables from gui and write to dictionary which
    can be passed around.
    '''
    stalta_params=dict()
    stalta_params['test_trace'] = int(self.entryttrc.get())
    stalta_params['stime_plot'] = float(self.entrytims.get())
    stalta_params['etime_plot'] = float(self.entrytime.get())
    stalta_params['shortwin'] = float(self.entrysta.get())
    stalta_params['longwin'] = float(self.entrylta.get())
    stalta_params['envelope'] = self.hilenv.get() # yes to use hilbert amplitude envelope
    stalta_params['subsamp_ratio'] = int(self.entrysubrat.get())
    stalta_params['tune'] = float(self.entrytune.get())
    stalta_params['stabilization'] = float(self.entrystab.get())
    stalta_params['verbose'] = self.verbose.get()
    stalta_params['qcplots'] = self.qcplt.get()
    stalta_params['peaks']=self.peaks.get()
    
    if stalta_params['peaks']=='y':
        stalta_params['picktype']='pk'
    else:
        stalta_params['picktype']='ipt'
    # variable to indicate tkinter useage
    stalta_params['tkinter_use']='y'
    
    return stalta_params


def pick_vsp(self,root):
    r''' Pick first breaks with variation of STA/LTA.
    pick_vsp calls recursive_sta_lta_py to make the inital picks

    :param a: Seismic Trace
    :type nsta: int
    :param nsta: Length of short time average window in samples
    :type nlta: int
    :param nlta: Length of long time average window in samples
    :rtype: NumPy :class:`~numpy.ndarray`
    :return: Characteristic function of recursive STA/LTA
    .. seealso:: [Withers1998]_ (p. 98) and [Trnkoczy2012]_
    
    :pickdata: full seismic matrix

    Returns : updated headers. Field times are overwritten

    Useage:
    pick_params = dict(
        test_trace = 13,
         shortwin = 20,
         longwin = 180,
         qcplots = 'y',
         stime_plot = 0,
         etime_plot = 3000,         
         subsamp_ratio = 1, ## increase is slower but get  fractional ms
         stabilization = .0001, # smaller if first break is relatively weak
         tune = 40,
         verbose='n', 
         table = 'n')
    picks_headers = tpick.pick_vsp(stack_bpf,zvsp_headers,fs,**pick_params) 

    There are 2 versions of pick_vsp, this one is necessary to pass self to the qc plots   
    '''
    import numpy as np
    import timepick.pickers as tpicks
    import timepick.pick_plot_gui as pickplot
    import iovsp.text_io as text_io

    pickdata=self.data
    raw_headers=self.headers[:,0:15] # assumes standard header file shape
    fs=self.fs

    # parameters from gui read from dictionary
    stalta_inputs=stalta_params(self)
    verbo=stalta_inputs['verbose']
    testtrace = stalta_inputs['test_trace']
    qcplots = stalta_inputs['qcplots'] 
       
    if (verbo == 'y')or(verbo == 'Y') :
        print("\u0332".join('\nPick VSP: Information :'))
        print ("\u0332".join('\nPick VSP Global Information :'))
        print (' raw_headers.shape :',raw_headers.shape)
        print (' fs input :',fs)

    shape = (pickdata.shape[0])
    pick = np.zeros(shape,dtype=np.float32)
    tunedpick = np.zeros(shape,dtype=np.float32)
     
    for i, trace in enumerate(pickdata[::1, :]):     
        pick[i], tunedpick[i],qcdata,pickparams,kwargs2 = tpicks.recursive_sta_lta_py(pickdata[i:i+1,:],i,fs,**stalta_inputs)

        if i==testtrace:
            if qcplots=='y':
                kwargs2['titletxt'] = 'Recursive STA\LTA QC Plots'
                pickplot.test_plot_canvas(self,root,pickparams,*qcdata, **kwargs2)

    new_headers = np.hstack((raw_headers, pick.reshape(-1,1),tunedpick.reshape(-1,1)))
    new_headers[:,16] = tunedpick # replace field picks with tuned picks    
   
    table=self.hedup.get()
    # print a table and copy STA LTA to field picks
    if (table=='y')or (table=='Y'):
        text_io.header_table(new_headers,numcols=17)
        new_headers[:,8]=new_headers[:,16]
    return new_headers

def tpicks(self,root):
    ''' Get automated picking paramaters
    1. set up tpick_param menu
    '''
    import timepick.pickers_gui as tpick

    tpick.sta_params_grid(self,root)

def tpicks_sta(self,root):
    ''' Get automated picking paramaters
    2. run STA LTA
    '''
    from timepick import pickers_gui
    from plotvsp import seisplots_gui
    # update header file with 2 extra columns for 
    # STA-LTA calculations
    self.headers = pickers_gui.pick_vsp(self,root)
    # use these new columns for the plots
    # make a plot with all data - no options
    self.plot_all = 'y'
    self.sta_pick = 'y'

    seisplots_gui.stack_plot(self)

    