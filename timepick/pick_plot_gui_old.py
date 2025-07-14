def sta_params_grid(self):

    # create a toplevel window (pop-out) to get STA LTA parametersand 
    # plot parameters
    # self.master is the root window in a form which can be passed to functions
    
    import tkinter as tk
    from tkinter import ttk
    import procvsp.utils_gui as utilg
    import timepick.pickers_gui as tpick

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
    labelfs.grid(row=0, column=15, ipadx=0,padx=5, pady=5, sticky='w')

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

    buttonttpick = ttk.Button(self.tpick_paramframe, text = "Run STA LTA", command = self.tpicks_sta)        
    buttonttpick.grid(row=0, column=14, columnspan=1, sticky='w',padx=10, pady=5 )

    buttontrplt = ttk.Button(self.tpick_paramframe, text = "Make Plot", command = self.plot_stack)        
    buttontrplt.grid(row=1, column=14, columnspan=1, sticky='w',padx=10, pady=5 )
    
    # create checkbuttons for plotting controls 
    
    self.hedup = tk.StringVar()
    self.qcplt = tk.StringVar()
    self.verbose = tk.StringVar()
       
    self.hedup.set('y')
    self.qcplt.set('n')
    self.verbose.set('n')

    ckhedupdate = ttk.Checkbutton(self.tpick_paramframe, text='Update Headers ?', variable=self.hedup, onvalue='y', offvalue = 'n')
    ckhedqcplot = ttk.Checkbutton(self.tpick_paramframe, text='QC plots ?', variable=self.qcplt, onvalue='y', offvalue = 'n')
    ckverbose = ttk.Checkbutton(self.tpick_paramframe, text='Verbose ?', variable=self.verbose, onvalue='y', offvalue = 'n')

    ckhedupdate.grid(row=0, column=12, ipadx=0, padx=5, pady=5, sticky='w')
    ckhedqcplot.grid(row=1, column=12, ipadx=0, padx=5, pady=5, sticky='w')    
    ckverbose.grid(row=0, column=13, ipadx=0, padx=5, pady=5, sticky='w')

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

def picker_multitrace(VSPdata, thead, fs, title_top, **kwargs):
    
    """Run the interactive picker on multiple traces
    
    Crossplot x (amplitude) and y (time). Add amplitude to receiver depth to 
    get trace deflection. Alternatively add amplitude to receiver number to get 
    trace deflection. Scaling in X direction (amplitude) is different in each 
    case
    
    Inputs
    VSPdata: VSP data array rows=traces
    thead: output from STA LTA time pickers
           Field TT column 8, STA/LTA TT column 15
           Peak tuned STA/LTA column 16 **used in code
    """
    import numpy as np    
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import scipy.signal as signal
    
    plt.rcParams.update({'font.size': 10})
    print("\u0332".join('\nPicker_Multi trace Global Information :'))
    print (' Input data shape :', VSPdata.shape, 'fs :',fs)
    
    rasterized = 'True' # True can speed up large data sets at expense of resolution


    scal = kwargs['Scalar']
    trange = kwargs['time_range']
    tracerange = kwargs['trace_range']
    pol= kwargs['polarity']
    skiplabel=kwargs['header_spacing']
    norm=kwargs['plot_norm']
    va=kwargs['var_area']
    tpi=kwargs['tr_per_in']
    #ips=kwargs['in_per_s']
    subsamp = kwargs['subsamp_ratio']

    ttime = thead[tracerange[0]:tracerange[1],16]
    ttime2 = thead[tracerange[0]:tracerange[1],15]
    # input transit time convert, to samples
    pick_in = thead[tracerange[0]:tracerange[1],16]*(fs//1000) # repicked by ?
    print (' thead.shape:', thead.shape)
    print (' pick_in :', pick_in)
    
    # trace deflection is based on sample value. Plots are spaced by receiver 
    # number or trace number. 
    # A scalar needs to be applied to make reasonable deflections, dependent 
    # on data amplitudes and plot spacing
    
    # scale the plot width by traces/inch, 12 being baseline plot width
    xscal = tpi/10 # 10 tpi standard
    #yscal = 5/ips # 5 ips standard    
    
    ssamp = int(trange[0]*fs/1000)
    esamp = int(trange[1]*fs/1000)
    
    print (' fs :',fs,'ssamp :', ssamp,' esamp :', esamp)
    VSPtrim = VSPdata[tracerange[0]:tracerange[1],ssamp:esamp]
    # make a local trace number vector
    trnum = np.arange(1,VSPtrim.shape[0]+1,1) # trace number data
    # output array for manipulations
    data2 = np.zeros(shape = (VSPtrim.shape[0], VSPtrim.shape[1]))    
    data1 = VSPtrim
    
    print ('Number of traces in plot :', VSPtrim.shape[0], 
           ' Number of samples per trace :', VSPtrim.shape[1])    
#    print (' thead shape :', thead.shape)

    # upsample the trace to get more precise pick        
    data1 = signal.resample(data1,data1.shape[1]*subsamp,axis=1, domain='time')
    fs_resamp = fs*subsamp

    print (' fs :',fs, ' fs_resamp :',fs_resamp)
    # upsampled time axis data
    ssamp_resamp = int(trange[0]*fs_resamp/1000)
    esamp_resamp = int(trange[1]*fs_resamp/1000)

    time = np.arange(ssamp_resamp*(1000/fs_resamp),esamp_resamp*(1000/fs_resamp),(1000/fs_resamp)) 
    print (' fs resampled :',fs_resamp)
    print (' time min , time max:',time.min(), time.max())    
    print (' data1 shape resampled :',data1.shape)


    if (norm == 'Y') or (norm =='y'):
        row_max = np.max(np.abs(data1), axis=1)
        where_0 = np.where(row_max == 0) # find traces of all 0s
        row_max[where_0] = 1 # set 0 to 1
        data2 = (data1 / row_max[:, np.newaxis])        
        datascaled = data2 * scal
#        print( 'row_max shape :', row_max.shape, ' row_max :', row_max)
    else:        
        datascaled = data1 * scal
                
    if (pol == 'r') or (pol =='R'):        
        datascaled = datascaled * -1
        
    print (' datascaled shape [0]',datascaled.shape[0], 
           ' datascaled shape [1]',datascaled.shape[1])
        
    skip =1
    
    # Get the value to add to the trace, receiver depth or trace number
    # When dscaler is addded to the trace, we can scale the plot by depth or trace number
    

    dscaler, pad = (trnum, 1)
    print (' dscaler.shape :', dscaler.shape, ' ttime.shape :', ttime.shape)
    print (' dscaler :', dscaler)
    print (' dscaler[0]-pad :',dscaler[0]-pad, ' dscaler[-1]+pad :',(dscaler[-1]+pad)*xscal)
    print (' np.min(trnum)-pad :',np.min(trnum)-pad, ' np.max(trnum)+pad :',(np.max(trnum)+pad)*xscal)

         
    ####### prepare header track data ##########################
    yaxis = trnum[::skiplabel]*0+.3    
    xaxis = trnum[::skiplabel]
    trace_num = trnum[::skiplabel]
    
    ######## set up colors, fonts ##############################
    markercolor = 'green'
    
    ####### set up plot axes ###################################         
    fig = plt.figure(figsize=(14,12))
    
    gs = gridspec.GridSpec(2, 1, height_ratios=[.2, 2], hspace = .02)    
    ax1a = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    # set up trace number header track
    for i, txt in enumerate(trace_num):
        ax1a.annotate(txt, (xaxis[i], yaxis[i]), rotation = 90)        
    #ax1a.spines['right'].set_visible(False)
    #ax1a.spines['top'].set_visible(False)
    #ax1a.spines['left'].set_visible(False)
    ax1a.set_ylabel('Trace number', rotation = 90 )
    ax1a.set_yticklabels([])    
    ax1a.set_xticks(trnum[:-1:1])
    ax1a.tick_params(axis = 'x', direction = 'in')
    ax1a.set_xticklabels([])
    ax1a.set_title('%s'%(title_top),fontsize=14)
    ax1a.set_xlim(np.min(trnum)-pad, (np.max(trnum) + pad)*xscal )

    ########## add the main track with traces ################    
    for i, trace in enumerate(datascaled[::skip, :]):
        x = trace + dscaler[i]
        line, = ax2.plot(x, time, 'k-', linewidth = .5, rasterized = rasterized)
        if (va=='y')or(va=='Y'):
            ax2.fill_betweenx(time, dscaler[i], x, where=(x > dscaler[i]), color='k', 
                              rasterized = rasterized)
        #ax2.plot(dscaler[i], ttime[i],marker='+', markersize=30, color=markercolor)
        #dots,=ax2.plot(dscaler[i], ttime[i],marker='+', markersize=30, color=markercolor)
#        ax2.plot(dscaler[i], ttime2[i],color='green', marker ="+",s=120  )
        ax2.set_xlim((dscaler[0]-pad), (dscaler[-1]+pad)*xscal)
        ax2.set_ylim(time.max(), time.min())
    ax2.plot(dscaler, ttime,marker='+', markersize=30, color=markercolor)
    #ax2.plot(dscaler, ttime2,marker='+', markersize=30, color='brown')
    dots,=ax2.plot(dscaler, ttime,marker='+', markersize=30, color=markercolor)
    ax2.set_ylabel(' Time (ms)')
    ax2.set_xticks([])              
    ax2.yaxis.grid()

    globaltime = ttime # prior to repick by lsta   
    globaldscaler = dscaler    

    def onclick(event):
        # :val Search for the globaldscaler value  (trace baseline) 
        # :tol Tolerance of the mouse click location
        # https://stackoverflow.com/questions/41022765/find-index-of-first-element-in-array-close-to-float-within-tolerance-with-numpy        
        # rtol gets bigger as the trace number gets bigger, so keep it small
        # atol is consistent around the trace number
        # the tolerance or epsilon = (atol + rtol * np.absolute(repick_resampx))
        
        print("\u0332".join('\nClick Event Information :'))

        if event.inaxes!= line.axes: return

        # create the new picks as global variables
        global repick_resampx, repick_resampy 

        # store the location of the button press event
        repick_resampx = event.xdata
        repick_resampy = event.ydata
        print (' repick_resampx:',repick_resampx,
               ' repick_resampy:',repick_resampy )
               
        # set the tolerances to reasonable values when the traces are spread by trace number       
        rtol = 1e-05 # relative tolerance
        atol = .2 # absolute tolerance - trace number +/- atol        
        i_array = np.where(np.isclose(globaldscaler, repick_resampx,atol=atol,rtol=rtol))
        my_i = None if (len(i_array[0])==0) else int(i_array[0])
        
        print (' i_array :',i_array,)# ' i_array[1] :',i_array[1])
        Index =np.array(i_array).reshape(-1)
        print (' Index.shape :', Index.shape,' Index :', Index, ' Index[-1] :', Index[-1])
        
        globaltime[Index[-1]] = repick_resampy

        # tolerance test
        epsilon = (atol + rtol * np.absolute(repick_resampx))
        print (' epsilon :', epsilon)

        # update the value of the pick
        dots.set_xdata(globaldscaler)
        dots.set_ydata(globaltime)
   
        # redraw only the dot
        dots.set_color('red') # changes marker color, adds line!
        ax2.draw_artist(dots, )
        
        # redraw the canvas with nothing changed except the pick
        fig.canvas.draw()
          

    # create a canvas for interactive plotting           
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
 
    #print (' pick_resampy :', repick_resampy)  
    #picktime = repick_resampy#*(1000/fs_resamp)
    #print (' picktime :', picktime)
    print (' globaltime :',globaltime)
    
    # update the tuned time header column
    thead[tracerange[0]:tracerange[1],16]  = globaltime 

    return thead

def test_plots(self,root,seis,fs, a, sta, lta,cft,index_max, index_tune, starts, ends, title_txt, prep_txt):

    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    
    import tkinter as tk    
    import tkinter as tk
    from tkinter import ttk
    import numpy as np
    
    # derivative in 100 sample sliding window
    seis_grad = np.gradient(seis,1000/fs)
    # derivative in 100 sample sliding window
    seis_grad_2nd = np.gradient(seis_grad,1000/fs)

#    slopew = slope(seis,50,-1,fs,)
#    print ( 'slopew.shape :', slopew.shape)

#    slopew = slopew[0,:]
    starttime=starts*(1000/fs)
    endtime=ends*(1000/fs)
    staltatime=index_max*1000/fs
    tunetime=index_tune*1000/fs

    time = np.arange(starttime,endtime,(1000/fs))
    
    print("\u0332".join('\nTest Plots Information :'))
    print (' time.shape : ',time.shape, ' seis.shape :',seis.shape, ' a.shape :', a.shape)
    
    ##############################################################################

    # create a toplevel window (pop-out) to plot response
    #
    staqcwindow = tk.Toplevel(root)
    staqcwindow.geometry('1300x800') # depends on fig. size 
    staqcwindow.wm_title("STA LTA QC")# double quotations!

    # this section is necessary to make scrollbars work
    # scrollbars need to be attached to a re-sizable Frame       
    staqcframe = ttk.Frame(staqcwindow)
    staqcframe.pack(fill="both", expand=True)            
    staqcframe.rowconfigure(0, weight=1)
    staqcframe.columnconfigure(0, weight=1)# this column will expand

    # create a tkinter canvas and scroll bars
    canvas_qcsta = tk.Canvas(staqcframe,borderwidth=1,relief=tk.RIDGE)#highlightthickness=1, highlightbackground="black")#plot_frame)
    canvas_qcsta.grid(row=0, column=0, sticky=tk.NSEW)#, padx = 0, pady =0,

    #Trace Plot Scroll bars
    xScrollbarsta = tk.Scrollbar(staqcframe, orient=tk.HORIZONTAL)
    yScrollbarsta = tk.Scrollbar(staqcframe)                
        
    xScrollbarsta.grid(row=1, column=0, sticky=tk.EW)
    yScrollbarsta.grid(row=0, column=1, sticky=tk.NS)

    canvas_qcsta.config(xscrollcommand=xScrollbarsta.set)#,xscrollincrement=0)
    canvas_qcsta.config(yscrollcommand=yScrollbarsta.set)#,yscrollincrement=0)
#    speed = 2

    xScrollbarsta.config(command=canvas_qcsta.xview)
    yScrollbarsta.config(command=canvas_qcsta.yview)
    
    ##############################################################################
    
    fig = Figure(figsize=(12,8)) # use Figure, not pyplot for tkinter to work
    gs = fig.add_gridspec(5, 1,height_ratios=[1,1,1,1,1])#, wspace = 0, hspace = .02)
    
    ax1 = fig.add_subplot(gs[0])

    ax1.plot(time,seis, 'k', label = 'seismic')
    ax1.scatter(staltatime, 0,color='red', marker ="+",s=120,label = 'STA/LTA' )
    ax1.scatter(tunetime, 0,color='green', marker ="+",s=120,label = 'STA/LTA tuned' )
    ax1.set_xlim(starttime, endtime)    
    ax1.legend()
    ax1.set_title('%s '%(title_txt))    

    ax2 = fig.add_subplot(gs[1])
    #plt.plot(time,seis_grad_2nd[starts:ends], 'k', label = 'seismic_gradient')
    ax2.plot(time,seis_grad_2nd, 'k', label = 'seismic_gradient')
    ax2.scatter(staltatime, 0,color='red', marker ="+",s=120,label = 'STA/LTA' )    
    ax2.scatter(tunetime, 0,color='green', marker ="+",s=120,label = 'STA/LTA tuned' )    
    ax2.set_xlim(starttime, endtime)    
    ax2.legend()

    ax3 = fig.add_subplot(gs[2])
    #plt.plot(time,a[starts:ends], 'k', label = 'a %s'%(prep_txt))
    ax3.plot(time,a, 'k', label = 'a %s'%(prep_txt))
    ax3.scatter(staltatime, 0,color='red', marker ="+",s=120,label = 'STA/LTA' )
    ax3.scatter(tunetime, 0,color='green', marker ="+",s=120,label = 'STA/LTA tuned' )    
    ax3.set_xlim(starttime, endtime)    
    ax3.legend()

    ax4 = fig.add_subplot(gs[3])
    #plt.plot(time,sta[starts:ends], 'k', label = 'sta')
    #plt.plot(time,lta[starts:ends], 'r', label = 'lta')     
    ax4.plot(time,sta, 'k', label = 'sta')
    ax4.plot(time,lta, 'r', label = 'lta')   
    ax4.scatter(staltatime, 0,color='red', marker ="+",s=120,label = 'STA/LTA' )
    ax4.scatter(tunetime, 0,color='green', marker ="+",s=120,label = 'STA/LTA tuned' )    
    ax4.set_xlim(starttime, endtime)    
    ax4.legend()

    ax5 = fig.add_subplot(gs[4])
    #plt.plot(time,cft[starts:ends], 'k', label = 'cft')
    ax5.plot(time,cft, 'k', label = 'cft')
    ax5.scatter(staltatime, 0,color='red', marker ="+",s=120,label = 'STA/LTA' )
    ax5.scatter(tunetime, 0,color='green', marker ="+",s=120,label = 'STA/LTA tuned' )    
    ax5.set_xlim(starttime, endtime)    
    ax5.legend()
#    plt.hlines([3.5, 0.5], 0, len(cft), color=['r', 'b'], linestyle='--')
#    plt.axis('tight')
#  remove whitespace at top and bottom of plot, adjust if a title is required
    fig.subplots_adjust( top=0.99, bottom=0.01)
    
    plt.show()

    figAgg = FigureCanvasTkAgg(fig, canvas_qcsta)
    mplCanvas = figAgg.get_tk_widget()
    canvas_qcsta.create_window(0, 0, window=mplCanvas, anchor=tk.NW)
    canvas_qcsta.config(scrollregion=canvas_qcsta.bbox(tk.constants.ALL))
    
    tk.Misc.lift(staqcframe)