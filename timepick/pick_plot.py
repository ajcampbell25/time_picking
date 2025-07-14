def picker_multitrace(VSPdata, thead, fs, **kwargs):
    
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
    from matplotlib.figure import Figure
    
    
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

    first=tracerange[0]
    last=tracerange[1]
    
    # if STA LTA is run, 2 extra time columns are appended to thead
    if thead.shape[1] !=15:
        ttime = thead[first:last,16]
    else:
        ttime = thead[first:last,8]

    tk_used = kwargs['tkinter_use']
    print (' tkinter used? :',tk_used)

    # scale the plot width by traces/inch, 10 being baseline plot width
    xscal = tpi/10 # 10 tpi standard

    ssamp = int(trange[0]*fs/1000)
    esamp = int(trange[1]*fs/1000)
    
    print (' fs :',fs,'ssamp :', ssamp,' esamp :', esamp)
    VSPtrim = VSPdata[first:last,ssamp:esamp]
    # make a local trace number vector
    trnum = np.arange(1,VSPtrim.shape[0]+1,1) # trace number data
    # output array for manipulations
    data2 = np.zeros(shape = (VSPtrim.shape[0], VSPtrim.shape[1]))    
    data1 = VSPtrim
    
    print ('Number of traces in plot :', VSPtrim.shape[0], 
           ' Number of samples per trace :', VSPtrim.shape[1])    
    print (' thead shape :', thead.shape)

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
    markercolor = 'red'
    
    ####### set up plot axes ###################################
    # 
    fig = plt.figure(figsize=(14,12))    
    if tk_used=='y':
        fig = Figure(figsize=(14,12))    
        gs = fig.add_gridspec(2, 1, height_ratios=[.2, 2], hspace = .02)
        axm=gs.subplots(sharex=True)    

    else:
        gs = fig.add_gridspec(2, 1,height_ratios=[.1,1], wspace = 0, hspace = .02)
        axm = gs.subplots(sharex=True)    


    # set up trace number header track
    for i, txt in enumerate(trace_num):
        axm[0].annotate(txt, (xaxis[i], yaxis[i]), rotation = 90)        
    axm[0].set_ylabel('Trace number', rotation = 90 )
    axm[0].set_yticklabels([])    
    axm[0].set_xticks(trnum[:-1:1])
    axm[0].tick_params(axis = 'x', direction = 'in')
    axm[0].set_xticklabels([])
    axm[0].set_title('Manual Time Picking',fontsize=14)
    axm[0].set_xlim(np.min(trnum)-pad, (np.max(trnum) + pad)*xscal )

    ########## add the main track with traces ################    
    for i, trace in enumerate(datascaled[::1, :]):
        x = trace + dscaler[i]
        line, = axm[1].plot(x, time, 'k-', linewidth = .5, rasterized = rasterized)
        if (va=='y')or(va=='Y'):
            axm[1].fill_betweenx(time, dscaler[i], x, where=(x > dscaler[i]), color='k', 
                              rasterized = rasterized)
        axm[1].set_xlim((dscaler[0]-pad), (dscaler[-1]+pad)*xscal)
        axm[1].set_ylim(time.max(), time.min())
    # plot ttime twice, one version stays on canvas, the other gets updated
    # by mouse clicks
    axm[1].plot(dscaler, ttime,marker='_', markersize=20, color='red')
    dots,=axm[1].plot(dscaler, ttime,marker='+', markersize=20)#, color=markercolor)
    axm[1].set_ylabel(' Time (ms)')
    axm[1].set_xticks([])              
    axm[1].yaxis.grid()

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
        dots.set_color('green') # changes marker color, adds line!
        axm[1].draw_artist(dots, )
        
        # redraw the canvas with nothing changed except the pick
        fig.canvas.draw()
        if thead.shape[1]!=15:
            thead[first:last,16]  = globaltime
            # copy STA LTA manually updated picks to field pick
            #if update_h=='y':
            #    thead[first:last,8]  = globaltime
        else:
            # field picks updated after manual pick - no STA LTA
            thead[first:last,8]  = globaltime
        print (' thead[first:last,8] :',thead[first:last,8])          
    
    # create a canvas for interactive plotting
    #plt.show()           
    fig.canvas.mpl_connect('button_press_event', onclick)

    if tk_used=='n':
        plt.show()
        return thead
    else:
        return fig,thead
    
def test_plots_guts(pick_params,*args,**kwargs):
    ''' separate canvas from plotting for re-usability of plotting code
    '''
    import numpy as np    
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.widgets import MultiCursor
    import numpy as np

    starts = kwargs['starts']
    ends = kwargs['ends']
    title_txt = kwargs['titletxt']
    prep_txt = kwargs['preptxt']
    #labelgrad= kwargs['labelgrad1']
    #labelgrad2 = kwargs['labelgrad2']
    labelgradtrc = kwargs['labelgradtrc1']
    labelgradtrc2 = kwargs['labelgradtrc2']
    labelcft = kwargs['labelcft']
    tk_used = kwargs['tkinter_use']

    seis = args[0]
    a3 = args[1]
    sta = args[2]
    lta = args[3]
    cft = args[4]
    trc_grad1 = args[5]
    trc_grad2 = args[6]

    fs=pick_params[0]
    index_max=pick_params[1]
    index_tune=pick_params[2]

    #create time vector
    starttime=int(starts)*1000/fs
    endtime=int(ends)*1000/fs
    time = np.arange(starttime,endtime,(1000/fs))

    #convert from sample number to time
    staltatime=index_max*1000/fs
    tunetime=index_tune*1000/fs

    print("\u0332".join('\nTest Plots Information :'))
    print (' fs :',fs)
    print (' time.shape : ',time.shape, ' seis.shape :',seis.shape, ' a3.shape :', a3.shape)
    print (' starts :',starts,' ends :',ends)
    print (' starttime :',starttime,' endtime :',endtime)
    print(' tk used  :',tk_used)

    rasterized='True'
    numplots = len(args)-1 #we plot 2 traces in track3, making for 1 less plot

    
    fig = plt.figure(figsize=(12,8))# use figure when using pyplot, not for tkinter
    if tk_used=='y':
        fig = Figure(figsize=(12,12)) # use Figure, not pyplot for tkinter to work

    gs = fig.add_gridspec(numplots, 1, wspace = 0, hspace = .4)
    axs = gs.subplots(sharex=True)
    fig.suptitle('%s '%(title_txt))

    for ax in axs:   
        ax.scatter(staltatime, 0,color='red', marker ="+",s=120,label = 'STA/LTA' )
        ax.scatter(tunetime, 0,color='green', marker ="+",s=120,label = 'STA/LTA tuned' )
        ax.grid(axis='x')
        ax.set_xlim(starttime, endtime)
          
    axs[0].plot(time,seis, 'k', label = 'seismic',rasterized=rasterized) 
    axs[0].set_ylim(np.min(seis),np.max(seis))

    axs[1].plot(time,a3, 'k', label = prep_txt,rasterized=rasterized)
    axs[1].set_ylim(np.min(a3),np.max(a3))

    axs[2].plot(time,sta, 'k', label = 'sta')
    axs[2].plot(time,lta, 'r', label = 'lta')   
    axs[2].set_ylim(np.min(sta),np.max(sta))

    axs[3].plot(time,cft, 'k', label = labelcft,rasterized=rasterized)
    axs[3].set_ylim(np.min(cft),np.max(cft))

    axs[4].plot(time,trc_grad1, 'k', label = labelgradtrc,rasterized=rasterized)
    axs[4].set_ylim(np.min(trc_grad1),np.max(trc_grad1))

    axs[5].plot(time,trc_grad2, 'k', label = labelgradtrc2,rasterized=rasterized)
    axs[5].set_ylim(np.min(trc_grad2),np.max(trc_grad2))
    for ax in axs:   
        ax.legend(loc='right')

    # Enable interactive cursors - *****unable to get this to work on tkinter canvas*****
    # use fig.get_axes() to gather all the plot axes instead of (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8)
    # use plt.cursor seems to work better for running in gui or jupyter
    if tk_used=='n':
        plt.cursor = MultiCursor(fig.canvas,fig.get_axes() , color='r',lw=0.5, horizOn=True, 
                                  vertOn=True, useblit=True)
    # run plt.show on the returned figure - useful for making this run in 
    # tkinter and other places    
    return fig,axs