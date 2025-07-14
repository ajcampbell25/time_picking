def sta_lta(pickdata, fs, **kwargs):# nstatime, nltatime, fs, testtrace):
    """ Classic STA LTA

    From https://github.com/jwellik/STALTA_tuner/blob/master/stalta_tuner/trigger.py
    
    :copyright:
    The ObsPy Development Team (devs@obspy.org)
    :license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html) 
    
    Computes the standard STA/LTA from a given input array a. The length of
    the STA is given by nsta in samples, respectively is the length of the
    LTA given by nlta in samples. Written in Python.
    .. note::
        There exists a faster version of this trigger wrapped in C
        called :func:`~obspy.signal.trigger.classic_sta_lta` in this module!
    :type a: NumPy :class:`~numpy.ndarray`
    :param a: Seismic Trace
    :type nsta: int
    :param nsta: Length of short time average window in samples
    :type nlta: int
    :param nlta: Length of long time average window in samples
    :rtype: NumPy :class:`~numpy.ndarray`
    :return: Characteristic function of classic STA/LTA
    
    :testtrace = if -1 use all traces else use select a trace

    """
    import numpy as np
    #from scipy import fftpack    
    import scipy.signal as signal
    import timepick.pick_plot as pickplot
    
    print("\u0332".join('\nclassic sta-lta: Information :'))
    
    starts = kwargs['stime_plot']
    ends = kwargs['etime_plot']
    testtrace = kwargs['test_trace']
    plots = kwargs['qcplots']
    nstatime = kwargs['shortwin']
    nltatime = kwargs['longwin']
    subsamp = kwargs['subsamp_ratio']
    tunetime = kwargs['tune']
    stab = kwargs['stabilization']
    verbo = kwargs['verbose']
    
    stab=.01 # add a stabiliztion factor to denominator of STA/LTA  

    starts1 = int(starts*(fs//1000))# plotting start in samples
    ends1 = int(ends*(fs//1000))

    select_data=pickdata[testtrace:testtrace+1,:].reshape(-1)
    
    a_resamp = signal.resample(select_data,select_data.shape[0]*subsamp,domain='time')
    print (' select_data.shape :', select_data.shape,' a_resamp.shape:', a_resamp.shape)    

    #a_norm=a/np.max(a)      
    #a_resamp = signal.resample(a_norm,a_norm.shape[0]*subsamp,domain='time')
    fs_resamp = fs*subsamp
    starts = int(starts1*(fs_resamp//1000))# plotting start in samples
    ends = int(ends1*(fs_resamp//1000))
    a_select= (a_resamp[starts:ends])#.T).reshape(-1)
    # square amplitudes to get energy 
    a = a_select**2
    print (' a.shape :', a.shape,' starts:', starts,' ends :', ends)    
    # some example pick on the hilbert envelope
    #    analytic_signal = signal.hilbert(a_resamp) 
    #    amp_env = np.abs(analytic_signal)    

    # convert time windows to samples (milliseconds) 
    nsta = int(nstatime*(fs_resamp//1000)) # sta time window in samples
    nlta = int(nltatime*(fs_resamp//1000))
    
    # The cumulative sum can be exploited to calculate a moving average (the
    # cumsum function is quite efficient)
    sta = np.cumsum(a)#,dtype=np.float64)
#    sta = np.cumsum(amp_env) # hilbert envelope picking

    # Convert to float
    sta = np.require(sta, dtype=float)

    # Copy for LTA
    lta = sta.copy()
#    print (' sta :', sta, 'lta :',lta)

    # Compute the STA and the LTA
    sta[nsta:] = sta[nsta:] - sta[:-nsta]
    sta /= nsta
    
    lta[nlta:] = lta[nlta:] - lta[:-nlta]
    print (' sta.shape :', sta.shape,' lta.shape :', lta.shape)
    lta /= nlta

    # Pad zeros
#    sta[:nlta//2 - 1] = 0

    # Avoid division by zero by setting zero values to tiny float
    dtiny = np.finfo(0.0).tiny
    idx = lta < dtiny
    lta[idx] = dtiny

    cft = sta / (lta +stab*lta.max())

    index_max = int(np.argmax(cft) + starts)    
    pick_time = index_max*(1000/fs_resamp) 
    index_tune, tuned_time = peak_tune(a_resamp, fs_resamp, index_max, tunetime,verbo)
    
    print (' testtrace : ',testtrace,'  index_max :', index_max, ' fs_resamp :',fs_resamp,
               ' pick_time:',pick_time, )
    if testtrace!=-1:
        print (' classic index_max :', index_max, ' pick_time:',pick_time) 
        
    if (plots=='y')or(plots=='Y'):
        titletxt = 'Classic STA\LTA'
        preptxt = 'resamp**2'
        pickplot.test_plots(a_select,fs_resamp, a, sta, lta,cft,index_max,index_tune, starts, ends,
               titletxt, preptxt)
    
    return pick_time

def recursive_sta_lta_py(pickdata,trace, fs,**kwargs):    

    """    Recursive STA/LTA written in Python.

    This routine is modified from the obspy version:
    https://docs.obspy.org/packages/autogen/obspy.signal.trigger.recursive_sta_lta_py.html

    Some modifications were guided by The Crewes routine Picker.m. The addition of the 
    stabilization factor helps where direct arrival is relatively weak.

    .. note::
        There exists a faster version of this trigger wrapped in C
        called :func:`~obspy.signal.trigger.recursive_sta_lta` in this module!
    :type a: NumPy :class:`~numpy.ndarray`
    :param a: Seismic Trace
    :type nsta: int
    :param nsta: Length of short time average window in samples
    :type nlta: int
    :param nlta: Length of long time average window in samples
    :rtype: NumPy :class:`~numpy.ndarray`
    :trace: iteration from pick_vsp for testtrace plotting 
    :stab: Stabilization factor (from Crewes Picker.m)
    :return: Characteristic function of recursive STA/LTA
    .. seealso:: [Withers1998]_ (p. 98) and [Trnkoczy2012]_
    
    :pickdata: full seismic matrix
    :trace: iteration number from pick_vsp, if it matches testtrace, 
            make a plot
    testtrace = output a test plot for this trace
    Useage:
    picking on a selection of traces, using a for loop. Kwargs come from pick_vsp
    for i, trace in enumerate(pickdata[::1, :]):
        pick[i], tunedpick[i] = tpicks.recursive_sta_lta_py(pickdata[i:i+1,:],i, fs, **kwargs) 
    """
    import numpy as np
    import scipy.signal as signal
    import timepick.pickers as tpicks
    import timepick.pick_plot as pickplot
                            
    starts = kwargs['stime_plot']
    ends = kwargs['etime_plot']
    testtrace = kwargs['test_trace']
    plots = kwargs['qcplots']
    nstatime = kwargs['shortwin'] # convert to nsta (samples)
    nltatime = kwargs['longwin']  # convert to nlta (samples)
    hilenv = kwargs['envelope'] # yes to use hilbert amplitude envelope
    subsamp = kwargs['subsamp_ratio']
    tunetime = kwargs['tune']
    stab = kwargs['stabilization']
    verbo=kwargs['verbose'] # print a lot of qc info to the terminal
    
    #stab=.001 # add a stabiliztion factor to denominator of STA/LTA
              # stab should be smaller than the other methods

    starts1 = int(starts*(fs//1000))# plotting start in samples
    ends1 = int(ends*(fs//1000))
 
    resamp1 = signal.resample(pickdata,pickdata.shape[1]*subsamp,axis=1, domain='time')   
    fs_resamp = fs*subsamp

    # start calculations   
    nsta = nstatime*(fs_resamp//1000) # convert time to number of samples
    nlta = nltatime*(fs_resamp//1000)    
    starts = int(starts1*(fs_resamp//1000)) # plotting start in samples
    ends = int(ends1*(fs_resamp//1000))
    a_plot= (resamp1[:,starts:ends].T).reshape(-1)
    a= (resamp1[:,starts:ends].T).reshape(-1)
    # reset start and end to resampled values
    kwargs['starts']=starts
    kwargs['ends']=ends

    ### some examples use envelope instead of trace
    kwargs['labelcft'] = 'cft from trace amp.'
    kwargs['preptxt'] = 'seismic squared'
    analytic_signal = signal.hilbert(a)
    amp_env = np.abs(analytic_signal)
    if (hilenv =='y')or(hilenv =='y'):
        a=amp_env
        kwargs['preptxt'] = 'hilbert envelope'
        kwargs['labelcft'] = 'cft from hilbert env.'
    
    if (verbo=='y')or(verbo=='Y'):
        print("\u0332".join('\nRecursive sta-lta: Information :'))
        print (' testtrace :',testtrace, ' pickdata.shape :',pickdata.shape)
        print (' plots y or n :',plots)
        print (' resamp.shape :',resamp1.shape, 'resampled limits a.shape : ',a.shape)
        print (' fs : ',fs,' starts :',starts1, ' ends :',ends1)
    
    ndat = len(a)
    # compute the short time average (STA) and long time average (LTA)
    # given by Evans and Allen
    csta = 1. / nsta
    clta = 1. / nlta
    sta3 = 0.
    lta3 = 0 #np.finfo(0.0).tiny  # avoid zero division
    a3 = np.square(a)
    # do not square if using hibert envelope
    if (hilenv =='y')or(hilenv =='y'):
        a3=amp_env

    sta4 = np.zeros(ndat, dtype=np.float64)
    lta4 = np.zeros(ndat, dtype=np.float64)

    icsta = 1 - csta
    iclta = 1 - clta
    for i in range(1, ndat):
        sta3 = csta * a3[i] + icsta * sta3
        lta3 = clta * a3[i] + iclta * lta3
        sta4[i] = sta3
        lta4[i] = lta3
    charfct2 = sta4 / (lta4 +(stab*lta4.max())) #from Crewes Picker.m        

    index_max = int(np.argmax(charfct2)+starts)
    pick_time = index_max*(1000/fs_resamp)
    
    if (verbo=='y')or(verbo=='Y'):
        print (' ndat : ',ndat, ' sta4.shape :',sta4.shape,' charfct2.shape :',charfct2.shape)
        print (' trace : ',trace,' recursive index_max :', index_max, ' fs_resamp :',fs_resamp,
               ' pick_time:',pick_time, )
    
    # derivative of trace amplitude
    trc_grad1 = np.gradient(a_plot,1000/fs_resamp)
    kwargs['labelgradtrc1'] = 'trace gradient'
 
    # 2nd derivative of trace amplitude
    trc_grad2 = np.gradient(trc_grad1,1000/fs_resamp)
    kwargs['labelgradtrc2'] = 'trace 2nd gradient' 

    # collect gradient into list for tuning   
    tunedata=(trc_grad1,trc_grad2)

    ### tune to peak value in window around pick_time
    index_tune, tuned_time = tpicks.peak_tune(a_plot, fs_resamp, index_max, tunetime,*tunedata, **kwargs)

    # collect fs and indices for plotting
    outparams=[fs_resamp,index_max,index_tune]

    testdata=(a_plot,a3, sta4, lta4,charfct2,trc_grad1,trc_grad2)

    return pick_time, tuned_time, testdata, outparams, kwargs

def peak_tune( seis,fs_resamp, index_max, tunetime,*td2, **kwargs):
    ''' find the index of the peak amplitude in a window.
    Window is centerd around index_max.
    Inputs
    : a: seismic trace
    : fs_resamp: sample rate
    : index_max: sample index around which window is applied 
                 Usually time pick index around which tuning to peak is desired
    : tunetime: window length
    : ptype: tune to peak or inflection point tangent nearest to STA LTA pick
           : ipt or zc is inflection point tangent or zero-crossing
           : anything else is tune to peak only
    Returns
    : index_tune STA-LTA pick tuned to largest peak in window, or to 
      first inflection point tangent before first peak
    : tuned_time STA-LTA pick tuned to largest peak,, or to 
      first inflection point tangent before first peak, converted to time
    '''
    import numpy as np

    grad1=td2[0] # seismic trace gradient
    grad2=td2[1] # gradient of gradient

    verbo=kwargs['verbose'] # print lots of stuff
    ptype = kwargs['picktype'] 

    #convert the window time to samples
    tunewin = int(tunetime*(fs_resamp//1000))    
    
    # prevent a crash if pick is close to 0
    winstart=(index_max - (tunewin//2))
    winend=(index_max + (tunewin//2))
    if verbo =='y':
        print("\u0332".join('\nPeak Tuner: Information :'))
        print (' tunewin :',tunewin,' winstart :',winstart,' winend :',winend)
        print (' seis.shape :',seis.shape)
    if winstart<0:
        winstart=0
        winend=0+tunewin

    # tune on peak of seismic - vibroseis sources
    index_tune = np.argmax(seis[winstart:winend])
    # get index relative to zero time sample
    index_tune += (index_max - (tunewin//2))
    # convert from sample number to time
    tuned_time = index_tune*(1000/fs_resamp)
    # inflection point tangent or zero-crossing picks, 
    # min phase source
    if (ptype =='ipt')or(ptype =='zc'):
        # tune on min of seismic gradient - equivalent to  SLB
        # inflection point tangent picking for airgun sources
        # zc should be zero-crossing, not done yet
        from scipy.signal import argrelmin
        # set a window relative to previously calculate peak amplitude
        winstart_min=(index_tune - (tunewin//2))
        # find indices of relative min in window
        if ptype =='ipt':
            min_index_tune = np.array(argrelmin(grad1[winstart_min:index_tune]))
        else:
            min_index_tune = np.array(argrelmin(grad2[winstart_min:index_tune]))
            # sometimes no minima are found - need to research
        if min_index_tune.size==0:
            min_index_tune = index_tune
        if min_index_tune.size>1:
            min_index_tune = min_index_tune[0,0]  
        # get indices of all minima in window relative to zero time sample
        min_index_tune += (index_tune - (tunewin//2))
        if verbo=='y':
            print (' min_index_tune', min_index_tune)
            print (' min_index_tune.shape', min_index_tune.shape)
            print (' min_index_tune.size', min_index_tune.size)
        # convert from sample number to time,
        # use only the first of the minima
        #if min_index_tune.size>1:
        #    tuned_time = min_index_tune[0,0]*(1000/fs_resamp)
        #else:
        tuned_time = min_index_tune*(1000/fs_resamp)
        index_tune = min_index_tune

    if verbo =='y':
        print('STA LTA_time :', index_max*(1000/fs_resamp))
        print('tuned_time :',tuned_time)

    return index_tune, tuned_time


def pick_vsp(pickdata, raw_headers, fs, **kwargs):
    ''' Pick first breaks with variation of STA/LTA.
    pick_vsp calls recursive_sta_lta_py to make the inital picks.

    Do not use this version with tkinter - handles the test plots differently.

    :param a: Seismic Trace
    :type nsta: int
    :param nsta: Length of short time average window in samples
    :type nlta: int
    :param nlta: Length of long time average window in samples
    :rtype: NumPy :class:`~numpy.ndarray`
    :return: Characteristic function of recursive STA/LTA
    .. seealso:: [Withers1998]_ (p. 98) and [Trnkoczy2012]_
    
    :pickdata: full seismic matrix

    Useage:
    pick_params = dict(
        test_trace = 13,
         shortwin = 20,
         longwin = 180,
         qcplots = 'y',
         stime_plot = 0,
         etime_plot = 3000,
         envelope='y'# use Hilbert amplitude envelope         
         subsamp_ratio = 1, ## increase is slower but get  fractional ms
         stabilization = .0001, # smaller if first break is relatively weak
         tune = 40,
         verbose='n', 
         table = 'n')
    picks_headers = tpick.pick_vsp(stack_bpf,zvsp_headers,fs,**pick_params)

    There are 2 versions of pick_vsp, this one is the basic, all-purpose one    
    '''
    import numpy as np
    import timepick.pickers as tpicks
    import timepick.pick_plot as pickplot
    import iovsp.text_io as text_io

    # print some error tracking info
    verbo=kwargs['verbose']
    #print a table of updated headers to screen    
    table=kwargs['table']
    # variable to indicate tkinter useage
    kwargs['tkinter_use']='n'

    testtrace = kwargs['test_trace']
    qcplots = kwargs['qcplots']

    shape = (pickdata.shape[0])
    shape2 = (pickdata.shape)
    pick = np.zeros(shape,dtype=np.float32)
    tunedpick = np.zeros(shape,dtype=np.float32)

    if (verbo == 'y')or(verbo == 'Y') :
        print("\u0332".join('\nPick VSP: Information :'))
        print (' kwargs etime : ',kwargs['etime_plot'])
        print (' table y or n : ',kwargs['table'])
    # run STA/LTA trace by trace   
    for i, trace in enumerate(pickdata[::1, :]):        
        pick[i], tunedpick[i], qcdata,pickparams,kwargs2 = tpicks.recursive_sta_lta_py(pickdata[i:i+1,:],i, fs, **kwargs) # carefulwith input array shape
        # plot the test traces
        if i==testtrace:
            if (qcplots=='y')or(qcplots=='Y'):
                kwargs2['titletxt'] = 'Recursive STA\LTA QC Plots'
                fig,ax = pickplot.test_plots_guts(pickparams,*qcdata, **kwargs2)
                fig.canvas.draw_idle()  # should be faster then plt.show but in reality is not noticeable

    # Add 2 new headers to input heeader file
    # 1. STA/LTA pick time 2. Tuned pick time    
    new_headers = np.hstack((raw_headers, pick.reshape(-1,1),tunedpick.reshape(-1,1)))
    if (verbo == 'y')or(verbo == 'Y') :
        print (' pick[10] :', pick[10],' tunedpick[10] :',tunedpick[10])
    #print a table of updated headers to screen    
    if (table=='y')or (table=='Y'):
        text_io.header_table(new_headers,numcols=16)

    return new_headers

