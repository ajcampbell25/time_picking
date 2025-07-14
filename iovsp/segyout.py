def write_segyio(data, texthead,tracehead,fs, name, **kwargs):
    
    '''Write a segy file using segyio routines
    
    Number of traces can come from binary header or trace headers or
    the shape of the data file

    We use the matrix shape in this example to figure out the number of 
    traces and the number of samples per trace
        
    data - 2d numpy array of VSP data [traces, samples]
    tracehead - 2d array of trace headers [traces, header values]
    fs - sample rate in Hz
    name - output file path
    
    1. Arbitrary byte positions not accessible for headers

    '''    
    import segyio

    ################# print some basic information ###############

    print("\u0332".join('\nWrite Segyio Stats :'))    
    print ('Data shape [0] :', data.shape[0],'Data shape [1] :', data.shape[1])    
    print('Trace header shape', tracehead.shape)    
    print ('time header :', tracehead[0:2,:])
    
    ################### geometry scalars  ########################

    scalcoord = 10
    scaldep = 10

    ################## byte locations and scalars ################    
    timbit = kwargs['timbit']
    timscal = kwargs['timscal']
    sxbit=kwargs['sxbit']
    sybit=kwargs['sybit']
    szbit=kwargs['szbit']
    rxbit=kwargs['rxbit']
    rybit=kwargs['rybit']
    rzbit=kwargs['rzbit']

    
    ################# read trace headers ##########################

    FFID = tracehead[:,14].astype(int)
    SRC = tracehead[:,13].astype(int)
    MD = (tracehead[:,1]*scaldep).astype(int)
    TVD = (tracehead[:,2]*scaldep).astype(int)
    RcvX = (tracehead[:,3]*scalcoord).astype(int)
    RcvY = (tracehead[:,4]*scalcoord).astype(int)
    SrcX = (tracehead[:,5]*scalcoord).astype(int)
    SrcY = (tracehead[:,6]*scalcoord).astype(int)
    SrcZ = (tracehead[:,7]*scalcoord).astype(int)
    TVD_SRD = (tracehead[:,9]*scaldep).astype(int)
    TVD_Src = (tracehead[:,10]*scaldep).astype(int)
    SrcZ_SRD = (tracehead[:,11]*scaldep).astype(int)
    Tobs = (tracehead[:,8]*timscal).astype(int)

    print ('MD shape :', MD.shape, ' MD dtype :' , MD.dtype)
    print (' MD [1:10] :', MD[0:10], ' TVD[0:10] :',TVD[0:10])
    
    ################# create text header ##########################
    ''' 
    text_header = {1: 'Synthetic Walkaway VSP', 
    2: 'SEG-Y read-write tests',
    3: 'VSProwess',
    4: ' ',
    5: 'WELL: Allan 1',
    6: ' ',
    7: 'Processed by VSP Consultants',
    8: ' ' ,
    9: 'Reference Elevation:  ft',
    10: ' ',
    11: 'Geophone component: 025-028 Z=1, X=2, Y=3' ,
    12: 'Source Easting:  073-076     Source Elev: 045-048' ,
    13: 'Source Northing: 077-080' ,
    14: 'Receiver Easting: 081-084    Measured Depth: 037-040' ,
    15: 'Receiver Northing: 085-088   Vertical Depth: 041-044' ,
    16: '  ',
    17: 'Uncorrected Pick Time: 107-108' ,
    18: 'TWO-Way Time : 109-110' ,
    19: ' ',
    20: ' ',
    21: 'Units  = Survey Feet' ,
    22: 'Wellhead Easting (ft): 0' ,
    23: 'Wellhead Northing (ft): 0' ,
    24: ' ',
    25: ' ****Processing Steps: ***********' ,
    26: 'Original file from VSProwess' ,
    27: 'Read and write with segyio',
    28: ' ' ,
    29: ' ',
    30: ' ', 
    31: ' ',
    32: ' ',
    33: ' ',
    34: ' ',
    35: ' ',
    36: ' ',
    37: ' ',
    38: 'VSProwess processing system'  }
    '''    
    txt_hed = segyio.tools.create_text_header(texthead)
    
    ################# set some global parameters ###################

    spec = segyio.spec()
    spec.samples = list(range(data.shape[1]))
    spec.sorting = 0 # 0 for unstructured ie. not a cube of data
    spec.format  = 5 #1 = IBM float, 5 = IEEE float
    spec.tracecount = data.shape[0]
    
    print( ' spec.tracecount :', spec.tracecount)
    
    ################# create output file  #########################
    
    with segyio.create(name, spec) as fout:
        fout.trace = data                      # populate traces
        fout.text[0] = txt_hed
        fout.bin.update(hns=len(spec.samples)) # update ample count binary header, su style keys 
        fout.bin.update(format=5)      #1 = IBM float, 5 = IEEE float        
        fout.bin.update(mfeet=2)        
        fout.bin.update(rev=1)                                       
#        fout.bin.update{segyio.BinField.MeasurementSystem: 2} #alternate methods
#        fout.bin.update{3255: 2}      # all binary keys use 3200 + 55 etc for byte positions

    ############## create output trace headers  ####################
        tr=0
        for i, trace in enumerate(data):
#            fout.header[tr] = {segyio.tracefield.TraceField.ReceiverGroupElevation: MD[tr]}
            fout.header[tr] = {1: i+1}    
            fout.header[tr] = {115: data.shape[1]}
            fout.header[tr] = {117:int((1/fs)*1000000)}
            fout.header[tr] = {9: FFID[tr]}
            fout.header[tr] = {13: SRC[tr]}
            fout.header[tr] = {rzbit: MD[tr]}
            fout.header[tr] = {41: TVD[tr]}
            fout.header[tr] = {69: -1*scaldep}
            fout.header[tr] = {71: -1*scalcoord}
            fout.header[tr] = {sxbit: SrcX[tr]}
            fout.header[tr] = {sybit: SrcY[tr]}
            fout.header[tr] = {szbit: SrcZ[tr]}
            fout.header[tr] = {rxbit: RcvX[tr]}
            fout.header[tr] = {rybit: RcvY[tr]}            
            fout.header[tr] = {timbit: int(Tobs[tr]/timscal)}
            fout.header[tr] = {197: Tobs[tr]}
                    
            tr +=1


text_header = {1: 'Generic ZVSP', 
    2: 'Time Pick and SEG-Y read-write tests in jupyter',
    3: 'VSProwess',
    4: ' ',
    5: 'WELL: Allan 1',
    6: ' ',
    7: 'Processed by VSP Consultants',
    8: ' ' ,
    9: 'Reference Elevation:  ft',
    10: ' ',
    11: 'Geophone component: 025-028 Z=1, X=2, Y=3' ,
    12: 'Source Easting:  073-076     Source Elev: 045-048' ,
    13: 'Source Northing: 077-080' ,
    14: 'Receiver Easting: 081-084    Measured Depth: 037-040' ,
    15: 'Receiver Northing: 085-088   Vertical Depth: 041-044' ,
    16: '  ',
    17: 'Uncorrected Pick Time: 107-108' ,
    18: 'TWO-Way Time : 109-110' ,
    19: ' ',
    20: ' ',
    21: 'Units  = Survey Feet' ,
    22: 'Wellhead Easting (ft): 0' ,
    23: 'Wellhead Northing (ft): 0' ,
    24: ' ',
    25: ' ****Processing Steps: ***********' ,
    26: 'Original file from VSProwess' ,
    27: 'Read and write with segyio',
    28: ' ' ,
    29: ' ',
    30: ' ', 
    31: ' ',
    32: ' ',
    33: ' ',
    34: ' ',
    35: ' ',
    36: ' ',
    37: ' ',
    38: 'VSProwess processing system'  }

