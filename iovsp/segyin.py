
def readsegyio3(inputfile,**kwargs):
    
    """ Load a segy file using segyio from Equinor

    Number of traces can come from binary header or trace headers or
    the shape of the data file
        
    1. open the segy file.
    2. The numpy data file is created in the first with loop.
    3. The data file shape is used to get number of traces and samples.
    4. The file is closed automatically at end of with...
    5. Then open file again to create the header arrays
        
    inputfile - segy file
    file_headers  - optionally print binary headers to terminal
    DF_ASL - drill floor elevation above Sea Level
    SrcElev - source elevation above sea level
    SRD_ASL - seismic reference datum above sea level
        
    Some useful headers are filled for velocity calculation
        
    For alternatives to my method can use panda dataframes for headers:
        
    https://github.com/equinor/segyio-notebooks/blob/master/notebooks/basic/02_segy_quicklook.ipynb
            
    TRACE_SEQUENCE_LINE= 1
    TRACE_SEQUENCE_FILE= 5
    FieldRecord= 9
    TraceNumber= 13
    EnergySourcePoint= 17
    CDP= 21
    CDP_TRACE= 25
    TraceIdentificationCode= 29
    NSummedTraces= 31
    NStackedTraces= 33
    DataUse= 35
    offset= 37
    ReceiverGroupElevation= 41
    SourceSurfaceElevation= 45
    SourceDepth= 49
    ReceiverDatumElevation= 53
    SourceDatumElevation= 57
    SourceWaterDepth= 61
    GroupWaterDepth= 65
    ElevationScalar= 69
    SourceGroupScalar= 71
    SourceX= 73
    SourceY= 77
    GroupX= 81
    GroupY= 85
    CoordinateUnits= 89
    WeatheringVelocity= 91
    SubWeatheringVelocity= 93
    SourceUpholeTime= 95
    GroupUpholeTime= 97
    SourceStaticCorrection= 99
    GroupStaticCorrection= 101
    TotalStaticApplied= 103
    LagTimeA= 105
    LagTimeB= 107
    DelayRecordingTime= 109
    MuteTimeStart= 111
    MuteTimeEND= 113
    TRACE_SAMPLE_COUNT= 115
    TRACE_SAMPLE_INTERVAL= 117
    GainType= 119
    InstrumentGainConstant= 121
    InstrumentInitialGain= 123
    Correlated= 125
    SweepFrequencyStart= 127
    SweepFrequencyEnd= 129
    SweepLength= 131
    SweepType= 133
    SweepTraceTaperLengthStart= 135
    SweepTraceTaperLengthEnd= 137
    TaperType= 139
    AliasFilterFrequency= 141
    AliasFilterSlope= 143
    NotchFilterFrequency= 145
    NotchFilterSlope= 147
    LowCutFrequency= 149
    HighCutFrequency= 151
    LowCutSlope= 153
    HighCutSlope= 155
    YearDataRecorded= 157
    DayOfYear= 159
    HourOfDay= 161
    MinuteOfHour= 163
    SecondOfMinute= 165
    TimeBaseCode= 167
    TraceWeightingFactor= 169
    GeophoneGroupNumberRoll1= 171
    GeophoneGroupNumberFirstTraceOrigField= 173
    GeophoneGroupNumberLastTraceOrigField= 175
    GapSize= 177
    OverTravel= 179
    CDP_X= 181
    CDP_Y= 185
    INLINE_3D= 189
    CROSSLINE_3D= 193
    ShotPoint= 197
    ShotPointScalar= 201
    TraceValueMeasurementUnit= 203
    TransductionConstantMantissa= 205
    TransductionConstantPower= 209
    TransductionUnit= 211
    TraceIdentifier= 213
    ScalarTraceHeader= 215
    SourceType= 217
    SourceEnergyDirectionMantissa= 219
    SourceEnergyDirectionExponent= 223
    SourceMeasurementMantissa= 225
    SourceMeasurementExponent= 229
    SourceMeasurementUnit= 231
    UnassignedInt1= 233
    UnassignedInt2= 237
    
    * VSProwess puts rcv MD in 37, TVD in 41
    """
    import segyio
    import numpy as np
    
    SrcElev = kwargs['SrcElev'] 
    DF_ASL = kwargs['DF_ASL']
    SRD_ASL = kwargs['SRD_ASL']
    ttbyte = int(kwargs['ttbyte'])
    ttscalar = int(kwargs['ttscalar'])
    file_headers = kwargs['file_headers']
    PTS = kwargs['PTS']

    ###### open the segy file and create numpy data array ################
    
    with segyio.open(inputfile, ignore_geometry=True) as f:
        data = segyio.tools.collect(f.trace[:])
        delta_t = segyio.tools.dt(f)

    nrcv, samples = data.shape

    print (' data.shape :', data.shape)
    
    ###### open the segy file and read trace headers ######################
    
    with segyio.open(inputfile, ignore_geometry=True) as f:
        trnum = np.array(f.attributes(1))
        ffid = np.array(f.attributes(9))
        src = np.array(f.attributes(17))
        nsamp = np.array(f.attributes(115))
        srate = np.array(f.attributes(117))
        zscale = np.array(f.attributes(69))
        mdpth = np.array(f.attributes(37)/abs(zscale))
        tvddpth = np.array(f.attributes(41)/abs(zscale))
        srd = np.array(f.attributes(53)/abs(zscale))
        scalcoord = np.array(f.attributes(71))
        xsrc = np.array(f.attributes(73)/abs(scalcoord))
        ysrc = np.array(f.attributes(77)/abs(scalcoord))
        sdpth = np.array(f.attributes(49)/abs(zscale))
        xrcv = np.array(f.attributes(81)/abs(scalcoord))
        yrcv = np.array(f.attributes(85)/abs(scalcoord))
    
        lagtime_A  = np.array(f.attributes(105))          #trunc. to nearest ms
        lagtime_B  = np.array(f.attributes(107))            # trunc.to nearest ms
        time_byte197  = np.array(f.attributes(197))/100     # divide if used
        ttime  = np.array(f.attributes(ttbyte))/ttscalar     # divide if used
        iline = np.array(f.attributes(189))
        
    ###### open the segy file and print binary and text headers ###########
    
    if (file_headers == 'y') or (file_headers == 'Y'):        
        with segyio.open(inputfile, ignore_geometry=True) as f:
            txt_hed = segyio.tools.wrap(f.text[0]) # [1...] are extended headers
            bin_hed = f.bin           # this is a dictionary with keys and values

            for k, v in bin_hed.items():
                keys = str(k)
                value = int(v)
                print ("\t{:<23} {:<10} ".format(keys, value))
           
            print ('\n',txt_hed)
                
    ############## check if times are in expected header  ######################
    # if nothing bytes 197-200 load times from lag time B
    if np.sum(ttime)==0:
        ttime=lagtime_B

    ############## check if trace number exists  ######################
    # if nothing generate number, starting at one, increment by 1, end at max 
    # trace num 
    if np.sum(trnum)==0:
        trnum=np.arange(1,data.shape[0]+1,1 )

    ############## shift elevations to reference datum ######################
    
    SrcZ = (sdpth *0) + SrcElev # careful, comment out if field is correct
    TVD_Src = tvddpth - (DF_ASL - SrcElev)        
    TVD_SRD = tvddpth - (DF_ASL - SRD_ASL)    
    SrcZ_SRD = SrcZ-SRD_ASL
    
    ############## merge arrays and make a nice table ######################
    
    if (PTS == 'y') or (PTS == 'Y'):        
        from tabulate import tabulate        
        thead_tabl = np.vstack((trnum,mdpth, tvddpth,xrcv, yrcv,  xsrc, ysrc, sdpth, 
                       ttime,TVD_SRD, TVD_Src, SrcZ_SRD, iline, ffid,src)).T        
        print (' table header file shape :', thead_tabl.shape)    

        cheaders2 = ["Trc\nNum", "Rcz\nMD", "Rcz\nTVD","Rcv X", 
                    "Rcv Y","Src X", "Src Y",
                    "Src Z", "Obs\nTime","TVD\nSRD", 
                    "TVD \nSrcZ","SrcZ\nSRD",
                    "ILN", "FFID","Src"]
        numfmt = (".0f",".1f", ".1f", ".1f", ".1f",".1f", ".1f",".1f", ".2f",".1f", ".1f",".1f", ".0f",".0f",".0f")                 
        table = tabulate(thead_tabl, headers = cheaders2,  floatfmt = numfmt)#,tablefmt="pretty")

        print(table)
    
    ############## make a file of trace headers #############################
    
    thead = np.vstack((trnum, mdpth, tvddpth,xrcv, yrcv,  xsrc, ysrc, sdpth, 
                       ttime,TVD_SRD, TVD_Src, SrcZ_SRD, iline, ffid,
                      src))
    ############## sort out sample rates ####################################
    
    numsamp = nsamp[0]
    samprate = srate[0]                    # sample interval in microseconds
    samprate_hz = 1000000/samprate
        
    #tindex = np.zeros(shape = (data.shape[0], int(numsamp*(samprate/1000))))
    # 
    # for k in range(0,data.shape[0]):
    #     tindex[k,:] = np.arange(0, numsamp*(samprate/1000) )  # samprate in ms
        
    ########### print some QC data to screen ################################
    
    print("\u0332".join('\nData Loading with segyio Stats :'))
    print ('\n data shape :', data.shape)    
    print (' data type :', data.dtype)    
    print (' trace header file shape :', thead.T.shape)    
    print (' samples :', data.shape[1],' traces :', data.shape[0], \
           ' fs samprate hz : ', samprate_hz, \
           'samprate milliseconds : ', samprate/1000, \
           '\n numsamp from headers : ', numsamp)    
    print (' first time header value (byte 197 divided by 100) : ',time_byte197[0:1], \
           '\n first lag time A value (byte 105) :', lagtime_A[0:1],
           '\n first lag time B value (byte 107) :', lagtime_B[0:1],
           '\n time saved to header file:', ttime[0:1])
    print (' source depth from header trace 1 :', sdpth[0:1])

    
    return data, numsamp, samprate, samprate_hz, thead.T
    