def segyout_top_grid(self,root):
    ''' create a toplevel window (pop-out) to get plot parameters
    # self.master is the root window in a form which can be passed to functions
    may need to get tt byte, ttscalar, DF_ASL, SrcElev, SRD_ASL
    '''
    import tkinter as tk
    from tkinter import ttk
    import procvsp.utils_gui as utilg
    #import timepick.gui_funcs as guifuncs
    from iovsp import segyout_gui


    # if segy in frame (or any other)already exists, delete it and 
    # create a new one in it's place
    
    utilg.grid_clear(self)

    # create a fancy frame for segy load parameter entry in column 0, row 1 of master window
    self.segyout_paramframe = ttk.LabelFrame(self.mframe, text = 'Start Bytes and Scalars ')#, height = 5, width = 800)
    self.segyout_paramframe.grid(row = 0, column = 0, padx = 2, pady=2, sticky='nw')

    # geometry updating variables frame - insert text for user input boxes
    labelttst = tk.Label(master = self.segyout_paramframe,text="TT start byte" )
    labelttsc = tk.Label(master = self.segyout_paramframe,text="Scalar to adjust TT" )
    labelsrcx= tk.Label(master = self.segyout_paramframe,text="Source X" )
    labelsrcy = tk.Label(master = self.segyout_paramframe,text="Source Y" )
    labelsrcz = tk.Label(master = self.segyout_paramframe,text="Source Z" )
    labelrcvx = tk.Label(master = self.segyout_paramframe,text="Receiver X" )
    labelrcvy = tk.Label(master = self.segyout_paramframe,text="Receiver Y" )
    labelrcvz = tk.Label(master = self.segyout_paramframe,text="Receiver Z (MD)" )

    labelttst.grid(row=0, column=0,padx=5, pady=5, sticky='w')
    labelttsc.grid(row=1, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labelsrcx.grid(row=0, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelsrcy.grid(row=1, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelsrcz.grid(row=0, column=4,ipadx=0,padx=5, pady=5, sticky='w')
    labelrcvx.grid(row=1, column=6,ipadx=0,padx=5, pady=5, sticky='w')
    labelrcvy.grid(row=0, column=6,ipadx=0,padx=5, pady=5, sticky='w')
    labelrcvz.grid(row=0, column=8,ipadx=0,padx=5, pady=5, sticky='w')
    
   # geometry updating frame- get user inputs, insert defaults
    self.entryttst_o = tk.Entry(master=self.segyout_paramframe, width=8)
    self.entryttsc_o = tk.Entry(master=self.segyout_paramframe, width=8)
    self.entrysrcx_o = tk.Entry(master=self.segyout_paramframe, width=8)
    self.entrysrcy_o = tk.Entry(master=self.segyout_paramframe, width=8)
    self.entrysrcz_o = tk.Entry(master=self.segyout_paramframe, width=8)
    self.entryrcvx_o = tk.Entry(master=self.segyout_paramframe, width=8)
    self.entryrcvy_o = tk.Entry(master=self.segyout_paramframe, width=8)
    self.entryrcvz_o = tk.Entry(master=self.segyout_paramframe, width=8)
   
    self.entryttst_o.grid(row=0, column=1, padx=5, pady=5, sticky='w')
    self.entryttsc_o.grid(row=1, column=1, padx=5, pady=5, sticky='w')
    self.entrysrcx_o.grid(row=0, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entrysrcy_o.grid(row=1, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entrysrcz_o.grid(row=0, column=5, ipadx=0,padx=5, pady=5, sticky='w')
    self.entryrcvx_o.grid(row=0, column=7, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryrcvy_o.grid(row=1, column=7, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryrcvz_o.grid(row=0, column=9, ipadx=0, padx=5, pady=5, sticky='w')

   # Set some default values which are updated by SEG-Y header values
    self.entryttst_o.insert(0, '107')
    self.entryttsc_o.insert(0, '100')
    self.entrysrcx_o.insert(0, '73')
    self.entrysrcy_o.insert(0, '77')
    self.entrysrcz_o.insert(0, '49')
    self.entryrcvx_o.insert(0, '81')
    self.entryrcvy_o.insert(0, '85')
    self.entryrcvz_o.insert(0, '37')

    # create checkbuttons for writing tables
    self.hdrschk_o = tk.StringVar()
    self.tblchk_o = tk.StringVar()
    
    self.hdrschk_o.set('y')
    self.tblchk_o.set('n')

    buttontxtedt = ttk.Button(self.segyout_paramframe, text = "Edit Text Header", command = lambda: segyout_gui.texth_edit(self,root))#self.pick_segy)        
    buttontxtedt.grid(row=0, column=11, columnspan=1, sticky='w',padx=10, pady=5 )   
    buttonsegout = ttk.Button(self.segyout_paramframe, text = "Write SEG-Y", command = lambda: segyout_gui.write_segy(self,root))#self.pick_segy)        
    buttonsegout.grid(row=1, column=11, columnspan=1, sticky='w',padx=10, pady=5 )   
 
    self.segyout_menu_exists = self.segyout_paramframe.winfo_exists()

def get_segyout_params(self):
    ''' Get SEG-Y output variables from gui and write to dictionary which
    can be passed around.
    '''    
    segyout_params=dict()

    segyout_params['timbit']=int(self.entryttst_o.get())
    segyout_params['timscal']=int(self.entryttsc_o.get())
    segyout_params['sxbit'] =int(self.entrysrcx_o.get())
    segyout_params['sybit'] =int(self.entrysrcy_o.get())
    segyout_params['szbit'] =int(self.entrysrcz_o.get())
    segyout_params['rxbit'] =int(self.entryrcvx_o.get())
    segyout_params['rybit'] =int(self.entryrcvy_o.get())
    segyout_params['rzbit'] =int(self.entryrcvz_o.get())

    # variable to indicate tkinter useage
    segyout_params['tkinter_use']='y'
    
    return segyout_params

def segyout(self):
    ''' set up and plot parameter entry grid.
    Action button invokes pick_segy'''

    from iovsp import segyin_gui

    segyin_gui.segy_top_grid(self)

def texth_edit(self,root):
        ''' Text box for the text header, whixh can be edited 
        prior to writing segy file.
        Has a save button for writing the text header to an external file.
        '''
        import tkinter as tk
        from tkinter import ttk
        from iovsp import segyout_gui
        from iovsp import text_io

        #print ( "txt_hdr_out :",txt_hdr_out)
        # create a toplevel window (pop-out) to print headers
        # then add a text widget to the window
        thedwindow = tk.Toplevel(root)
        thedwindow.wm_title("Text Header")# double quotations!
 
        # this section is necessary to make scrollbars work
        # scrollbars need to be attached to a re-sizable Frame       
        self.thedframe = ttk.Frame(thedwindow)
        self.thedframe.pack(fill="both", expand=True)            
        self.thedframe.rowconfigure(0, weight=1)
        self.thedframe.columnconfigure(0, weight=1)# this column will expand
        
        txt_hed = tk.Text(self.thedframe)

        txt_hed.grid(row = 0, column = 0, sticky='nsew')# or tk.HORIZONTAL
        txt_hed.columnconfigure(0,  weight=1)
        txt_hed.rowconfigure(0,  weight=1)
        txt_hed.config(font=("consolas", 10), wrap="none", borderwidth=3, relief="sunken")

        ##### Text Header Scroll bars ######
        lfxscrollbar = ttk.Scrollbar(self.thedframe, command=txt_hed.xview, orient = tk.HORIZONTAL)
        lfyscrollbar = ttk.Scrollbar(self.thedframe, command=txt_hed.yview, orient = tk.VERTICAL)
        lfxscrollbar.grid(row=1, column=0,sticky='ew')
        lfyscrollbar.grid(row=0, column=1,sticky='ns')

        txt_hed.configure(xscrollcommand = lfxscrollbar.set)
        txt_hed.configure(yscrollcommand = lfyscrollbar.set)   
        
        # get the external text header file
        infile=('iovsp\\text_header.py')
        text_header_read = text_io.import_ascii(infile)
        # insert text header file to tkinter text box
        # Editing can take place
        self.title_sgout='SEG-Y Text Header'
        txt_hed.insert(tk.END,text_header_read)
        save_button = tk.Button(self.thedframe, text="Save Changes", command=lambda: segyout_gui.save_dictionary(self,txt_hed))
        save_button.grid(row = 1, column = 0, ) 
    
def save_dictionary(self,text_win):
    ''' Save text header from tkinter text box as an ascii file.
        No choice of output names, write segy assumes hardwired name.
    ''' 
    # from the text window, extract the text header after edit
    filetext = text_win.get("1.0", "end-1c")

    save_text='iovsp\\text_header_edit.py'
    if save_text:
        with open(save_text, "w") as f:
            f.write(filetext)
    self.new_txt_hdr=filetext

def write_segy(self,root):
    ''' Started from segyout parameter grid action button. 
    Uses a windows file selection widget to get the output path.
    Uses the updated text header created by save_dictionary.
    '''    
    from iovsp import segyout
    from tkinter import filedialog as fd
    from pathlib import Path
    
    # test if an updated text header exists.
    # if not, use the default text header
    try:
        from iovsp.text_header_edit import text_header
        print('\n updated file found')
    except ModuleNotFoundError:
        from iovsp.text_header import text_header
        print ('\n using default file')
    #
 
    # use a windows dialog box to set the output path and file name
    file_path=fd.asksaveasfilename(initialdir = "c:\\Users\\acampbell45\\Documents",
                                    defaultextension=".sgy",
                                    title = "Save SEG-Y file",
                                    filetypes = (("seg-y files","*.sgy"),("all files","*.*")))
    
    # get byte positions, scalars, from gui
    segyout_params=get_segyout_params(self)
    segyout.write_segyio(self.data,text_header, self.headers,self.fs,file_path,**segyout_params)
    

