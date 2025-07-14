def drop_down_menu(self,root):
    ''' Create a pull-down menu in the main window.
        The menu will be placed in row 0 of mframe
    '''
    import tkinter as tk

    from iovsp import segyin_gui
    from iovsp import segyout_gui

    from procvsp import spec_gui
    from plotvsp import seisplots_gui
    from timepick import pickers_gui
    from timepick import pick_plot_gui
    from timepick import headers_gui


    self.menubar = tk.Menu(self.mframe)
    root.config(menu= self.menubar)
    # Need to use lambda to reate an anonymous function which calls commands,
    # or the commands will run immediately.
    # This seems to occur when passing root to the function drop_down_menu
    self.filemenu = tk.Menu(self.menubar,tearoff=0)
    self.menubar.add_cascade(label = "File", menu = self.filemenu)
    self.filemenu.add_command(label = "Open", 
                              command = lambda: segyin_gui.segyread(self))
    self.filemenu.add_command(label = "Spectra", 
                              command = lambda: spec_gui.spectrum(self,root))
    self.filemenu.add_command(label = "Spectra F-Z", 
                              command = lambda: spec_gui.spectrumFZ(self,root)) 
    self.filemenu.add_command(label = "Filter", 
                              command = lambda: spec_gui.bpf_filt(self,root)) 
    self.filemenu.add_command(label = "Plot", 
                              command = lambda: seisplots_gui.plot_stack(self))
    self.filemenu.add_separator()               
    self.filemenu.add_command(label = "Exit", 
                              command = self.exit_function)

    self.headmenu = tk.Menu(self.menubar,tearoff=0)
    self.menubar.add_cascade(label = "Headers", menu = self.headmenu)
    self.headmenu.add_command(label = "Text", 
                              command = lambda: headers_gui.file_heads(self,root))
    self.headmenu.add_command(label = "Trace", 
                              command = lambda: headers_gui.trc_heads(self,root))
    self.headmenu.add_separator()
    self.headmenu = tk.Menu(self.menubar,tearoff=0)
    self.menubar.add_cascade(label = "Time Pick", menu = self.headmenu)
    self.headmenu.add_command(label = "STA LTA", 
                              command = lambda: pickers_gui.tpicks(self,root))
    self.headmenu.add_command(label = "Manual Pick", 
                              command = lambda: pick_plot_gui.man_pick(self,root))
    self.headmenu.add_command(label = "Save SEG-Y", 
                              command = lambda: segyout_gui.segyout_top_grid(self,root))
    self.headmenu.add_command(label = "Save ASCII", 
                              command = lambda: headers_gui.trc_heads(self,root))
    self.aboutmenu = tk.Menu(self.menubar,tearoff=0)
    self.aboutmenu.add_command(label="About", command="")

def plot_area(self):
    r''' create a frame (window) to place plots inside.
    We put the label frames for menus into row 0 of mframe
    We put plot_frame into  row 1 of mframe 
    '''
    import tkinter as tk
    
    self.plot_frame = tk.Frame(self.mframe, borderwidth=2,background="green")

    self.plot_frame.grid(row = 1, column = 0, padx=0, pady=0, sticky='nsew',)#,rowspan = 6, columnspan = 9
    self.plot_frame.grid_columnconfigure(0, weight=1)
    self.plot_frame.grid_rowconfigure(0, weight=1)


def grid_init(self):
    ''' set the default values to 0
    parameter grids do not exist yet
    '''
    self.plot_menu_exists = 0
    self.bpf_menu_exists = 0
    self.segyin_menu_exists = 0
    self.segyout_menu_exists = 0
    self.spec_menu_exists = 0
    self.tpick_menu_exists = 0
    self.mpick_menu_exists = 0
    self.fzspec_menu_exists = 0

def grid_clear(self):
    ''' remove existing parameter grids prior to inserting
    new parameter grids
    '''
    if self.plot_menu_exists ==1:
        self.plot_paramframe.destroy()
        self.plot_menu_exists = 0
    if self.segyin_menu_exists ==1:
        self.segyin_paramframe.destroy()
        self.segyin_menu_exists = 0
    if self.segyout_menu_exists ==1:
        self.segyout_paramframe.destroy()
        self.segyout_menu_exists = 0
    if self.spec_menu_exists == 1:
        self.spec_paramframe.destroy()
        self.spec_menu_exists = 0
    if self.bpf_menu_exists == 1:
        self.bpf_paramframe.destroy()
        self.bpf_menu_exists = 0
    if self.tpick_menu_exists == 1:
        self.tpick_paramframe.destroy()
        self.tpick_menu_exists = 0
    if self.mpick_menu_exists ==1:
        self.mpick_paramframe.destroy()
        self.mpick_menu_exists==0
    if self.fzspec_menu_exists ==1:
        self.fzspec_paramframe.destroy()
        self.fzspec_menu_exists==0

def create_toplevel(root,frame_name):
    import tkinter as tk
    from tkinter import ttk

    ''' Create a toplevel window (pop-out) to plot various things.
    Make the window expandable by adding a pack frame
    '''
    topwindow = tk.Toplevel(root)
    
    topwindow.geometry('1300x800') # depends on fig. size 
    topwindow.wm_title(frame_name)# double quotations!

    # this section is necessary to make scrollbars work
    # scrollbars need to be attached to a re-sizable Frame       
    topframe = ttk.Frame(topwindow)
    topframe.pack(fill="both", expand=True)            
    topframe.rowconfigure(0, weight=1)
    topframe.columnconfigure(0, weight=1)# this column will expand

    return topframe

def add_scrollbars(frame_in):
    ''' create a canvas, than add scroll bars to the canvas.
    This canvas is plotted on the expandable frame in a top-level window.
    '''
    import tkinter as tk
    # create a tkinter canvas and scroll bars
    canvas_scroll = tk.Canvas(frame_in,borderwidth=1,relief=tk.RIDGE)
    canvas_scroll.grid(row=0, column=0, sticky=tk.NSEW)

    #Trace Plot Scroll bars
    xScrollbarsta = tk.Scrollbar(frame_in, orient=tk.HORIZONTAL, width=25)
    yScrollbarsta = tk.Scrollbar(frame_in, width=25)                
        
    xScrollbarsta.grid(row=1, column=0, sticky=tk.EW)
    yScrollbarsta.grid(row=0, column=1, sticky=tk.NS)

    canvas_scroll.config(xscrollcommand=xScrollbarsta.set)
    canvas_scroll.config(yscrollcommand=yScrollbarsta.set)
#    speed = 2

    xScrollbarsta.config(command=canvas_scroll.xview)
    yScrollbarsta.config(command=canvas_scroll.yview)

    return canvas_scroll

def fig_canvas(fig,canvas_in):
    ''' put the figure on to the canvas
    '''
    import tkinter as tk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    figAgg = FigureCanvasTkAgg(fig, canvas_in)
    mplCanvas = figAgg.get_tk_widget()
    canvas_in.create_window(0, 0, window=mplCanvas, anchor=tk.NW)
    canvas_in.config(scrollregion=canvas_in.bbox(tk.constants.ALL))

    return canvas_in

def fig_exist():
    '''Check if there are any open figures
    If yes, close them.
    Prevents hangs on exiting the gui.
    '''
    import matplotlib.pyplot as plt
    # Check if any figure exists
    if not plt.get_fignums():
        print("No figure exists.")
    else:
        print("At least one figure exists.")
        plt.close('all')

