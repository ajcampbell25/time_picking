r''' Read, plot and time pick VSP data

Operates on one file or multiple files

time_pick_gui_v0 - adapt gui to time picking
time_pick_gui_v1 - adapt gui to time picking
    - import all functions instead of just algorthmic ones
    - refactor functions to accomodate above 

V0 has following options:
1. Load a SEG-Y file
2. Look at 1D spectra
3. Run a Butterworth band-pass filter
4. Dump SEG-Y binary, text and trace headers
5. Run STA-LTA time picking
6. Manual time picking with mouse clicking

'''
from tkinter import messagebox
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
import tkinter.font as tkFont
import procvsp.utils_gui as utilg

class GUI(tk.Frame):
    def __init__(self,master):
        super().__init__()
        
        ''' Insert a frame into the root window
        # The pack and configure statements seem necessary to get the plot_frame 
        # to expand properly and scrollbars to work properly
        # The dropdown menu occupies row 0
        # the plotting area occupies row 1
        '''        
        self.master = master

        default_font = tkFont.nametofont("TkDefaultFont")
        default_font.configure(size=8)
        root.option_add("*Font", default_font)
        
        self.mframe = ttk.Frame(master)
        self.mframe.pack(side="top", fill="both", expand=True)     
        self.mframe.rowconfigure(1, weight=1)# make row 1 expand to fit all available room
        self.mframe.columnconfigure(0, weight=1)# this column will expand
        
        # prevent parameter grids from overwriting each other        
        utilg.grid_init(self)
        # add a plotting area
        utilg.plot_area(self)
        # add dropdown menu
        utilg.drop_down_menu(self,root)
     
    def dummy(self):
        ''' useful for testing, does nothing
        '''
        pass

    def output_path(self):        
        import io

        self.outdir = fd.askdirectory(initialdir = "c:\\Users\\acampbell45\\Documents",
                    title = "Select path for outputs")
        print (' self.outdir path for outputs :', self.outdir)
        
    def exit_function(self):
        import procvsp.utils_gui as utilg
        utilg.fig_exist()
        if messagebox.askyesno("Close the window", "Do you want to close the window?", icon='warning'):
            root.destroy()
        else:
            pass
    
    def about(self):
        messagebox.showinfo("About SEG-Y Load and Time Picker", "Written by Allan Campbell")

#class Text(tk.Frame):
#    def __init__(self, master):
#        tk.Frame.__init__(self, master)

if __name__=='__main__':
    root = tk.Tk()
    root.wm_geometry("1400x900")
    root.wm_title('SEG_Y Reader and Time Picker')
    ttk.Style().theme_use('xpnative')
    ttk.Style().configure('TPanedwindow', background='black')

    app = GUI(root)
    app.mainloop()
