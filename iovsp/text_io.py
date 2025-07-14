def import_ascii(file_path):
    ''' useful for importing ascii such as text header
    
    Useage:
    # test reading from file
    infile=('iovsp\\text_header_edit.py')

    text_header_testread = segyout.import_ascii(infile)
    print ('text_header_testread :',text_header_testread)
    '''
    with open(file_path, 'r') as file:
        content = file.read()

    return content

def export_ascii(file_path_out, outtext):
    ''' write out header listings etc.
    '''
    with open (file_path_out, 'w+') as file:
        file.write(outtext)
        file.close()

def header_table(theaders, numcols):
    ''' print a table of headers including sta lta results
    Standard file from SEG-Y has 14 columns
    STA LTA has 16 columns, 
    15 = sta lta results
    16 = tuned sta lta results
    '''
    from tabulate import tabulate
    import numpy as np     
    
    tabl_head=["Trc\nNum", "Rcz\nMD", "Rcz\nTVD","Rcv X", 
                    "Rcv Y","Src X", "Src Y",
                    "Src Z", "Obs\nTime","TVD\nSRD", 
                    "TVD \nSrcZ","SrcZ\nSRD",
                    "ILN", "FFID","Src"]
    numfmt = (".0f",".1f", ".1f", ".1f", ".1f",".1f", ".1f",
                  ".1f", ".2f",".1f", ".1f",".1f", ".0f",".0f",".0f")
    if numcols!=15:
        # STA LTA adds 2 columns to the header file
        tabl_head=np.hstack((tabl_head,"STA-LTA","Tuned TT"))
        numfmt=np.hstack((numfmt,".2f",".2f",))        
    table = tabulate(theaders,tabl_head,floatfmt = numfmt )#, tablefmt="fancy_grid")
    print (table)
    return table