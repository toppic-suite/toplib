# convert db file to msalign file and tsv file
import sys
import os
import pandas as pd
# import numpy as np
import sqlite3 
# from tqdm import tqdm
from db_add_file import get_msfilename
from db_query2 import rep2msalign, rep2tsv
# from lib_representative_extraction import lib_extraction_output

    
def get_spectra_representatives(toplib_filepath):
    # read toplib data
    conn = sqlite3.connect(toplib_filepath)
    query = "SELECT * FROM target_decoy_spectra_representatives"
    target_decoy_rep = pd.read_sql_query(query, con=conn)  
    conn.close()
    mass_inte_ch_rep = target_decoy_rep[target_decoy_rep['flag']==1]
    mass_inte_ch_rep = mass_inte_ch_rep.sort_values(by='spectrum_id')
    # print(len(mass_inte_ch_rep))
    # print(mass_inte_ch_rep)
    return mass_inte_ch_rep
       
        
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_db_file>")
        sys.exit()
    else:
        curr_path = os.getcwd()
        directory = "TopLib"
        db_filename = sys.argv[1]
        db_filename = os.path.join(curr_path, directory, db_filename)     
        file_flag = os.path.isfile(db_filename)
        if file_flag==1:
            rep_df = get_spectra_representatives(db_filename)
            # write msalign file
            filename = get_msfilename(db_filename)
            msalign_wfile = filename + ".msalign"
            rep2msalign(rep_df, msalign_wfile)            
            # write tsv file
            tsv_wfile = filename + ".tsv"
            rep2tsv(rep_df, tsv_wfile)
            # write all msalign file
            # msalign_wfile2 = filename + "_all.msalign"
            # lib_extraction_output(rep_df, msalign_wfile2)
        else:
            print("Database .db file does not exist in the TopLib folder.")
        


