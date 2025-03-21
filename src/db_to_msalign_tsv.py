# convert db file to msalign file and tsv file
import sys
import os
import pandas as pd
import sqlite3 
from db_add_file import get_msfilename
from db_query import rep2msalign, rep2tsv

    
def get_spectra_representatives(toplib_filepath):
    # read toplib data
    conn = sqlite3.connect(toplib_filepath)
    query = "SELECT * FROM target_decoy_spectra_representatives"
    target_decoy_rep = pd.read_sql_query(query, con=conn)  
    conn.close()
    mass_inte_ch_rep = target_decoy_rep[target_decoy_rep['flag']==1]
    mass_inte_ch_rep = mass_inte_ch_rep.sort_values(by=['file_name','spectrum_id'])
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
            msalign_wfile = filename + "_representatives.msalign"
            rep2msalign(rep_df, msalign_wfile)            
            # write tsv file
            tsv_wfile = filename + "_identifications.tsv"
            rep2tsv(rep_df, tsv_wfile)
        else:
            print("Database .db file does not exist in the TopLib folder.")
        


