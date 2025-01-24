from datetime import datetime
import sys
import pandas as pd
import sqlite3 

def get_rep_data(lib_filename): 
    conn = sqlite3.connect(lib_filename) 
    query = "SELECT * FROM target_decoy_spectra_representatives"
    target_decoy_spectra_rep = pd.read_sql_query(query, con=conn)
    print(len(target_decoy_spectra_rep))
    conn.close()
    return target_decoy_spectra_rep


def reviewed_status(target_decoy_spectra_rep, lib_filename):
    reviewed_time = datetime.now().strftime('%Y/%m/%d, %H:%M:%S')
    target_decoy_spectra_rep.loc[target_decoy_spectra_rep['num_unexpected_modifications']==0,'reviewed'] = 'True'
    target_decoy_spectra_rep.loc[target_decoy_spectra_rep['num_unexpected_modifications']==0,'reviewed_time'] = reviewed_time
    proteoforms_val = target_decoy_spectra_rep.loc[target_decoy_spectra_rep['num_unexpected_modifications']==0,'proteoform'].values
    target_decoy_spectra_rep.loc[target_decoy_spectra_rep['num_unexpected_modifications']==0,'reviewed_proteoform'] = proteoforms_val
    conn = sqlite3.connect(lib_filename)
    target_decoy_spectra_rep.to_sql(name = 'target_decoy_spectra_representatives', con = conn, index=False, if_exists='replace') 
    conn.close()

     
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_db_filename>")
        sys.exit()
    else:
        lib_filename = sys.argv[1] # ex. sw480_rep1_combined_ms2_single_representatives.db
        spectra_ms_rep = get_rep_data(lib_filename)
        reviewed_status(spectra_ms_rep, lib_filename)
    
