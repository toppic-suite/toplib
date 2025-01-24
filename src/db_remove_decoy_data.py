import sys
import pandas as pd
import sqlite3 


def get_rep_target_data(lib_filename):
    # get the representatives data    
    conn = sqlite3.connect(lib_filename) 
    query = "SELECT * FROM target_decoy_spectra_representatives"
    target_decoy_spectra_rep = pd.read_sql_query(query, con=conn)
    target_spectra_rep = target_decoy_spectra_rep[target_decoy_spectra_rep['flag']==1]
    # remove normalized intensity
    columns = [f"{prefix}_{i}" for prefix in ['norm_intensity'] for i in range(50)]
    target_spectra_rep = target_spectra_rep.drop(columns=columns, axis=1) 
    print(len(target_spectra_rep))
    conn.close()
    return target_spectra_rep


def store_target_rep(target_spectra_rep, lib_new_filename):
    conn = sqlite3.connect(lib_new_filename)
    target_spectra_rep.to_sql(name = 'target_decoy_spectra_representatives', con = conn, index=False, if_exists='replace') 
    conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_db_filename> <input_db_new_filename>")
        sys.exit()
    else:
        lib_filename = sys.argv[1] # ex. sw480_rep1_combined_ms2_single_representatives_new4.db
        lib_new_filename = sys.argv[2] 
        target_spectra_rep = get_rep_target_data(lib_filename)
        store_target_rep(target_spectra_rep, lib_new_filename)

        