import sys
import os
import pandas as pd
import sqlite3 
from ms_library_building import ms_rep_library_building


def get_user_input():
    print("Please enter your file ID (e.g., 1, 2..):")    
    input_fileID = sys.stdin.readline()
    file_id = input_fileID.strip()
    file_id = int(file_id)            
    user_para_dict = {
                    'file_id': file_id
                    }
    return user_para_dict


def get_fileID(file_id, conn):
    cursor = conn.cursor()
    sql = "SELECT file_id, file_name, cluster_flag FROM files WHERE file_id = ?"
    val = (int(file_id), )
    cursor.execute(sql, val)
    res = cursor.fetchone()    
    if res is not None:
        file_id, file_name, cluster_flag = res[0], res[1], res[2]
    else:
        print('No such file exists!')
        file_id, file_name, cluster_flag = None, None, None
    return file_id, file_name, cluster_flag


def acquire_Libspectra(file_id, file_name, conn):
    print('Reading library spectral data...')
    query = "SELECT * FROM spectra WHERE file_id= '%s'" % file_id
    precursor_df = pd.read_sql_query(sql=query, con=conn)
    # add a file name column
    precursor_df['file_name'] = file_name 
    spec_id = precursor_df['spectrum_id'].values
    spec_id = ', '.join(map(str, spec_id))
    query2 = "SELECT * FROM masses WHERE spectrum_id IN ({})".format(spec_id)        
    mass_df = pd.read_sql_query(sql=query2, con=conn)    
    # add mass data to spectra table
    mass_v_df = mass_df.groupby('spectrum_id').agg(list).reset_index()
    lib_ms_df = precursor_df.merge(mass_v_df, on='spectrum_id', how='inner')
    return lib_ms_df


def cluster_flag_update(file_id, conn):
    # update files table record: cluster flag, representative flags 
    cursor = conn.cursor()
    sql = "UPDATE files SET cluster_flag = ?, single_representative = ?, average_representative = ? WHERE file_id = ?"
    val = (int(1), 'True', 'True', file_id)
    cursor.execute(sql, val)
    conn.commit()  
    print('files table record updated')
    
    
def cluster_id_update(lib_ms_df, filtered_ms_df, file_id, conn):
    # update cluster id
    lib_ms_df.loc[lib_ms_df['spectrum_id'].isin(filtered_ms_df['spectrum_id']),'cluster_id'] = filtered_ms_df['cluster_id'].values
    lib_ms_df.loc[~lib_ms_df['spectrum_id'].isin(filtered_ms_df['spectrum_id']),'cluster_id'] = -1
    spectrum_id = lib_ms_df['spectrum_id'].values
    cluster_id = lib_ms_df['cluster_id'].values
    cursor = conn.cursor()
    for row in range(len(cluster_id)):
        sql = "UPDATE spectra SET cluster_id = ? WHERE file_id = ? AND spectrum_id = ?"
        val = (int(cluster_id[row]), file_id, int(spectrum_id[row]))
        cursor.execute(sql, val)
    conn.commit()  
    print('cluster_id added!') 
    

def representative_table_gen(target_decoy_spectra_reps, consensus_method, conn):
    # generate representative table in library
    col_drop_names = ['file_name','scan', 'retention_time', 'precursor_feature_id', 'proteoform_id', 'protein_accession', 'protein_description', 'first_residue', 'last_residue',
                       'protein_sequence', 'proteoform', 'e_value', 'spectrum_level_q_value', 'proteoform_level_q_value','cluster_id','flag']    
    target_decoy_spectra_reps = target_decoy_spectra_reps.drop(col_drop_names, axis=1)
    rep_name = consensus_method.lower()
    target_decoy_spectra_reps.to_sql(name = rep_name + '_representatives', con = conn, index=False, if_exists='append') 
    print(rep_name + '_representatives generated for this file!\n')

    
if __name__ == "__main__":
    if len(sys.argv) == 1:
        user_inputs = get_user_input()
        user_file_id = user_inputs['file_id']
        curr_path = os.getcwd()
        directory = "TopLib"
        toplib_filepath = os.path.join(curr_path, directory, "toplib.db")   
        if os.path.isfile(toplib_filepath):
            # Connecting to sqlite
            conn = sqlite3.connect(toplib_filepath)
            # check user's input
            file_id, file_name, cluster_flag = get_fileID(user_file_id, conn)
            if file_id:
                if cluster_flag:
                    print('Spectra representatives of this file already generated in the library!')
                else:
                    lib_ms_df = acquire_Libspectra(file_id, file_name, conn)
                    # generate single representatives
                    data_mode = 'library'
                    print('Building single representative...')
                    ms_df, spectra_single_rep = ms_rep_library_building(lib_ms_df, file_name, 'Single', data_mode)
                    representative_table_gen(spectra_single_rep, 'single', conn)
                    # generate average representatives
                    print('Building average representative...')
                    ms_df, spectra_average_rep = ms_rep_library_building(lib_ms_df, file_name, 'Average', data_mode)
                    representative_table_gen(spectra_average_rep, 'average', conn)
                    # add cluster id to spectra table
                    cluster_id_update(lib_ms_df, ms_df, file_id, conn)
                    # update files table
                    cluster_flag_update(file_id, conn)
            else:
                print('Check your file ID!')
            conn.close()  
        else:
            print("Database .db file does not exist in the TopLib folder.")
    else:
        print("Usage: python script.py")
        sys.exit()



    
   
