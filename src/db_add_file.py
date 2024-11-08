import sys
import os
import argparse
import pandas as pd
import spectra_masses_gen as sm_data
import sqlite3 
from ms_library_building import ms_rep_library_building


def user_command():
    parser = argparse.ArgumentParser(description="Top-down library file parameters specification")
    parser.add_argument("msalign_filename", type=str,
                        help="An msalign file name (*.msalign)")
    parser.add_argument("tsv_filename", type=str,
                        help="A tsv file name (*.tsv)")
    parser.add_argument("experiment_id", type=int,  
                        help="Experiment ID (ex.,1)")
    args = parser.parse_args()
    msalign_filename = args.msalign_filename
    tsv_filename = args.tsv_filename
    experiment_id = int(args.experiment_id)
    return msalign_filename, tsv_filename, experiment_id 


# insert a file record to files table
def file_table_insert(file_name, experiment_id, conn):
    cursor = conn.cursor()
    # check if the record already exists
    sql = "SELECT * FROM files WHERE file_name = ? AND experiment_id = ?"
    val = (str(file_name), int(experiment_id))     
    cursor.execute(sql, val)
    res = cursor.fetchone()   
    if res is None:
        sql = "INSERT INTO files (file_name, experiment_id, cluster_flag, single_representative, average_representative) VALUES (?, ?, ?, ?, ?)"
        val = (str(file_name), int(experiment_id), int(0), 'Flase', 'Flase')
        cursor.execute(sql, val)
        conn.commit()
        print("File record add!")
        # get file_id
        sql = "SELECT file_id FROM files WHERE file_name = ? AND experiment_id = ?"
        val = (str(file_name), int(experiment_id))
        cursor.execute(sql, val)
        file_res = cursor.fetchone()
        file_id = file_res[0]
    else:
        print("File record already exists")
        file_id = res[0]
        print(f"current file id is : {file_id}")           
    return file_id


def get_experimentID(experiment_id, conn):
    cursor = conn.cursor()
    sql = "SELECT experiment_id FROM experiments WHERE experiment_id = ?"
    val = (int(experiment_id), )
    cursor.execute(sql, val)
    res = cursor.fetchone()
    if res is not None:
        experiment_id = res[0]
    else:
        print('No such method exists!')
        experiment_id = None
    return experiment_id 


def get_msfilename(path):
    # remove the path and extension name
    ms_filename_all = os.path.basename(path)
    ms_filename = ms_filename_all.rsplit('.',1)[0]
    return ms_filename


def get_fileName(file_id, conn):
    cursor = conn.cursor()
    sql = "SELECT file_name FROM files WHERE file_id = ?"
    val = (int(file_id), )
    cursor.execute(sql, val)
    res = cursor.fetchone()
    file_name = res[0]
    return file_name 


def acquire_Libspectra(file_id, conn):
    print('Reading library spectral data...')
    query = "SELECT * FROM spectra WHERE file_id= '%s'" % file_id
    precursor_df = pd.read_sql_query(sql=query, con=conn)
    # add a file name column
    file_name = get_fileName(file_id, conn)
    precursor_df['file_name'] = file_name 
    spec_id = precursor_df['spectrum_id'].values
    spec_id = ', '.join(map(str, spec_id))
    query2 = "SELECT * FROM masses WHERE spectrum_id IN ({})".format(spec_id)        
    mass_df = pd.read_sql_query(sql=query2, con=conn)    
    # combine mass data and spectra data
    mass_v_df = mass_df.groupby('spectrum_id').agg(list).reset_index()
    lib_ms_df = precursor_df.merge(mass_v_df, on='spectrum_id', how='inner')
    return lib_ms_df


def get_clusterFlag(file_id, conn):
    cursor = conn.cursor()
    sql = "SELECT cluster_flag FROM files WHERE file_id = ?"
    val = (int(file_id), )
    cursor.execute(sql, val)
    res = cursor.fetchone()
    cluster_flag = res[0]
    return cluster_flag


def cluster_flag_update(file_id, conn):
    # update files table record: cluster flag, representative flags 
    cursor = conn.cursor()
    sql = "UPDATE files SET cluster_flag = ?, single_representative = ?, average_representative = ? WHERE file_id = ?"
    val = (int(1), 'True', 'True', file_id)
    cursor.execute(sql, val)
    conn.commit()  
    print('files table record updated')
    
        
def cluster_id_update(lib_ms_df, clustered_ms_df, file_id, conn):
    # update cluster id
    lib_ms_df.loc[lib_ms_df['spectrum_id'].isin(clustered_ms_df['spectrum_id']),'cluster_id'] = clustered_ms_df['cluster_id'].values
    lib_ms_df.loc[~lib_ms_df['spectrum_id'].isin(clustered_ms_df['spectrum_id']),'cluster_id'] = -1
    spectrum_id = lib_ms_df['spectrum_id'].values
    cluster_id = lib_ms_df['cluster_id'].values
    cursor = conn.cursor()
    for row in range(len(cluster_id)):
        sql = "UPDATE spectra SET cluster_id = ? WHERE file_id = ? AND spectrum_id = ?"
        val = (int(cluster_id[row]), file_id, int(spectrum_id[row]))
        cursor.execute(sql, val)
    conn.commit()  
    print('cluster_id added!') 


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_msalign_file> <input_tsv_file> <input_experiment_id>")
        sys.exit()
    else:
        msalign_filename, tsv_filename, user_experiment_id = user_command() 
        curr_path = os.getcwd()
        directory = "TopLib"
        toplib_filepath = os.path.join(curr_path, directory, "toplib.db")   
        if os.path.isfile(toplib_filepath):
            # Connecting to sqlite
            conn = sqlite3.connect(toplib_filepath)
            # check user's input
            experiment_id = get_experimentID(user_experiment_id, conn)
            # file name check
            filenames = [msalign_filename, tsv_filename]
            if all(os.path.isfile(f) for f in filenames):
                # add items to files tables
                file_name = get_msfilename(msalign_filename)
                if experiment_id: 
                    file_id = file_table_insert(file_name, experiment_id, conn)
                    # insert data to tables: spectra and masses
                    sm_data.get_spectra_masses_tables(file_id, msalign_filename, tsv_filename, conn)
                    cluster_flag = get_clusterFlag(file_id, conn)
                    if cluster_flag:
                        print('Spectra representatives of this file already generated in the library!')
                    else:
                        lib_ms_df = acquire_Libspectra(file_id, conn)
                        print('Building single representative...')
                        clustered_ms_df = ms_rep_library_building(lib_ms_df.copy(), 'single', 'library')
                        print('Building average representative...')
                        clustered_ms_df = ms_rep_library_building(lib_ms_df.copy(), 'average', 'library')
                        # add cluster id to spectra table
                        cluster_id_update(lib_ms_df, clustered_ms_df, file_id, conn)
                        # update files table
                        cluster_flag_update(file_id, conn)                                                                
                else:
                    print('Check your experiment ID, no such experiment id in the library!')
            else:
                print('At least one of input files incorrect and please check file name and path!')
            conn.close()
        else:
            print("Database .db file does not exist in the TopLib folder.")


