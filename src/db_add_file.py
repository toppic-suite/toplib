import sys
import os
import spectra_masses_gen as sm_data
import sqlite3 


def get_user_input():
    print("Please enter your msalign file name (with extension .msalign):")    
    input_msalignfile = sys.stdin.readline()
    msalign_file = input_msalignfile.strip()
    print("Please enter your identification file name (with extension .tsv):")    
    input_tsvfile = sys.stdin.readline()
    tsv_file = input_tsvfile.strip()
    
    print('Please enter your sample ID:')
    input_sampleID = sys.stdin.readline()
    sample_id = input_sampleID.strip()
    sample_id = int(sample_id)
    
    print('Please enter your method ID:')
    input_methodID = sys.stdin.readline()
    method_id = input_methodID.strip()
    method_id = int(method_id)            
    user_para_dict = {
                    'sample_id': sample_id,
                    'method_id': method_id,
                    'file_name': {'msalign_file': msalign_file, 'tsv_file': tsv_file}
                    }
    return user_para_dict


# insert a file record to files table
def file_table_insert(file_name, sample_id, method_id, conn):
    cursor = conn.cursor()
    # check if the record already exists
    sql = "SELECT * FROM files WHERE file_name = ? AND sample_id = ? AND method_id = ?"
    val = (str(file_name), int(sample_id), int(method_id))     
    cursor.execute(sql, val)
    res = cursor.fetchone()   
    if res is None:
        sql = "INSERT INTO files (file_name, sample_id, method_id, cluster_flag, single_representative, average_representative) VALUES (?, ?, ?, ?, ?, ?)"
        val = (str(file_name), int(sample_id), int(method_id), int(0), 'Flase', 'Flase')
        cursor.execute(sql, val)
        conn.commit()
        print("File record add!")
        # get file_id
        sql = "SELECT file_id FROM files WHERE file_name = ? AND sample_id = ? AND method_id = ?"
        val = (str(file_name), int(sample_id), int(method_id))
        cursor.execute(sql, val)
        file_res = cursor.fetchone()
        file_id = file_res[0]
    else:
        print("File record already exists")
        file_id = res[0]
        print(f"current file id is : {file_id}")           
    return file_id

 
def get_sampleID(sample_id, conn):
    cursor = conn.cursor()
    sql = "SELECT sample_id FROM samples WHERE sample_id = ?"
    val = (int(sample_id), )
    cursor.execute(sql, val)
    res = cursor.fetchone()
    if res is not None:
        sample_id = res[0]
    else:
        print('No such sample id exists!')
        sample_id = None
    return sample_id 


def get_methodID(method_id, conn):
    cursor = conn.cursor()
    sql = "SELECT method_id FROM methods WHERE method_id = ?"
    val = (int(method_id), )
    cursor.execute(sql, val)
    res = cursor.fetchone()
    if res is not None:
        method_id = res[0]
    else:
        print('No such method exists!')
        method_id = None
    return method_id 


def get_msfilename(path):
    # remove the path and extension name
    ms_filename_all = os.path.basename(path)
    ms_filename = ms_filename_all.rsplit('.',1)[0]
    return ms_filename


def filename_check(file_names):
    if all(os.path.isfile(f) for f in file_names):
        file_flag = 1
    else:
        file_flag = 0
    return file_flag


if __name__ == "__main__":
    if len(sys.argv) == 1:
        user_inputs = get_user_input()
        user_sample_id = user_inputs['sample_id']
        user_method_id = user_inputs['method_id']
        msalign_filename = user_inputs['file_name']['msalign_file']
        tsv_filename = user_inputs['file_name']['tsv_file']
        curr_path = os.getcwd()
        directory = "TopLib"
        toplib_filepath = os.path.join(curr_path, directory, "toplib.db")   
        if os.path.isfile(toplib_filepath):
            # Connecting to sqlite
            conn = sqlite3.connect(toplib_filepath)
            # check user's input
            sample_id = get_sampleID(user_sample_id, conn)
            method_id = get_methodID(user_method_id, conn)  
            # file name check
            msalign_filename = os.path.join(curr_path, directory, msalign_filename)  
            tsv_filename = os.path.join(curr_path, directory, tsv_filename) 
            filenames = [msalign_filename, tsv_filename]
            file_flag = filename_check(filenames)
            if file_flag == 1:
                # add items to files tables
                file_name = get_msfilename(msalign_filename)
                if sample_id and method_id: 
                    file_id = file_table_insert(file_name, sample_id, method_id, conn)
                    # insert data to tables: spectra and masses
                    sm_data.get_spectra_masses_tables(file_id, msalign_filename, tsv_filename, conn)
                else:
                    print('Check your sample ID and method ID, at least one of them not exist!')
            else:
                print('At least one of input files incorrect and please check file name and path!')
            conn.close()
        else:
            print("Database .db file does not exist in the TopLib folder.")
    else:
        print("Usage: python script.py")
        sys.exit()

