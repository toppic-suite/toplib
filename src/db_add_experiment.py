import sys
import os
import sqlite3 


def get_user_input():
    print('Please enter your experiment name:')
    input_exp_name = sys.stdin.readline()
    exp_name = input_exp_name.strip()
    print('Please enter your experiment description:')
    input_exp_des = sys.stdin.readline()
    exp_description = input_exp_des.strip()   
    print('Please enter your sample ID:')
    input_sampleID = sys.stdin.readline()
    sample_id = input_sampleID.strip()
    sample_id = int(sample_id)
    print('Please enter your method ID:')
    input_methodID = sys.stdin.readline()
    method_id = input_methodID.strip()
    method_id = int(method_id)                    
    user_para_dict = {
                    'experiment_name': exp_name, 'experiment_description': exp_description,
                    'sample_id': sample_id,
                    'method_id': method_id
                    }
    return user_para_dict


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


# insert a sample record to sample table
def experiment_table_insert(experiment_name, experiment_description, sample_id, method_id, conn):
    cursor = conn.cursor()
    # check if the record already exists
    sql = "SELECT * FROM experiments WHERE experiment_name = ? AND experiment_description = ? AND sample_id = ? AND method_id = ?"
    val = (str(experiment_name), str(experiment_description), sample_id, method_id)
    cursor.execute(sql, val)
    res = cursor.fetchone()
    if res is None:
        sql = "INSERT INTO experiments (experiment_name, experiment_description, sample_id, method_id) VALUES (?, ?, ?, ?)"
        val = (str(experiment_name), str(experiment_description), sample_id, method_id)
        cursor.execute(sql, val)
        conn.commit()
        print("Experiment record added!")
    else:
        print("Experiment record already exists")
        experiment_id = res[0]
        print(f"current experiment id is : {experiment_id}")



if __name__ == "__main__":
    if len(sys.argv) == 1:
        user_inputs = get_user_input()
        experiment_name = user_inputs['experiment_name']
        experiment_description = user_inputs['experiment_description']
        user_sample_id = user_inputs['sample_id']
        user_method_id = user_inputs['method_id']
        curr_path = os.getcwd()
        directory = "TopLib"
        toplib_filepath = os.path.join(curr_path, directory, "toplib.db")   
        if os.path.isfile(toplib_filepath):
            # Connecting to sqlite
            conn = sqlite3.connect(toplib_filepath) 
            sample_id = get_sampleID(user_sample_id, conn)
            method_id = get_methodID(user_method_id, conn)
            if sample_id and method_id:
                experiment_table_insert(experiment_name, experiment_description, sample_id, method_id, conn)
            else:
                print('Check your sample ID and method ID, at least one of them not exist!')
            conn.close()
        else:
            print("Database .db file does not exist in the TopLib folder.")
    else:
        print("Usage: python script.py")
        sys.exit()
        