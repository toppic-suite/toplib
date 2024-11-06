import sys
import os
import sqlite3 


def get_user_input():
    print('Please enter your specie ID: 1 = human; 2 = mouse.')
    input_species = sys.stdin.readline()
    species_code = input_species.strip()
    if species_code == '1':
        species_name = 'human'
    elif species_code == 2:
        species_name = 'mouse'
    else:
        print("No such species available and please choose one of species: 1 = human or 2 = mouse.")
        return 
    
    print('Please enter your sample name:')
    input_sample = sys.stdin.readline()
    sample_name = input_sample.strip()
    print('Please enter your sample description:')
    input_sample_des = sys.stdin.readline()
    sample_description = input_sample_des.strip()   
    print('Please enter your project ID:')
    input_project_id = sys.stdin.readline()
    project_id = input_project_id.strip()
    project_id = int(project_id)
    user_para_dict = {
                    'sample_name': sample_name, 'sample_description': sample_description,
                    'species_name': species_name,
                    'project_id': project_id
                    }
    return user_para_dict


def get_projectID(project_id, conn):
    # get the project_id
    cursor = conn.cursor()
    sql = "SELECT project_id FROM projects WHERE project_id = ?"
    val = (int(project_id), )
    cursor.execute(sql, val)
    res = cursor.fetchone() 
    if res is not None:
        project_id = res[0]
    else:
        print('no such project id exists!')
        project_id = None
    return project_id 


def get_speciesID(species_name, conn):
    # get species id
    cursor = conn.cursor()
    sql = "SELECT species_id FROM species WHERE species_name = ?"
    val = (str(species_name), )
    cursor.execute(sql, val)
    res = cursor.fetchone() 
    species_id = res[0]
    return species_id


# insert a sample record to sample table
def sample_table_insert(sample_name, sample_description, species_id, project_id, conn):
    cursor = conn.cursor()
    # check if the record already exists
    sql = "SELECT * FROM samples WHERE sample_name = ? AND sample_description = ? AND species_id = ? AND project_id = ?"
    val = (str(sample_name), str(sample_description), species_id, project_id)
    cursor.execute(sql, val)
    res = cursor.fetchone()
    if res is None:
        sql = "INSERT INTO samples (sample_name, sample_description, species_id, project_id) VALUES (?, ?, ?, ?)"
        val = (str(sample_name), str(sample_description), species_id, project_id)
        cursor.execute(sql, val)
        conn.commit()
        print("Sample record added!")
    else:
        print("Sample record already exists")
        sample_id = res[0]
        print(f"current sample id is : {sample_id}")



if __name__ == "__main__":
    if len(sys.argv) == 1:
        user_inputs = get_user_input()
        sample_name = user_inputs['sample_name']
        sample_description = user_inputs['sample_description']
        user_species_name = user_inputs['species_name']
        user_project_id = user_inputs['project_id']
        curr_path = os.getcwd()
        directory = "TopLib"
        toplib_filepath = os.path.join(curr_path, directory, "toplib.db")   
        if os.path.isfile(toplib_filepath):
            # Connecting to sqlite
            conn = sqlite3.connect(toplib_filepath) 
            project_id = get_projectID(user_project_id, conn)
            species_id = get_speciesID(user_species_name, conn)
            if project_id:
                sample_table_insert(sample_name, sample_description, species_id, project_id, conn)
            else:
                print('Current project ID not found!')
            conn.close()
        else:
            print("Database .db file does not exist in the TopLib folder.")
    else:
        print("Usage: python script.py")
        sys.exit()
        