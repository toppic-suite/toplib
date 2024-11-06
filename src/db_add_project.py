import sys
import os
import sqlite3 


def get_user_input():
    print("Please enter your project name:")    
    input_project_name = sys.stdin.readline()
    project_name = input_project_name.strip()
    print("Please enter your project description:")    
    input_project_des = sys.stdin.readline()
    project_description = input_project_des.strip() 
    user_para_dict = {'project_name': project_name,
                      'project_description': project_description}
    return user_para_dict


# insert a project record to projects table
def project_table_insert(project_name, project_des, conn):
    # check if the record already exists
    cursor = conn.cursor()
    sql = "SELECT * FROM projects WHERE project_name = ?"
    val = (project_name, )
    cursor.execute(sql, val)
    res = cursor.fetchone()
    # if the record does not exist, insert a new record
    if res is None:
        sql = "INSERT INTO projects (project_name, project_description) VALUES (?, ?)"
        val = (str(project_name), str(project_des))
        cursor.execute(sql, val)
        conn.commit()
        print("Project record added!")
    else:
        print("Project record already exists")
        project_id = res[0]
        print(f"current project id is : {project_id}")
       

        
if __name__ == "__main__":
    if len(sys.argv) == 1:
        user_inputs = get_user_input()
        project_name = user_inputs['project_name']
        project_description = user_inputs['project_description']
        curr_path = os.getcwd()
        directory = "TopLib"
        toplib_filepath = os.path.join(curr_path, directory, "toplib.db")   
        if os.path.isfile(toplib_filepath):
            # Connecting to sqlite
            conn = sqlite3.connect(toplib_filepath) 
            # insert items to projects table
            project_table_insert(project_name, project_description, conn)
            conn.close()
        else:
            print("Database .db file does not exist in the TopLib folder.")
    else:
        print("Usage: python script.py")
        sys.exit()
        
        