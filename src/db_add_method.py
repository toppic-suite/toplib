import sys
import sqlite3 
import os
import pandas as pd


def get_user_input():    
    print("Please enter your instrument:\n1 = Q Exactive HF; 2 = Orbitrap Fusion Lumos; 3 = LTQ Orbitrap Elite.")       
    input_instrument = sys.stdin.readline()
    instrument_code = input_instrument.strip()
    if instrument_code == '1':
        instrument_name = 'Q Exactive HF'
    elif instrument_code == '2':
        instrument_name = 'Orbitrap Fusion Lumos'
    elif instrument_code == '3':
        instrument_name = 'LTQ Orbitrap Elite'
    else:
        print('No such instrument avaiable and please choose one of instruments: 1 or 2 or 3.')
        return
        
    print("Please enter your dissociation method: 1 = HCD; 2 = CID; 3 = ECD; 4 = EID; 5 = ETD.")    
    input_dissociation = sys.stdin.readline()
    dissociation_code = input_dissociation.strip()
    if dissociation_code == '1':
        dissociation_method = 'HCD'
    elif dissociation_code == '2':
        dissociation_method = 'CID'
    elif dissociation_code == '3':
        dissociation_method = 'ECD'
    elif dissociation_code == '4':
        dissociation_method = 'EID'
    elif dissociation_code == '5':
        dissociation_method = 'ETD'
    else:
        print("No such dissociation method avaiable and please choose one of methods: 1 or 2 or 3 or 4 or 5.")
        return
   
    print("Please enter your energy value (normalized collision energy): Ex., 0.2.")    
    input_energy = sys.stdin.readline()
    energy_value = input_energy.strip()
    energy_value = float(energy_value)
        
    print("Please enter your resolution:")    
    input_resolution = sys.stdin.readline()
    resolution = input_resolution.strip()
    resolution = float(resolution)
    
    user_para_dict = {
                    'instrument_name': instrument_name, 
                    'dissociation_method': dissociation_method,
                    'energy_value': energy_value, 'resolution': resolution
                    }
    return user_para_dict


def get_vendorID(conn):
    vendor_df = pd.read_sql(sql="SELECT * FROM vendors", con=conn)
    vendor_id = vendor_df['vendor_id'].values
    return vendor_id


def get_dissociationID(dissociation_method, conn):
    cursor = conn.cursor()            
    sql = "SELECT dissociation_id FROM dissociations WHERE dissociation_method = ?"
    val = (str(dissociation_method), )
    cursor.execute(sql, val)
    res = cursor.fetchone() 
    dissociation_id = res[0]
    return dissociation_id  
            

def get_instrumentID(instrument_name, conn):
    vendor_id = get_vendorID(conn)
    cursor = conn.cursor()
    sql = "SELECT instrument_id FROM instruments WHERE instrument_name = ? AND vendor_id = ?"
    val = (str(instrument_name), int(vendor_id))
    cursor.execute(sql, val)
    res = cursor.fetchone() 
    instrument_id = res[0]
    return instrument_id 


# insert a record to methods table
def methods_table_insert(instrument_id, diss_id, energy_val, resolution, conn):
    # check if the record already exists
    cursor = conn.cursor()
    sql = "SELECT * FROM methods WHERE instrument_id = ? AND dissociation_id = ? AND energy_value = ? AND resolution =?"
    val = (int(instrument_id), int(diss_id), energy_val, resolution)
    cursor.execute(sql, val)
    res = cursor.fetchone() 
    if res is None:
        sql = "INSERT INTO methods (instrument_id, dissociation_id, energy_value, resolution) VALUES (?, ?, ?, ?)"
        val = (int(instrument_id), int(diss_id), energy_val, resolution)
        cursor.execute(sql, val)
        conn.commit()
        print("Method record add!")
    else:
        print("Method record already exists")
        method_id = res[0]
        print(f"current method id is : {method_id}")
        


if __name__ == "__main__":
    if len(sys.argv) == 1:
        user_inputs = get_user_input()
        instrument_name = user_inputs['instrument_name']
        dissociation_method = user_inputs['dissociation_method']
        energy_value = user_inputs['energy_value']
        resolution = user_inputs['resolution']
        curr_path = os.getcwd()
        directory = "TopLib"
        toplib_filepath = os.path.join(curr_path, directory, "toplib.db")   
        if os.path.isfile(toplib_filepath):
            # Connecting to sqlite
            conn = sqlite3.connect(toplib_filepath)
            dissociation_id = get_dissociationID(dissociation_method, conn)
            instrument_id = get_instrumentID(instrument_name, conn)
            # insert record to methods table
            methods_table_insert(instrument_id, dissociation_id, energy_value, resolution, conn)
            conn.close()                  
        else:
            print("Database .db file does not exist in the TopLib folder.")
    else:
        print("Usage: python script.py")
        sys.exit()
        
        
        