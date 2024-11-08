import sys
import os
import sqlite3 
import db_add_species_diss_instru_vendor as add_sdiv


def db_tables_gen(filepath):
    # generate tables in database
    print("Generating a toplib library...")
    # Connecting to sqlite
    conn = sqlite3.connect(filepath) 
    cursor = conn.cursor()
    # create metadata tables 
    # create project table
    cursor.execute("CREATE TABLE IF NOT EXISTS projects (project_id INTEGER PRIMARY KEY AUTOINCREMENT, project_name VARCHAR(50), project_description VARCHAR(50))")
    # create species table
    cursor.execute("CREATE TABLE IF NOT EXISTS species (species_id INTEGER PRIMARY KEY AUTOINCREMENT, species_name VARCHAR(50), taxon_id INTEGER)")
    # create sample table
    cursor.execute("CREATE TABLE IF NOT EXISTS samples (sample_id INTEGER PRIMARY KEY AUTOINCREMENT, sample_name VARCHAR(50), sample_description VARCHAR(50), species_id INTEGER, project_id INTEGER, FOREIGN KEY(species_id) REFERENCES species(species_id) ON DELETE CASCADE, FOREIGN KEY(project_id) REFERENCES projects(project_id) ON DELETE CASCADE)")
    # create dissociations table
    cursor.execute("CREATE TABLE IF NOT EXISTS dissociations (dissociation_id INTEGER PRIMARY KEY AUTOINCREMENT, dissociation_method VARCHAR(50))")
    # create vendors table
    cursor.execute("CREATE TABLE IF NOT EXISTS vendors (vendor_id INTEGER PRIMARY KEY AUTOINCREMENT, vendor_name VARCHAR(50))")
    # create instruments table 
    cursor.execute("CREATE TABLE IF NOT EXISTS instruments (instrument_id INTEGER PRIMARY KEY AUTOINCREMENT, vendor_id INTEGER, instrument_name VARCHAR(50), FOREIGN KEY(vendor_id) REFERENCES vendors(vendor_id) ON DELETE CASCADE)")
    # create methods table
    cursor.execute("CREATE TABLE IF NOT EXISTS methods (method_id INTEGER PRIMARY KEY AUTOINCREMENT, instrument_id INTEGER, dissociation_id INTEGER, energy_value DOUBLE, resolution DOUBLE, FOREIGN KEY(instrument_id) REFERENCES instruments(instrument_id) ON DELETE CASCADE, FOREIGN KEY(dissociation_id) REFERENCES dissociations(dissociation_id) ON DELETE CASCADE)")
    # create experiment table
    cursor.execute("CREATE TABLE IF NOT EXISTS experiments (experiment_id INTEGER PRIMARY KEY AUTOINCREMENT, experiment_name VARCHAR(50), experiment_description VARCHAR(50), sample_id INTEGER, method_id INTEGER, FOREIGN KEY(method_id) REFERENCES methods(method_id) ON DELETE CASCADE, FOREIGN KEY(sample_id) REFERENCES samples(sample_id) ON DELETE CASCADE)")
    # create file table
    cursor.execute("CREATE TABLE IF NOT EXISTS files (file_id INTEGER PRIMARY KEY AUTOINCREMENT, file_name VARCHAR(50), experiment_id INTEGER, cluster_flag INTEGER, single_representative VARCHAR(50), average_representative VARCHAR(50), FOREIGN KEY(experiment_id) REFERENCES experiments(experiment_id) ON DELETE CASCADE)")
    
    # create spectrum tables
    # create spectra table
    cursor.execute("CREATE TABLE IF NOT EXISTS spectra (spectrum_id INTEGER PRIMARY KEY AUTOINCREMENT, file_id INTEGER, scan INTEGER, retention_time DOUBLE, precursor_mass DOUBLE, precursor_intensity DOUBLE, precursor_charge INTEGER, num_mass INTEGER, precursor_feature_id INTEGER, proteoform_id INTEGER, cluster_id INTEGER, protein_accession VARCHAR(50), first_residue INTEGER, last_residue INTEGER, proteoform VARCHAR(1024), protein_sequence VARCHAR(1024), protein_description VARCHAR(1024), e_value DOUBLE, spectrum_level_q_value DOUBLE, proteoform_level_q_value DOUBLE, FOREIGN KEY(file_id) REFERENCES files(file_id) ON DELETE CASCADE)")
    # create masses table 
    cursor.execute("CREATE TABLE IF NOT EXISTS masses (mass_id INTEGER PRIMARY KEY AUTOINCREMENT, spectrum_id INTEGER, mass DOUBLE, intensity DOUBLE, charge INTEGER, FOREIGN KEY(spectrum_id) REFERENCES spectra(spectrum_id) ON DELETE CASCADE)")
    
    species_taxon_lst = {'common_name': ['human','mouse'],
                         'taxon_id': [9606,10090]}
    dissociation_lst = ['HCD','CID','ECD','EID','ETD']
    instrument_lst = ['Q Exactive HF','Orbitrap Fusion Lumos','LTQ Orbitrap Elite']
    vendor_lst = ['Thermo']
    lib_parameters = {
                'species_taxon_id': species_taxon_lst,
                'instrument_name': instrument_lst,
                'dissociation_method': dissociation_lst,
                'vendor_name': vendor_lst
                }
    add_sdiv.species_diss_instru_vendor_tables(lib_parameters, conn)    
    conn.close() 
    

def create_database(db_name):        
    # create a new database
    curr_path = os.getcwd()
    directory = "TopLib"
    path = os.path.join(curr_path, directory)
    try:
        os.makedirs(path, exist_ok = True)
        print("Directory '%s' created successfully" % directory)
        filepath = os.path.join(curr_path, directory, db_name)   
        # create a ms library
        db_tables_gen(filepath)     
        
    except OSError as error:
        print("Directory '%s' can not be created" % directory)


if __name__ == "__main__":   
    if len(sys.argv)==1:
        db_names = ['toplib.db']
        for db_name in db_names:
            create_database(db_name)
            print(db_name + ' has been created!')
    else:
        print("Usage: python script.py")
        sys.exit()
        

            
    
