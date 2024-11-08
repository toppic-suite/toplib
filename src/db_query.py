import sys
import os
import argparse
import sqlite3 
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as shc 
from db_add_sample import get_projectID
from tqdm import tqdm


def user_command():
    parser = argparse.ArgumentParser(description="Extracting library spectra parameters specification")
    parser.add_argument("project_id", type=int,
                        help="A project id (ex., 1)")
    parser.add_argument("rep_method", type=str,
                        help="Type of representative spectra (average or single)")
    args = parser.parse_args()
    project_id = int(args.project_id)
    rep_method = args.rep_method.lower()
    return project_id, rep_method 


def get_sampleID_by_projectID(project_id, conn):
    query = "SELECT * FROM samples WHERE project_id= '%s'" % project_id
    sample_df = pd.read_sql_query(sql=query, con=conn)
    sample_ids = sample_df['sample_id'].values
    return sample_ids 


def get_experimentID_by_sampleID(sample_ids, conn):
    sample_ids_str = ', '.join(map(str, sample_ids))
    query = "SELECT * FROM experiments WHERE sample_id IN ({})".format(sample_ids_str)    
    experiment_df = pd.read_sql_query(sql=query, con=conn) 
    exp_ids = experiment_df['experiment_id'].values
    return exp_ids


def get_fileID_by_experimentID(experiment_ids, conn):
    experiment_ids_str = ', '.join(map(str, experiment_ids))
    query = "SELECT * FROM files WHERE experiment_id IN ({})".format(experiment_ids_str)    
    file_df = pd.read_sql_query(sql=query, con=conn) 
    file_ids = file_df['file_id'].values
    return file_ids


def get_spectra_data(file_id, conn):
    file_id = ', '.join(map(str, file_id))
    query = "SELECT * FROM spectra WHERE file_id IN ({})".format(file_id)        
    spectra_df = pd.read_sql_query(sql=query, con=conn)
    return spectra_df


def get_fileName_by_fileID(file_id, conn):
    file_id = ', '.join(map(str, file_id))
    query = "SELECT * FROM files WHERE file_id IN ({})".format(file_id)        
    file_df = pd.read_sql_query(sql=query, con=conn)
    file_name = file_df['file_name'].values
    return file_name


def get_rep_by_fileID(file_ids, consensus_method, conn):    
    file_ids_str = ', '.join(map(str, np.array(file_ids)))
    if consensus_method == 'single':
        query = "SELECT * FROM single_representatives WHERE file_id IN ({})".format(file_ids_str)        
        rep_df = pd.read_sql_query(sql=query, con=conn)    
    else:
        query = "SELECT * FROM average_representatives WHERE file_id IN ({})".format(file_ids_str)    
        rep_df = pd.read_sql_query(sql=query, con=conn)    
        
    # add proteoform information
    spectra_df = get_spectra_data(file_ids, conn)
    file_names = get_fileName_by_fileID(file_ids, conn)
    file_id_to_name = dict(zip(file_ids, file_names))
    spectra_df['file_name'] = spectra_df['file_id'].map(file_id_to_name)    
    spectra_df = spectra_df[spectra_df['spectrum_id'].isin(rep_df['spectrum_id'].values)]
    spectra_df = spectra_df.drop(['file_id','precursor_mass','precursor_intensity','precursor_charge'], axis=1)    
    rep_df = spectra_df.merge(rep_df, on='spectrum_id', how='inner')
    # update proteoform ID
    rep_df_new = proteoformID_update(rep_df)
    return rep_df_new


def proteoformID_update(rep_df, precursor_tol=2.2):
    # update proteoform id for combined representatives
    rep_ext_df = rep_df.dropna(subset='proteoform')
    start_id = 0
    for protein in rep_ext_df['protein_accession'].unique():
        sub_df = rep_ext_df[rep_ext_df['protein_accession'] == protein]
        sub_idx = sub_df.index.values
        mz = sub_df['precursor_mass'].values
        if len(mz)>1:
            c_mz = mz[:, None] # column vectors
            Y = pdist(c_mz, 'cityblock')
            Z = shc.linkage(Y, method='complete')
            labels = shc.fcluster(Z, t=precursor_tol, criterion='distance')-1        
        else:
            labels = np.array([0])         
        new_ids = labels + start_id
        start_id = max(new_ids) + 1
        rep_df.loc[sub_idx,'proteoform_id'] = new_ids
    return rep_df


def rep2msalign(rep_df, msalign_wfile):
    curr_path = os.getcwd()
    directory = "TopLib"
    wfile = os.path.join(curr_path, directory, msalign_wfile)  
    # extract spectral data and convert to a msalign file
    columns1 = list(map(lambda x: 'mass_' + str(x), np.arange(50)))
    columns2 = list(map(lambda x: 'intensity_' + str(x), np.arange(50)))
    columns3 = list(map(lambda x: 'charge_' + str(x), np.arange(50))) 
    # writing to a msalign file
    with tqdm(total=len(rep_df), desc="Processing") as pbar:
        with open(wfile, 'w') as as_file:
            for ss in range(len(rep_df)):
                pbar.update(1)
                spectrumID_lib = rep_df['spectrum_id'].iloc[ss]
                scans_lib = rep_df['scan'].iloc[ss]
                filename_lib = rep_df['file_name'].iloc[ss]
                retention_time_lib = rep_df['retention_time'].iloc[ss]
                preMass_lib = rep_df['precursor_mass'].iloc[ss]
                preInte_lib = rep_df['precursor_intensity'].iloc[ss]
                preCharge_lib = rep_df['precursor_charge'].iloc[ss]
                pre_featureID_lib = rep_df['precursor_feature_id'].iloc[ss]
                mass_lib = rep_df[columns1].iloc[ss].values
                mass_lib = mass_lib[~np.isnan(mass_lib)]
                inte_lib = rep_df[columns2].iloc[ss].values
                inte_lib = inte_lib[~np.isnan(inte_lib)]
                ch_lib = rep_df[columns3].iloc[ss].values
                ch_lib = ch_lib[~np.isnan(ch_lib)]
                ch_lib = ch_lib.astype('int32')
                # generate a dataframe
                data = {'lib_mass': mass_lib,
                        'lib_intensity': inte_lib,
                        'lib_charge': ch_lib}
                df = pd.DataFrame(data)
                # write in
                as_file.write('BEGIN IONS' + "\n")
                as_file.write("FILE_NAME={}\n".format(filename_lib))
                as_file.write("SPECTRUM_ID={}\n".format(spectrumID_lib))
                as_file.write("SCANS={}\n".format(scans_lib))
                as_file.write("RETENTION_TIME={}\n".format(retention_time_lib))
                as_file.write("PRECURSOR_MASS={:.5f}\n".format(preMass_lib))
                as_file.write("PRECURSOR_INTENSITY={:.2f}\n".format(preInte_lib))
                as_file.write("PRECURSOR_CHARGE={}\n".format(preCharge_lib))
                as_file.write("PRECURSOR_FEATURE_ID={}\n".format(pre_featureID_lib))
                for m in range(len(df)):
                    as_file.write("{:.5f}\t{:.2f}\t{}\n".format(df['lib_mass'].values[m], df['lib_intensity'].values[m], df['lib_charge'].values[m]))
                as_file.write('END IONS' + "\n")
                as_file.write("\n")            

            
def rep2tsv(rep_df, tsv_wfile):
    rep_df = rep_df.dropna(subset='proteoform')
    curr_path = os.getcwd()
    directory = "TopLib"
    wfile = os.path.join(curr_path, directory, tsv_wfile)   
    columns = [f"{prefix}_{i}" for prefix in ['mass', 'intensity', 'charge', 'norm_intensity'] for i in range(50)]
    rep_df = rep_df.drop(columns=columns, axis=1)  
    rep_df = rep_df.drop(columns=['representative_id','file_id'], axis=1)  
    # rename columns name and keep the same column's name as tsv file obtained from toppic
    rep_df = rep_df.rename(columns={
        'spectrum_id': 'Spectrum ID',
        'proteoform_id': 'Proteoform ID',
        'protein_accession': 'Protein accession',
        'protein_description': 'Protein description',
        'first_residue': 'First residue',
        'last_residue': 'Last residue',
        'protein_sequence': 'Database protein sequence',
        'proteoform': 'Proteoform',
        'e_value': 'E-value',
        'spectrum_level_q_value': 'Spectrum-level Q-value',
        'proteoform_level_q_value': 'Proteoform-level Q-value'
        }) 
    # reorder
    col_order = ['Spectrum ID','file_name','scan', 'retention_time','precursor_mass', 'precursor_intensity', 'precursor_charge',
                  'precursor_feature_id', 'Proteoform ID', 'Protein accession', 'Protein description', 'First residue', 'Last residue',
                  'Database protein sequence', 'Proteoform', 'E-value', 'Spectrum-level Q-value', 'Proteoform-level Q-value']
    rep_df = rep_df[col_order]             
    rep_df.to_csv(wfile, sep='\t', index=False)    
    


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <project_id> <representative_method>")
        sys.exit()
    else:
        user_project_id, rep_method = user_command()
        curr_path = os.getcwd()
        directory = "TopLib"
        toplib_filepath = os.path.join(curr_path, directory, "toplib.db")   
        if os.path.isfile(toplib_filepath):
            # Connecting to sqlite
            conn = sqlite3.connect(toplib_filepath)
            # check user's input
            project_id = get_projectID(user_project_id, conn)
            if project_id:
                sample_ids = get_sampleID_by_projectID(project_id, conn)
                experiment_ids = get_experimentID_by_sampleID(sample_ids, conn)
                file_ids = get_fileID_by_experimentID(experiment_ids, conn)
                rep_df = get_rep_by_fileID(file_ids, rep_method, conn) 
                #convert to msalign file
                msalign_wfile = 'lib_spectra_ms2.msalign'
                rep2msalign(rep_df, msalign_wfile)
                print(f'msalign file: {msalign_wfile} generated!')
                # convert to tsv file
                tsv_wfile = 'lib_spectra_ms2_toppic_prsm_single_filtered.tsv'
                rep2tsv(rep_df, tsv_wfile)
                print(f'tsv file: {tsv_wfile} generated!')
            else:
                print('No such project ID and please check your project ID!')
            conn.close()  
        else:
            print("Database .db file does not exist in the TopLib folder.")

