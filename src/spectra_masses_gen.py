import pandas as pd
import numpy as np
import msalign_file as mf


def get_spectra_masses_tables(file_id, msalign_filename, prsm_filename, conn):
    # insert items to spectra table
    def spectra_table_insert(file_id, data, conn):
        cursor = conn.cursor()
        sql = "SELECT * FROM spectra WHERE file_id = ?"
        val = (int(file_id), )
        cursor.execute(sql, val)
        # if the record does not exist, insert a new record
        if cursor.fetchone() is None:
            print('Adding spectra data...')
            for row in data.itertuples():
                cursor.execute('''
                INSERT INTO spectra (file_id, scan, retention_time, precursor_mass, precursor_intensity, precursor_charge, num_mass, precursor_feature_id, proteoform_id, cluster_id, protein_accession, first_residue, last_residue, proteoform, protein_sequence, protein_description, e_value, spectrum_level_q_value, proteoform_level_q_value)
                VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                ''',
                (row.file_id,
                row.scan,
                row.retention_time,
                row.precursor_mass,
                row.precursor_intensity,
                row.precursor_charge,
                row.num_mass,
                row.precursor_feature_id,
                row.proteoform_id, 
                row.cluster_id,
                row.protein_accession,
                row.first_residue,
                row.last_residue,
                row.proteoform,
                row.protein_sequence,
                row.protein_description,
                row.e_value,
                row.spectrum_level_q_value,
                row.proteoform_level_q_value)
                )
            conn.commit()
            print("Spectra record added successfully")
            spectra_flag = 1
        else:
            print("Spectra record already exists")
            spectra_flag = 0
        return spectra_flag
    
    
    # get ids of spectra table
    def get_specID(file_id, conn):
        cursor = conn.cursor()
        # get spectrum_id
        sql = "SELECT spectrum_id FROM spectra WHERE file_id=?"
        val = (file_id, ) 
        cursor.execute(sql, val)
        fid_record = cursor.fetchall()
        spectrum_id_start = fid_record[0][0] - 1
        return spectrum_id_start
    
    
    # add the protein accession, first and last residue and proteoform
    def get_proteoform(prsm_filename, ms_df):
        tsv_df = pd.read_csv(prsm_filename, sep='\t')
        idx_ext_all = []
        tsv_idx_all = []
        for i in range(len(tsv_df)):
            spec_id = tsv_df['Spectrum ID'].iloc[i]
            idx_ext = ms_df[ms_df['spectrum_id']==spec_id].index.values
            if len(idx_ext):
                idx_ext_all.extend(idx_ext)
                tsv_idx_all.append(i)  
        
        ms_df_new = ms_df.copy()
        ms_df_new.loc[idx_ext_all,'proteoform_id'] = tsv_df['Proteoform ID'].iloc[tsv_idx_all].values
        ms_df_new.loc[idx_ext_all,'protein_accession'] = tsv_df['Protein accession'].iloc[tsv_idx_all].values
        ms_df_new.loc[idx_ext_all,'protein_description'] = tsv_df['Protein description'].iloc[tsv_idx_all].values
        ms_df_new.loc[idx_ext_all,'first_residue'] = tsv_df['First residue'].iloc[tsv_idx_all].values
        ms_df_new.loc[idx_ext_all,'last_residue'] = tsv_df['Last residue'].iloc[tsv_idx_all].values
        ms_df_new.loc[idx_ext_all,'protein_sequence'] = tsv_df['Database protein sequence'].iloc[tsv_idx_all].values 
        ms_df_new.loc[idx_ext_all,'proteoform'] = tsv_df['Proteoform'].iloc[tsv_idx_all].values
        ms_df_new.loc[idx_ext_all,'e_value'] = tsv_df['E-value'].iloc[tsv_idx_all].values
        ms_df_new.loc[idx_ext_all,'spectrum_level_q_value'] = tsv_df['Spectrum-level Q-value'].iloc[tsv_idx_all].values
        ms_df_new.loc[idx_ext_all,'proteoform_level_q_value'] = tsv_df['Proteoform-level Q-value'].iloc[tsv_idx_all].values
        return ms_df_new    
    
    
    # insert a new mass file into an existing mass table
    def masses_table_insert(data, conn):
        cursor = conn.cursor()
        # insert the value and run the query
        print('Adding fragment masses data...')
        for row in data.itertuples():
            cursor.execute('''
            INSERT INTO masses (spectrum_id, mass, intensity, charge)
            VALUES(?, ?, ?, ?);
            ''',
            (row.spectrum_id,
            row.mass,
            row.intensity,
            row.charge)
            )
        conn.commit()
        print("Mass record inserted successfully")
    
    print('Reading msalign file...')
    # generate the spectra and masses data from msalign file
    ms_df, out_f2 = mf.read_msalign(msalign_filename)
    ms_df = ms_df[ms_df['spectrum_id'].isin(out_f2['spectrum_id'].values)]
    # add file id to spectra
    ms_df['file_id'] = file_id
    # add cluster flag to spectra
    ms_df['cluster_flag'] = np.zeros((len(ms_df),), dtype=int) 
    ms_df['cluster_id'] = None
   
    # add proteoform information
    out_f1 = get_proteoform(prsm_filename, ms_df)
    out_f1 = out_f1.replace(np.nan, None) # convert to none
    
    # insert spectra table
    spectra_flag = spectra_table_insert(file_id, out_f1, conn)
    # update the spectrum id    
    spec_id_start = get_specID(file_id, conn)
    out_f2['Spectrum ID'] = out_f2.groupby(['spectrum_id']).ngroup() + 1
    out_f2 = out_f2.drop(['spectrum_id'], axis=1) 
    out_f2.rename(columns = {'Spectrum ID':'spectrum_id'}, inplace = True)
    out_f2.loc[:,'spectrum_id'] = out_f2.loc[:, 'spectrum_id'] + spec_id_start
    # insert masses table
    if spectra_flag:
        masses_table_insert(out_f2, conn)
    
