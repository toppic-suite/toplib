import sys
import os
import argparse
import pandas as pd
import numpy as np
import numba as nb
# import random
import sqlite3 
from tqdm import tqdm
import time
from datetime import datetime


def get_db_data(lib_filename):
    columns_mass = [f"mass_{i}" for i in range(50)]
    # columns_inte = [f"intensity_{i}" for i in range(50)]
    columns_ch   = [f"charge_{i}" for i in range(50)]
    columns_norm = [f"norm_intensity_{i}" for i in range(50)]
    
    meta_cols = [
        'representative_id', 
        'precursor_mass', 
        'precursor_intensity', 
        'precursor_charge', 
        'dataset_id', 
        'mzml_file_name', 
        'mzml_ms2_scan', 
        # 'mass', 'intensity', 
        # 'charge', 
        # 'mzml_instrument', 
        'mzml_ms2_retention_time', 
        'mzml_activation', 
        'mzml_collision_energy', 
        'msalign_file_name', 
        'msalign_feature_intensity', 
        'msalign_feature_score', 
        'msalign_feature_apex_time', 
        'msalign_number_of_fragment_ions', 
        'toppic_prsm_id', 
        'toppic_adjusted_precursor_mass', 
        'toppic_proteoform_id', 
        'toppic_proteoform_intensity', 
        'toppic_number_of_protein_hits', 
        'toppic_protein_accession', 
        'toppic_protein_description', 
        'toppic_first_residue_position', 
        'toppic_last_residue_position', 
        'toppic_special_amino_acids',
        'toppic_database_sequence', 
        'toppic_proteoform_mass',
        'toppic_protein_n_terminal_form',
        'toppic_fixed_ptms', 
        'toppic_number_of_unexpected_modifications',
        'toppic_unexpected_modifications',
        'toppic_number_of_variable_modifications',
        'toppic_variable_modifications',
        'toppic_miscore',
        'toppic_number_of_matched_experimental_fragment_ions',
        'toppic_number_of_matched_theoretical_fragment_masses',
        'toppic_e_value', 
        'toppic_spectrum_level_q_value', 
        'toppic_proteoform_level_q_value',
        'toppic_proteoform',
        'toppic_previous_residue',
        'toppic_next_residue',
        'proteoform_id',
        'flag'
        ]
    
    select_cols = meta_cols + columns_mass + columns_ch + columns_norm
    
    query = f"""
    SELECT {",".join(select_cols)}
    FROM target_decoy_spectra_representatives
    """
    
    conn = sqlite3.connect(lib_filename)
    
    chunksize = 5000   
    
    chunk_iter = pd.read_sql_query(
        query,
        con=conn,
        chunksize=chunksize
    )
    
    
    mass_chunks = []
    # inte_chunks = []
    ch_chunks   = []
    norm_chunks = []
    
    for chunk_id, chunk in enumerate(chunk_iter, start=1):
        if chunk_id % 10 == 0:
            print(f"Processing chunk {chunk_id}")
    
        mass_df = chunk[meta_cols + columns_mass]
        # inte_df = chunk[columns_inte]
        ch_df   = chunk[columns_ch]
        norm_df = chunk[columns_norm]
    
        
        mass_chunks.append(mass_df)
        # inte_chunks.append(inte_df)
        ch_chunks.append(ch_df)
        norm_chunks.append(norm_df)
    
    conn.close()
    
    target_decoy_mass_rep = pd.concat(mass_chunks, ignore_index=True)
    # target_decoy_inte_rep = pd.concat(inte_chunks, ignore_index=True)
    target_decoy_ch_rep   = pd.concat(ch_chunks, ignore_index=True)
    target_decoy_norm_rep = pd.concat(norm_chunks, ignore_index=True)
    print(f"total rows: {len(target_decoy_mass_rep)}")
    
    return target_decoy_mass_rep, target_decoy_ch_rep, target_decoy_norm_rep


def tsv_file_preprocessing(tsv_filename):
    df = pd.read_csv(tsv_filename, sep='\t', low_memory=False)   
    # df = rename_cols(df)
    sel_cols =['DATASET_id', 'MZML_instrument', 'MZML_file_name', 'MZML_ms2_scan', 'MZML_ms2_retention_time', 'MZML_activation', 'MZML_collision_energy',
               'MSALIGN_file_name', 'MSALIGN_precursor_charge','MSALIGN_precursor_monoisotopic_mass', 'MSALIGN_precursor_intensity',
               'MSALIGN_feature_id', 'MSALIGN_feature_intensity',
               'MSALIGN_feature_score', 'MSALIGN_feature_apex_time',
               'MSALIGN_number_of_fragment_ions', 'TOPPIC_prsm_id',
               'TOPPIC_adjusted_precursor_mass', 'TOPPIC_proteoform_id',
               'TOPPIC_proteoform_intensity', 'TOPPIC_number_of_protein_hits',
               'TOPPIC_protein_accession', 'TOPPIC_protein_description',
               'TOPPIC_first_residue_position', 'TOPPIC_last_residue_position',
               'TOPPIC_special_amino_acids', 'TOPPIC_database_sequence',
               'TOPPIC_proteoform_mass', 'TOPPIC_protein_n-terminal_form',
               'TOPPIC_fixed_modifications', 'TOPPIC_number_of_unexpected_modifications',
               'TOPPIC_unexpected_modifications', 'TOPPIC_number_of_variable_modifications',
               'TOPPIC_variable_modifications', 'TOPPIC_miscore',
               'TOPPIC_number_of_matched_experimental_fragment_ions',
               'TOPPIC_number_of_matched_theoretical_fragment_masses',
               'TOPPIC_e-value', 'TOPPIC_spectrum-level_q-value',
               'TOPPIC_proteoform-level_q-value', 'TOPPIC_proteoform',
               'TOPPIC_previous_residue', 'TOPPIC_next_residue']
               
    df_ext = df[sel_cols]
    return df_ext


def load_msalign_data(msalign_file):
    rows = []
    dataset_id = None
    mzml_filename = None
    scan = None

    mass_all = []
    intensity_all = []
    charge_all = []
    confidence_all = []

    with open(msalign_file, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()

            if line == "BEGIN IONS":
                # reset for new spectrum
                dataset_id = None
                mzml_filename = None
                scan = None

                mass_all = []
                intensity_all = []
                charge_all = []
                confidence_all = []

            elif line.startswith("DATASET_ID="):
                dataset_id = line.split("=", 1)[1]

            elif line.startswith("MZML_FILE_NAME="):
                mzml_filename = line.split("=", 1)[1]

            elif line.startswith("MS2_SCAN="):
                scan = int(line.split("=", 1)[1])

            elif line == "END IONS":
                if scan is not None:
                    rows.append({
                        "dataset_id": dataset_id,
                        "mzml_file_name": mzml_filename,
                        "scan": scan,
                        "mass_all": mass_all,
                        "intensity_all": intensity_all,
                        "charge_all": charge_all,
                        "confidence_all": confidence_all
                    })

                # clean up explicitly
                dataset_id = None
                mzml_filename = None
                scan = None

            else:
                # Peak line
                parts = line.split()
                if len(parts) >= 4:
                    # mandatory fields
                    mass_all.append(float(parts[0]))
                    intensity_all.append(float(parts[1]))
                    charge_all.append(int(parts[2]))
                    confidence_all.append(float(parts[3]))
                   
    msalign_df = pd.DataFrame(rows)
    rename_map = {
        "dataset_id": "DATASET_id",
        "mzml_file_name": "MZML_file_name",
        "scan": "MZML_ms2_scan",
        "mass_all": "mass",
        "intensity_all": "intensity",
        "charge_all": "charge",
        "confidence_all": "confidence"
    }
    msalign_df = msalign_df.rename(columns=rename_map)
    
    return msalign_df


def first_precursor_mass_extract(df):
    df['precursor_mass'] = df['MSALIGN_precursor_monoisotopic_mass'].apply(lambda x: x.split(':')[0])
    df['precursor_mass'] = df['precursor_mass'].astype(float)
     
    df['precursor_intensity'] = df['MSALIGN_precursor_intensity'].apply(lambda x: x.split(':')[0])
    df['precursor_intensity'] = df['precursor_intensity'].astype(float)
    
    df['precursor_charge'] = df['MSALIGN_precursor_charge'].apply(lambda x: x.split(':')[0])
    df['precursor_charge'] = df['precursor_charge'].astype(int)
    
    return df


def tsv_msalign_merge(msalign_df, tsv_df):
    msalign_df_merged = msalign_df.merge(
        tsv_df,
        on=['DATASET_id', 'MZML_file_name', 'MZML_ms2_scan'],
        how='left',
        suffixes=('', '_new')
    )
    msalign_df_merged = first_precursor_mass_extract(msalign_df_merged)
    return msalign_df_merged


def peak_list_extract_based_precursor(ms_df):
    for i in range(len(ms_df)):
        precursor_mass = ms_df.at[i, 'precursor_mass']

        mass = np.asarray(ms_df.at[i, 'mass'])
        mask = mass <= precursor_mass

        ms_df.at[i, 'mass'] = mass[mask].tolist()
        ms_df.at[i, 'intensity'] = np.asarray(ms_df.at[i, 'intensity'])[mask].tolist()
        ms_df.at[i, 'charge'] = np.asarray(ms_df.at[i, 'charge'])[mask].tolist()
        ms_df.at[i, 'confidence'] = np.asarray(ms_df.at[i, 'confidence'])[mask].tolist()

    return ms_df


def ms_incorr_remove(ms_df, num_pks=50):
    ms_org_df = ms_df.copy()
    mass_arr= ms_df['mass'].to_numpy()
    intensity_arr = ms_df['intensity'].to_numpy()
    charge_arr = ms_df['charge'].to_numpy()
    conf_arr = ms_df['confidence'].to_numpy()
    precursor_intensity = ms_df['precursor_intensity'].values  
    ms_index = ms_df.index.values
    num_spec = len(ms_df)
    valid_indices = []
    for i in range(num_spec):
        mass_list = np.array(mass_arr[i])
        inte_list = np.array(intensity_arr[i])
        charge_list = np.array(charge_arr[i])
        conf_list = np.array(conf_arr[i])
        # Sort by intensity, keep top peaks
        sorted_inte_idx = np.argsort(inte_list)[::-1][0:num_pks]
        mass_arr[i] = mass_list[sorted_inte_idx]
        intensity_arr[i] = inte_list[sorted_inte_idx]
        charge_arr[i] = charge_list[sorted_inte_idx]
        conf_arr[i] = conf_list[sorted_inte_idx]
        # Filter out spectra with <=1 peak and zero precursor intensity
        if len(mass_list[sorted_inte_idx]) > 1 and precursor_intensity[i] != 0:
            valid_indices.append(i)
        
    valid_ms_index = ms_index[valid_indices]
    ms_ext_df = ms_df.loc[valid_ms_index,:]
    ms_org_ext_df = ms_org_df.loc[valid_ms_index,:]
    ms_ext_df['mass'] = mass_arr[valid_indices]
    ms_ext_df['intensity'] = intensity_arr[valid_indices]
    ms_ext_df['charge'] = charge_arr[valid_indices]
    ms_ext_df['confidence'] = conf_arr[valid_indices]

    ms_ext_df.reset_index(drop=True, inplace=True)
    ms_org_ext_df.reset_index(drop=True, inplace=True)
    return ms_ext_df, ms_org_ext_df


def ms_spectrum_preprocess(ms_all):
    # sort, log2
    sorted_mass_all = []
    sorted_inte_all = []
    sorted_ch_all = []
    sorted_conf_all = []
    sorted_norm_inte_all = []
    for i in range(len(ms_all)):
        # sort
        sorted_mass_idx = np.argsort(ms_all['mass'].iloc[i])
        sorted_mass = np.array(ms_all['mass'].iloc[i])[sorted_mass_idx]
        sorted_inte = np.array(ms_all['intensity'].iloc[i])[sorted_mass_idx]
        # log and normalization
        sorted_log_inte = np.log2(sorted_inte)
        sorted_norm_inte = sorted_log_inte / np.linalg.norm(sorted_log_inte)
        sorted_ch = np.array(ms_all['charge'].iloc[i])[sorted_mass_idx]
        sorted_conf = np.array(ms_all['confidence'].iloc[i])[sorted_mass_idx]
        
        sorted_inte_all.append(sorted_inte)
        sorted_mass_all.append(sorted_mass)
        sorted_ch_all.append(sorted_ch)
        sorted_conf_all.append(sorted_conf)
        sorted_norm_inte_all.append(sorted_norm_inte)
        
    data = {'mass': sorted_mass_all,
            'intensity': sorted_inte_all,
            'charge': sorted_ch_all,
            'confidence': sorted_conf_all,
            'norm_intensity': sorted_norm_inte_all}
    sorted_mass_inte_ch_df = pd.DataFrame(data)
    return sorted_mass_inte_ch_df


@nb.njit
def peak_50_cosine_mass_ppm(sorted_mz1, sorted_int1, sorted_ch1, sorted_mz2, sorted_int2, sorted_ch2, err, ppm):
    cos_unembed = 0
    for i in range(len(sorted_mz1)):
        for j in range(len(sorted_mz2)):
            tol2 = max(err, ppm * sorted_mz1[i] * 1e-6)
            th1 = 0.0 + tol2
            th2 = 1.00235 + tol2
            th3 = 1.00235 - tol2
            if sorted_ch1[i] == sorted_ch2[j]:
                if abs(sorted_mz2[j] - sorted_mz1[i]) <= th1 or (
                        (abs(sorted_mz2[j] - sorted_mz1[i]) <= th2) and (abs(sorted_mz2[j] - sorted_mz1[i]) >= th3)):
                    cos_unembed += sorted_int1[i] * sorted_int2[j]
                    if cos_unembed > 1:
                        return 1.0 - (cos_unembed - sorted_int1[i] * sorted_int2[j])

    return 1.0 - cos_unembed


@nb.njit
def precursor_mass_query_Da(pre_mass_lib, pre_ch_lib, spectra_mass_rep_index, pre_mass_query_all, pre_ch_query_all, spectra_query_index, pre_tol_val, pre_ch_use):
    # spectra_mass_rep: spectra mass representative 
    # spectra_query: queried spectra
    idx1_all =[]
    idx2_all =[]
    for i in range(len(pre_mass_query_all)):
        pre_mass_query = pre_mass_query_all[i]
        pre_ch_query = pre_ch_query_all[i]
        idx1 = []
        idx2 = []
        for j in range(len(pre_mass_lib)):
            if pre_ch_use == 'yes':
                if pre_ch_query == pre_ch_lib[j]:
                    if abs(pre_mass_lib[j] - pre_mass_query) <= pre_tol_val:
                        idx1.append(spectra_mass_rep_index[j])
                        if spectra_query_index[i] not in idx2:
                            idx2.append(spectra_query_index[i])
            else:
                if abs(pre_mass_lib[j] - pre_mass_query) <= pre_tol_val:
                    idx1.append(spectra_mass_rep_index[j])
                    if spectra_query_index[i] not in idx2:
                        idx2.append(spectra_query_index[i])        

        if idx1:
            idx1_all.append(idx1)
            idx2_all.extend(idx2)
    return idx1_all, idx2_all


@nb.njit
def precursor_mass_query(pre_mass_lib,pre_ch_lib, spectra_mass_rep_index, pre_mass_query_all, pre_ch_query_all, spectra_query_index, ppm, pre_ch_use):
    """
    perform precursor mass and charge match
    """
    err = 0.0
    idx1_all =[]
    idx2_all =[]
    for i in range(len(pre_mass_query_all)):
        pre_mass_query = pre_mass_query_all[i]
        pre_ch_query = pre_ch_query_all[i]
        idx1 = []
        idx2 = []
        for j in range(len(pre_mass_lib)):
            tol2 = np.maximum(err, ppm*pre_mass_lib[j]*1e-6)
            th1 = 0.0+tol2
            # th2 = 1.00235+tol2
            # th3 = 1.00235-tol2
            if pre_ch_use == 'yes':
                if pre_ch_query == pre_ch_lib[j]:
                    if abs(pre_mass_lib[j] - pre_mass_query) <= th1:# or ((abs(pre_mass_lib[j] - pre_mass_query) <= th2) and (abs(pre_mass_lib[j] - pre_mass_query) >= th3)):
                        idx1.append(spectra_mass_rep_index[j])
                        if spectra_query_index[i] not in idx2:
                            idx2.append(spectra_query_index[i])
            else:
                if abs(pre_mass_lib[j] - pre_mass_query) <= th1:# or ((abs(pre_mass_lib[j] - pre_mass_query) <= th2) and (abs(pre_mass_lib[j] - pre_mass_query) >= th3)):
                    idx1.append(spectra_mass_rep_index[j])
                    if spectra_query_index[i] not in idx2:
                        idx2.append(spectra_query_index[i])

        if idx1:
            idx1_all.append(idx1)
            idx2_all.extend(idx2)
    return idx1_all, idx2_all


def dist_cal(data1_mass, data1_inte, data1_ch, data2_mass, data2_inte, data2_ch, thre, idx1, idx2, ppm):
    """
    calculate distances and group them to generate matched id and cosine distance
    """
    idx12 = []
    res_cos = []
    err = 0.0
    for i in range(len(idx2)):
        cos_all = []
        for j in range(len(idx1[i])):
            cos_v = peak_50_cosine_mass_ppm(data1_mass[idx1[i][j]], data1_inte[idx1[i][j]], data1_ch[idx1[i][j]], 
                                            data2_mass, data2_inte, data2_ch, err, ppm)
            cos_all.append(cos_v)
        temp_dis = np.array(cos_all)
        idx_min = np.argmin(temp_dis)
        if min(temp_dis)<= thre:
            idx12.append([idx1[i][idx_min], idx2[i]])
            res_cos.append(temp_dis[idx_min])
    return idx12, res_cos


def fdr_cal(combined_target_decoy_res):
    # calculate FDR
    ntarget_thre_all = []
    ndecoy_thre_all = []
    fdr_all = []
    spec_id_all = []
    thre_cutoff = np.arange(0.1,1.0, 0.001)
    for i in range(len(thre_cutoff)):
        flag_thre = combined_target_decoy_res.loc[combined_target_decoy_res['calculated_cosine_dist']<=thre_cutoff[i],'flag'].values
        spec_id = combined_target_decoy_res.loc[(combined_target_decoy_res['calculated_cosine_dist']<=thre_cutoff[i]),'spectrum_id'].values
        ntarget_thre = len([x for x in flag_thre if x==1])
        ndecoy_thre = len([x for x in flag_thre if x==0])
        if ntarget_thre ==0 and ndecoy_thre ==0:
            continue
        ntarget_thre_all.append(ntarget_thre)
        ndecoy_thre_all.append(ndecoy_thre)
        fdr = ndecoy_thre / (ntarget_thre + ndecoy_thre)
        fdr_all.append(fdr)
        spec_id_all.append(spec_id)

    # find 1% fdr    
    fdr_res = np.array(fdr_all)
    fdr_idx = np.argmin(abs(fdr_res*100-1.0))
    return spec_id_all[fdr_idx]  

def write_result_fun(wfile, df):
    """
    report the querying results to a tsv file
    """
    df.to_csv(wfile, sep='\t', index=False)
    # df.to_excel(wfile, index=False)
    
    
def ms_query(target_decoy_mass_rep, target_decoy_inte_rep_convert, target_decoy_ch_rep, 
             ms_df_query, mass_df_query, tol_type, pre_tol_val, frag_tol_val, pre_charge_use):
    # ms spectrum query
    columns1 = list(map(lambda x: 'mass_' + str(x), np.arange(50)))
    columns2 = list(map(lambda x: 'norm_intensity_' + str(x), np.arange(50)))
    columns3 = list(map(lambda x: 'charge_' + str(x), np.arange(50))) 
    num_spec = len(ms_df_query)
    spec_id_list = ms_df_query['spectrum_id'].values 
    thre = 1.0        
    res_cos_all = []
    idx12_all = []
    spectra_query_id_all = []
    spectra_query_filename_all = []
    # mass_query_all = []
    # intensity_query_all = []
    # charge_query_all = []
    flag_all = []
    with tqdm(total=num_spec, desc="Processing") as pbar:
        for i in range(num_spec):
            pbar.update(1)
            spectra_query = ms_df_query.loc[ms_df_query['spectrum_id']==spec_id_list[i],['spectrum_id','MZML_file_name','precursor_mass','precursor_intensity','precursor_charge']]
            spec_query_id = spectra_query['spectrum_id'].values[0]
            spec_query_filename = spectra_query['MZML_file_name'].values[0]
            spectra_query.reset_index(drop=True, inplace=True)
            mass_query_v = mass_df_query[mass_df_query['spectrum_id']==spec_id_list[i]]
            
            # step1: precursor match querying
            spectra_mass_rep_index = target_decoy_mass_rep.index.values
            pre_mass_lib = target_decoy_mass_rep['precursor_mass'].values
            pre_ch_lib = target_decoy_mass_rep['precursor_charge'].values

            spectra_query_index = spectra_query.index.values
            pre_mass_query_all = spectra_query['precursor_mass'].values
            pre_ch_query_all = spectra_query['precursor_charge'].values
            # print(mass_query_v)
            if tol_type == 'Da':
                idx1, idx2 = precursor_mass_query_Da(pre_mass_lib, pre_ch_lib, spectra_mass_rep_index, pre_mass_query_all, pre_ch_query_all, spectra_query_index, pre_tol_val, pre_charge_use)
            else:
                idx1, idx2 = precursor_mass_query(pre_mass_lib,pre_ch_lib, spectra_mass_rep_index, pre_mass_query_all, pre_ch_query_all,spectra_query_index, pre_tol_val, pre_charge_use)
            # step2: calculate distance
            if len(idx1):
                idx12, res_cos = dist_cal(target_decoy_mass_rep[columns1].to_numpy(), target_decoy_inte_rep_convert[columns2].to_numpy(), target_decoy_ch_rep[columns3].to_numpy(), 
                                    mass_query_v['mass'].values[0], mass_query_v['norm_intensity'].values[0], mass_query_v['charge'].values[0], thre, np.array(idx1), np.array(idx2), frag_tol_val)

                if idx12:
                    flag=target_decoy_mass_rep['flag'].iloc[idx12[0][0]]
                    flag_all.append(flag)
                    res_cos_all.extend(res_cos)
                    idx12_all.append(idx12[0][0])
                    spectra_query_id_all.append(spec_query_id) 
                    spectra_query_filename_all.append(spec_query_filename)
                    # mass_query_all.append(mass_query_v['mass'].values[0])
                    # intensity_query_all.append(mass_query_v['intensity'].values[0])
                    # charge_query_all.append(mass_query_v['charge'].values[0])
                   

    data = {'spectrum_id': spectra_query_id_all,
            'file_name': spectra_query_filename_all,
            'calculated_cosine_dist': res_cos_all,
            'lib_index': idx12_all,
            'flag': flag_all}
            # 'mass_query': mass_query_all,
            # 'intensity_query': intensity_query_all,
            # 'charge_query': charge_query_all}
    combined_target_decoy_res = pd.DataFrame(data)

    combined_target_decoy_res = combined_target_decoy_res.sort_values(by='spectrum_id')
    lib_index_all = combined_target_decoy_res['lib_index'].values
    combined_target_decoy_res['representative_id'] = target_decoy_mass_rep.loc[lib_index_all,'representative_id'].values
    combined_target_decoy_res['similarity'] = 1-combined_target_decoy_res['calculated_cosine_dist'].values
    combined_target_decoy_res['file_name(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'mzml_file_name'].values
#     combined_target_decoy_res['spectrum_id(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'spectrum_id'].values
    combined_target_decoy_res['scan(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'mzml_ms2_scan'].values
    combined_target_decoy_res['retention_time(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'mzml_ms2_retention_time'].values
    combined_target_decoy_res['proteoform_id(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'proteoform_id'].values
    combined_target_decoy_res['protein_accession(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_protein_accession'].values
    combined_target_decoy_res['precursor_mass(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'precursor_mass'].values
    combined_target_decoy_res['precursor_charge(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'precursor_charge'].values
#     combined_target_decoy_res['precursor_feature_id(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'precursor_feature_id'].values
    combined_target_decoy_res['protein_description(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_protein_description'].values
    combined_target_decoy_res['protein_sequence(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_database_sequence'].values
    combined_target_decoy_res['first_residue(lib)'] = target_decoy_mass_rep.loc[lib_index_all,'toppic_first_residue_position'].values
    combined_target_decoy_res['last_residue(lib)'] = target_decoy_mass_rep.loc[lib_index_all,'toppic_last_residue_position'].values
    combined_target_decoy_res['proteoform(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_proteoform'].values
    combined_target_decoy_res['e_value(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_e_value'].values
    combined_target_decoy_res['spectrum_level_q_value(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_spectrum_level_q_value'].values
    combined_target_decoy_res['proteoform_level_q_value(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_proteoform_level_q_value'].values
    combined_target_decoy_res['previous_residue(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_previous_residue'].values
    combined_target_decoy_res['next_residue(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_next_residue'].values
#     combined_target_decoy_res['sequence_length(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'toppic_sequence_length'].values
#     combined_target_decoy_res['paper_group_id(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'paper_group_id'].values
    return combined_target_decoy_res


def ms_library_search(lib_file, tsv_query_file, msalign_query, wfile):
    target_decoy_ms_rep1, target_decoy_ch_rep, target_decoy_inte_rep_convert = get_db_data(lib_file) 
    
    query_tsv_df = tsv_file_preprocessing(tsv_query_file)
    query_msalign_df = load_msalign_data(msalign_query)
    
    query_msalign_df1 = query_msalign_df.copy()
    query_msalign_df1['spectrum_id'] = query_msalign_df1.groupby(['MZML_file_name','MZML_ms2_scan']).ngroup()
    
    index_all = query_msalign_df1.index.values
    empty_idx = []
    for i in range(len(query_msalign_df1)):
        if len(query_msalign_df1['mass'].iloc[i])==0:
            empty_idx.append(index_all[i])
    
    query_msalign_df1_filtered = query_msalign_df1.drop(index=empty_idx)
    query_msalign_df1_filtered.reset_index(drop=True, inplace=True)
    
    query_ms_df = tsv_msalign_merge(query_msalign_df1_filtered, query_tsv_df)
    print(f"ms_df : {len(query_ms_df)}")
    
    # remove incorrect spectra
    drop_idx2 = query_ms_df[query_ms_df['precursor_intensity']==0].index.values
    query_ms_df.drop(drop_idx2, inplace = True)
    query_ms_df.reset_index(drop=True, inplace=True)
    drop_idx = query_ms_df[query_ms_df['MSALIGN_number_of_fragment_ions']<=1].index.values
    query_ms_df.drop(drop_idx, inplace = True)
    query_ms_df.reset_index(drop=True, inplace=True)
    
    query_ms_df = query_ms_df.rename(columns={"Project_ID": 'DATASET_id',
                               "MZML_MS2_scan": 'MZML_ms2_scan',
                               "MZML_MS2_retention_time": "MZML_ms2_retention_time",
                               "MSALIGN_feature_ID": "MSALIGN_feature_id",
                               "TOPPIC_PrSM_ID": "TOPPIC_prsm_id",
                               "TOPPIC_proteoform_ID": "TOPPIC_proteoform_id",
                               "TOPPIC_database_protein_sequence": "TOPPIC_database_sequence",
                               "TOPPIC_protein_N-terminal_form": "TOPPIC_protein_n-terminal_form", 
                               "TOPPIC_fixed_PTMs": "TOPPIC_fixed_modifications",
                               "TOPPIC_variable_PTMs": "TOPPIC_variable_modifications",
                               "TOPPIC_MIScore": "TOPPIC_miscore",
                               "TOPPIC_number_of_variable_PTMs": "TOPPIC_number_of_variable_modifications",         
                               "TOPPIC_E-value": "TOPPIC_e-value", 
                               "TOPPIC_spectrum-level_Q-value": "TOPPIC_spectrum-level_q-value",
                               "TOPPIC_proteoform-level_Q-value": "TOPPIC_proteoform-level_q-value"}
                            )
    
    # filter based on precursor mass
    query_ms_df = peak_list_extract_based_precursor(query_ms_df)
    # print(f"remove precursor masses: {len(query_ms_df)}")
    # get spectra data and 50 peaks
    query_ms_df, query_ms_org_df = ms_incorr_remove(query_ms_df) 
    # print(f"query_ms_df: {len(query_ms_df)}")
    # convert to vector, log and normalization
    query_sorted_mass_inte_ch_df = ms_spectrum_preprocess(query_ms_df)
    query_sorted_mass_inte_ch_df['spectrum_id'] = query_ms_df['spectrum_id'].values
    
    # query
    tol_val = 10 
    tol_type = 'ppm' 
    frag_tol_val=10
    pre_charge_use='no'
    combined_target_decoy_res = ms_query(target_decoy_ms_rep1, target_decoy_inte_rep_convert, target_decoy_ch_rep, 
                                         query_ms_df, query_sorted_mass_inte_ch_df, tol_type, tol_val, frag_tol_val, pre_charge_use)  
    
    # FDR calculation
    spec_id_1fdr = fdr_cal(combined_target_decoy_res)
    lib_search = combined_target_decoy_res.loc[(combined_target_decoy_res['spectrum_id'].isin(spec_id_1fdr)) & (combined_target_decoy_res['flag']==1),'spectrum_id'].values
    lib_ident1 = combined_target_decoy_res[combined_target_decoy_res['spectrum_id'].isin(lib_search)]
    lib_ident1 = lib_ident1.sort_values(by='spectrum_id')
    lib_ident1 = lib_ident1.drop(['lib_index', 'flag','calculated_cosine_dist'], axis=1) 
    lib_ident1.reset_index(inplace = True, drop=True)
    sim_cutoff = 0.3
    lib_ident1_sim = lib_ident1[lib_ident1['similarity']>=sim_cutoff]                     
    
    print(f"1FDR: {len(lib_ident1)}")
    print(f'Total number of spectra identified by library with similarity>{sim_cutoff}: {len(lib_ident1_sim)}')
    # num_ptm = len(lib_ident1_sim['representative_id'].unique())
    num_clusters = len(lib_ident1_sim['representative_id'].unique())
    # print(f'number of proteoforms identified by library: {num_ptm}')
    print(f'number of identified clusters: {num_clusters}')  
    
    
    # curr_path = r"H:\sw480_datasets"
    # wfile = os.path.join(curr_path, w_filename)
    print("Query results have been stored under toplib_output folder!")
    # write_result_fun(wfile, combined_target_decoy_res)
    write_result_fun(wfile, lib_ident1)
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="top-down MS spectral library query")
    parser.add_argument(
        "-d", "--db",
        required=True,
        help="Input db file(library)"
    )

    parser.add_argument(
        "-m", "--msalign",
        required=True,
        help="Input msalign file(query)"
    )

    parser.add_argument(
        "-t", "--tsv",
        required=True,
        help="Input TSV file(query)"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file for TopLib database search results"
    )

    args = parser.parse_args()
    if not os.path.isfile(args.msalign):
        print(f"Error: msalign file not found: {args.msalign}")
        sys.exit(1)

    if not os.path.isfile(args.tsv):
        print(f"Error: TSV file not found: {args.tsv}")
        sys.exit(1)

    ms_library_search(
        lib_file=args.db,
        tsv_query_file=args.tsv,
        msalign_query=args.msalign,
        wfile=args.output
    )
    
