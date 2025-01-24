import sys
import os
import pandas as pd
import numpy as np
import numba as nb
from distance_calculation import peak_50_cosine_mass_ppm
import sqlite3 
from tqdm import tqdm
from db_query import proteoformID_update, rep2msalign, rep2tsv


# combined sw480-rep1, 2 and3
def get_rep_data(lib_filename):
    # get the representatives data    
    conn = sqlite3.connect(lib_filename) 
    query = "SELECT * FROM target_decoy_spectra_representatives"
    target_decoy_spectra_rep = pd.read_sql_query(query, con=conn)
    target_decoy_spectra_rep = target_decoy_spectra_rep.drop(['cluster_id'], axis=1)   
    conn.close()
    return target_decoy_spectra_rep


@nb.njit
def precursorMass_matched(pre_mass_lib, pre_ch_lib, spectra_rep_index, pre_mass_query_all, pre_ch_query_all, spectra_query_index):
    idx1_all =[]
    idx2_all =[]
    for i in range(len(pre_mass_query_all)):
        pre_mass_query = pre_mass_query_all[i]
        pre_ch_query = pre_ch_query_all[i]
        idx1 = []
        idx2 = []
        for j in range(len(pre_mass_lib)):
            if pre_ch_query == pre_ch_lib[j]:
                if abs(pre_mass_lib[j] - pre_mass_query) <= 2.2:
                    idx1.append(spectra_rep_index[j])
                    if spectra_query_index[i] not in idx2:
                        idx2.append(spectra_query_index[i])      

        if idx1:
            idx1_all.append(idx1)
            idx2_all.extend(idx2)
    return idx1_all, idx2_all


def dist_matched(data1_mass, data1_inte, data1_ch, data2_mass, data2_inte, data2_ch, idx1, idx2):
    """
    calculate distances and find matched id
    """
    idx11 = []
    idx22 = []
    err = 0.0
    ppm = 10
    thre = 0.7
    for i in range(len(idx2)):
        cos_all = []
        for j in range(len(idx1[i])):
            cos_v = peak_50_cosine_mass_ppm(data1_mass[idx1[i][j]], data1_inte[idx1[i][j]], data1_ch[idx1[i][j]], 
                                            data2_mass[idx2[i]], data2_inte[idx2[i]], data2_ch[idx2[i]], err, ppm)
            cos_all.append(cos_v)
        temp_dis = np.array(cos_all)
        idx_ser = list(np.where(temp_dis <= thre)[0])
        if len(idx_ser):
            idx11.append(list(np.array(idx1[i])[idx_ser]))
            idx22.append(idx2[i])
    return idx11, idx22


def matched_representatives(spectra_rep1, spectra_rep2):
    # find matched spectra between two files
    columns1 = list(map(lambda x: 'mass_' + str(x), np.arange(50)))
    columns3 = list(map(lambda x: 'charge_' + str(x), np.arange(50))) 
    columns4 = list(map(lambda x: 'norm_intensity_' + str(x), np.arange(50)))  
    target_spectra_rep1 = spectra_rep1[spectra_rep1['flag']==1]
    target_spectra_rep2 = spectra_rep2[spectra_rep2['flag']==1]
    decoy_spectra_rep1 = spectra_rep1[spectra_rep1['flag']==0]
    decoy_spectra_rep2 = spectra_rep2[spectra_rep2['flag']==0]
        
    pre_mass1 = target_spectra_rep1['precursor_mass'].values
    pre_ch1 = target_spectra_rep1['precursor_charge'].values
    spectra_rep1_index = target_spectra_rep1.index.values
    
    pre_mass2 = target_spectra_rep2['precursor_mass'].values
    pre_ch2 = target_spectra_rep2['precursor_charge'].values
    spectra_rep2_index = target_spectra_rep2.index.values

    idx1, idx2 = precursorMass_matched(pre_mass1, pre_ch1, spectra_rep1_index, 
                                       pre_mass2, pre_ch2, spectra_rep2_index)
    idx11, idx22 = dist_matched(target_spectra_rep1[columns1].to_numpy(), target_spectra_rep1[columns4].to_numpy(), target_spectra_rep1[columns3].to_numpy(), 
                                target_spectra_rep2[columns1].to_numpy(), target_spectra_rep2[columns4].to_numpy(), target_spectra_rep2[columns3].to_numpy(), idx1, idx2)
    # find unmatched represenative id
    unmatched_idx = [x for x in spectra_rep2_index if x not in idx22]    
    # add unmatched spectra
    spectra_add = target_spectra_rep2.loc[unmatched_idx,:]
    # decoy
    decoy_add = decoy_spectra_rep2[(decoy_spectra_rep2['representative_id'].isin(spectra_add['representative_id'].values))].copy()
    # update representative ids
    spectra_add.reset_index(drop=True, inplace=True)
    decoy_add.reset_index(drop=True, inplace=True)
    
    spectra_add['representative_id'] = spectra_add.index.values + len(target_spectra_rep1) + 1 
    decoy_add['representative_id'] = spectra_add['representative_id'].values    
    # combine
    target_decoy_spectra_combined = pd.concat([target_spectra_rep1, spectra_add, decoy_spectra_rep1, decoy_add], axis=0)
    target_decoy_spectra_combined.reset_index(drop=True, inplace=True)
    # update proteoform id
    target_spectra_combined = target_decoy_spectra_combined[target_decoy_spectra_combined['flag']==1].copy()
    target_spectra_combined_new = proteoformID_update(target_spectra_combined)
    target_decoy_spectra_combined.loc[target_decoy_spectra_combined['flag']==1,'proteoform_id'] = target_spectra_combined_new['proteoform_id'].values
    target_decoy_spectra_combined.loc[target_decoy_spectra_combined['flag']==0,'proteoform_id'] = target_spectra_combined_new['proteoform_id'].values

    return target_decoy_spectra_combined
    

def store_rep_data(lib_filename, rep_df_ID):
    # output representatives to a db file
    curr_path = os.getcwd()
    directory = "TopLib"
    wfile = os.path.join(curr_path, directory, lib_filename)   
    conn = sqlite3.connect(wfile)
    rep_df_ID.to_sql(name = 'target_decoy_spectra_representatives', con = conn, index=False, if_exists='replace') 
    conn.close()
    

if __name__ == "__main__":
    if len(sys.argv) <3:
        print("At leaset two library files, usage: python script.py <db1_filename> <db2_filename> <db3_filename>...")
        sys.exit()
    else: 
        # file name check
        curr_path = os.getcwd()
        directory = "TopLib"
        input_filenames = sys.argv[1:]
        filenames = []
        for i in range(len(input_filenames)):
            toplib_filepath = os.path.join(curr_path, directory, input_filenames[i])
            filenames.append(toplib_filepath)
            
        if all(os.path.isfile(f) for f in filenames):
            rep_df = get_rep_data(filenames[0])
            for f in tqdm(filenames[1:], desc='Combining spectra files'):
                add_rep_df = get_rep_data(f)
                print(len(rep_df))
                rep_df = matched_representatives(rep_df, add_rep_df)
            print('Processing completed!')
            print(f'total spectra representative: {len(rep_df)}')   
            rep_df_ID=rep_df.dropna(subset='proteoform')
            print(f'total spectra representative (ID): {len(rep_df_ID)}') 
            # print(rep_df2[['representative_id','precursor_mass','precursor_charge','flag','proteoform_id']])
            # save as db file
            # lib_filename = 'sw480_rep_combined_ms_average_pre_mass_55.db'
            lib_filename = 'sw480_rep_combined_ms2_single_representatives.db'
            # print(rep_df_ID)
            store_rep_data(lib_filename, rep_df_ID)
            # write to msalign file and tsv file
            target_rep_df = rep_df[rep_df['flag']==1]
            print(len(target_rep_df))
            target_rep_df_ID = rep_df_ID[rep_df_ID['flag']==1]
            print(len(target_rep_df_ID))
            # print(target_rep_df.dropna(subset='proteoform'))
            msalign_wfile = r"C:\Users\kunza\Documents\GitHub\TopLib_web_amazon\static\sw480_rep_combined_single_representatives.msalign"
            rep2msalign(target_rep_df_ID, msalign_wfile)  
            # write tsv file
            tsv_wfile = r"C:\Users\kunza\Documents\GitHub\TopLib_web_amazon\static\sw480_rep_combined_single_representatives.tsv"
            rep2tsv(target_rep_df_ID, tsv_wfile)
        else:
            print('At least one of input files is incorrect and please check file name and path!')   
            
            
