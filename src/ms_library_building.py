import sys
import os
import argparse
import pandas as pd
import numpy as np
import numba as nb
from distance_calculation import dist_func_cal
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as shc 
import msalign_file as mf
import random
import sqlite3 
from tqdm import tqdm
import time
    

def user_command():
    parser = argparse.ArgumentParser(description="Top-down library parameters specification")
    parser.add_argument("msalign_filename", type=str,
                        help="An msalign file name (*.msalign)")
    parser.add_argument("tsv_filename", type=str,
                        help="A tsv file name (*.tsv)")
    parser.add_argument("rep_method", type=str,
                        help="Representative method (average or single)")
    parser.add_argument("-c", "--PPIR_cutoff", type=float, default = 0, 
                        help="PPIR cutoff. Default value: 0")
    args = parser.parse_args()
    msalign_filename = args.msalign_filename
    tsv_filename = args.tsv_filename
    rep_method = args.rep_method.lower()
    cutoff = args.PPIR_cutoff
    return msalign_filename, tsv_filename, rep_method, cutoff 


def get_spectra_by_file(msalign_file_name, tsv_file_name):
    df_all = pd.read_csv(tsv_file_name, sep='\t')
    ms_df, mass_df = mf.read_msalign(msalign_file_name)
    # add protemform info
    idx_ext_all = []
    tsv_idx_all = []
    for i in range(len(df_all)):
        spec_id = df_all['Spectrum ID'].iloc[i]
        idx_ext = ms_df[ms_df['spectrum_id']==spec_id].index.values
        if len(idx_ext):
            idx_ext_all.extend(idx_ext)
            tsv_idx_all.append(i)
        
    ms_df_new = ms_df.copy()
    ms_df_new.loc[idx_ext_all,'proteoform_id'] = df_all['Proteoform ID'].iloc[tsv_idx_all].values
    ms_df_new.loc[idx_ext_all,'protein_accession'] = df_all['Protein accession'].iloc[tsv_idx_all].values
    ms_df_new.loc[idx_ext_all,'protein_description'] = df_all['Protein description'].iloc[tsv_idx_all].values
    ms_df_new.loc[idx_ext_all,'first_residue'] = df_all['First residue'].iloc[tsv_idx_all].values
    ms_df_new.loc[idx_ext_all,'last_residue'] = df_all['Last residue'].iloc[tsv_idx_all].values
    ms_df_new.loc[idx_ext_all,'protein_sequence'] = df_all['Database protein sequence'].iloc[tsv_idx_all].values 
    ms_df_new.loc[idx_ext_all,'proteoform'] = df_all['Proteoform'].iloc[tsv_idx_all].values
    ms_df_new.loc[idx_ext_all,'e_value'] = df_all['E-value'].iloc[tsv_idx_all].values
    ms_df_new.loc[idx_ext_all,'spectrum_level_q_value'] = df_all['Spectrum-level Q-value'].iloc[tsv_idx_all].values
    ms_df_new.loc[idx_ext_all,'proteoform_level_q_value'] = df_all['Proteoform-level Q-value'].iloc[tsv_idx_all].values
    
    column_map = {
        'num_unexpected_modifications': '#unexpected modifications',
        'unexpected_modifications': 'unexpected modifications',
        'fixed_ptms': 'Fixed PTMs',
        'adjusted_precursor_mass': 'Adjusted precursor mass',        
        'prsm_id': 'Prsm ID',
        'fragmentation': 'Fragmentation',        
        'num_original_peaks': '#peaks',
        'proteoform_intensity': 'Proteoform intensity', 
        'feature_intensity': 'Feature intensity',
        'feature_score': 'Feature score',
        'feature_apex_time':'Feature apex time',
        'num_protein_hits':'#Protein hits',
        'proteoform_mass': 'Proteoform mass',
        'protein_n_terminal_form': 'Protein N-terminal form',
        'num_variable_ptms': '#variable PTMs',
        'variable_ptms': 'variable PTMs',
        'num_matched_peaks': '#matched peaks',
        'num_matched_fragment_ions': '#matched fragment ions',
        'special_amion_acids': 'Special amino acids',
        'miscore': 'MIScore'
        }

    for col_in_ms_df, col_in_df in column_map.items():
        if col_in_df in df_all.columns:
            ms_df_new.loc[idx_ext_all, col_in_ms_df] = df_all[col_in_df].iloc[tsv_idx_all].values
    # add dummy file_id column for this file
    ms_df_new['file_id'] = 1
    return ms_df_new


def ms_first_intensity_remove(ms_df, cutoff):
    ms_idx = ms_df.index.values
    idx_sel = []
    for i in range(len(ms_df)):
        pre_inte = ms_df['precursor_intensity_list'].iloc[i]
        if pre_inte[0]:
            first_inte = pre_inte[0] / sum(pre_inte)
            if first_inte >= cutoff:
                idx_sel.append(ms_idx[i])

    ms_df_sel = ms_df.loc[idx_sel,:]
    ms_df_sel.reset_index(drop=True, inplace=True)
    return ms_df_sel
    

def ms_rep_library_building(lib_ms_df, rep_method, data_mode):
    def ms_spectrum_preprocess(ms_all,alg):
        # sort, log2
        sorted_mass_all = []
        sorted_int_all = []
        sorted_ch_all = []
        for i in range(len(ms_all)):
            # sort
            sorted_mass = np.sort(ms_all['mass'].iloc[i])
            sorted_mass_idx = np.argsort(ms_all['mass'].iloc[i])
            sorted_int = np.array(ms_all['intensity'].iloc[i])[sorted_mass_idx]
            if alg =='log':
                sorted_int = np.log2(sorted_int)
                sorted_int = sorted_int / np.linalg.norm(sorted_int)
            sorted_ch = np.array(ms_all['charge'].iloc[i])[sorted_mass_idx]
            sorted_int_all.append(sorted_int)
            sorted_mass_all.append(sorted_mass)
            sorted_ch_all.append(sorted_ch)
    
        data = {'mass': sorted_mass_all,
                'intensity': sorted_int_all,
                'charge': sorted_ch_all}
        sorted_mass_int_ch_df = pd.DataFrame(data)
        return sorted_mass_int_ch_df
    
    
    def ms_incorr_remove(ms_df, num_pks=50):
        ms_org_df = ms_df.copy()
        mass_arr= ms_df['mass'].to_numpy()
        intensity_arr = ms_df['intensity'].to_numpy()
        charge_arr = ms_df['charge'].to_numpy()
        precursor_intensity = ms_df['precursor_intensity'].values
        ms_index = ms_df.index.values
        num_spec = len(ms_df)
        valid_indices = []
        for i in range(num_spec):
            mass_list = np.array(mass_arr[i])
            inte_list = np.array(intensity_arr[i])
            charge_list = np.array(charge_arr[i])
            # Sort by intensity, keep top peaks
            sorted_inte_idx = np.argsort(inte_list)[::-1][0:num_pks]
            mass_arr[i] = mass_list[sorted_inte_idx]
            intensity_arr[i] = inte_list[sorted_inte_idx]
            charge_arr[i] = charge_list[sorted_inte_idx]
            
            # Filter out spectra with <=1 peak and zero precursor intensity
            if len(mass_list[sorted_inte_idx]) > 1 and precursor_intensity[i] != 0:
                valid_indices.append(i)
    
        valid_ms_index = ms_index[valid_indices]
        ms_ext_df = ms_df.loc[valid_ms_index,:]
        ms_org_ext_df = ms_org_df.loc[valid_ms_index,:]
        ms_ext_df['mass'] = mass_arr[valid_indices]
        ms_ext_df['intensity'] = intensity_arr[valid_indices]
        ms_ext_df['charge'] = charge_arr[valid_indices]
        
        ms_ext_df.reset_index(drop=True, inplace=True)
        ms_org_ext_df.reset_index(drop=True, inplace=True)
        return ms_ext_df, ms_org_ext_df
        
    
    def clustering_index_50peaks(mz_all, int_all, ch_all, err, ppm, thre, idx1, metric):
        """
        hierarchical clustering: 
        thre: cutoff threshold (distance),
        idx1: id list after precursor mass and charge  
        return: id list after clustering 
        """
        cluster_index_all = []
        for i in range(len(idx1)):
            if len(idx1[i])>1:
                Y = dist_func_cal(mz_all[idx1[i]].values, int_all[idx1[i]].values, ch_all[idx1[i]].values, err, ppm, metric) # pairwise distances                
                Z = shc.linkage(Y, method='average')
                labels = shc.fcluster(Z, t=thre, criterion='distance')-1
                c_index = cluster_index_gen(labels)
                cluster_index_all.extend([sorted(list(np.array(idx1[i])[x])) for x in c_index])
            else:
                cluster_index_all.extend([idx1[i]])
        return cluster_index_all


    def cluster_index_gen(data_cluster):
        """
        find out the indices after clustering
        """
        cluster_index_all = []
        for i in range(len(np.unique(data_cluster))):
            cluster_index = [k for k,v in enumerate(data_cluster) if v==i]
            cluster_index_all.append(cluster_index)
        # Remove empty List from List
        while [] in cluster_index_all:
            cluster_index_all.remove([])
        # if -1 exists--> single class
        cluster_single_index = [k for k,v in enumerate(data_cluster) if v==-1]
        if cluster_single_index:
            for i in range(len(cluster_single_index)):
                cluster_index_all.append([cluster_single_index[i]])
        return cluster_index_all


    def precursor_charge_matched_index(precursor_all):
        # group-by charge
        ch_idx = []
        g_p1=precursor_all.groupby(['precursor_charge'])
        charge_list = np.unique(precursor_all['precursor_charge'].values)
        for j in range(len(charge_list)):
            # idx_ch=g_p1.get_group(charge_list[j]).index.values.tolist()
            idx_ch = g_p1.get_group((charge_list[j],)).index.values.tolist()
            ch_idx.append(idx_ch)
        return ch_idx


    def clustering_precursor_mass_index(precursor_all, precursor_tol):
        cluster_index = precursor_charge_matched_index(precursor_all)
        cluster_index_all = []
        for i in range(len(cluster_index)):
            if len(cluster_index[i])>1:
                mz = precursor_all.loc[cluster_index[i],'precursor_mass'].values
                c_mz = mz[:, None] # column vectors
                Y = pdist(c_mz, 'cityblock')
                Z = shc.linkage(Y, method='complete')
                labels = shc.fcluster(Z, t=precursor_tol, criterion='distance')-1
                c_index = cluster_index_gen(labels)
                cluster_index_all.extend([sorted(list(np.array(cluster_index[i])[x])) for x in c_index])
            else:
                cluster_index_all.extend([cluster_index[i]])
        return cluster_index_all

    
    def cluster_representative_50pks(cluster_index, sorted_mass_int_ch_df, precursor_df):
        mass_rep_cluster_all = []
        inte_rep_cluster_all = []
        ch_rep_cluster_all = []
        precursor_mass_rep_all = []
        precursor_inte_rep_all = []
        precursor_charge_all = []
        for i in range(len(cluster_index)):
            if len(cluster_index[i])>1:
                c_index = cluster_index[i]
                mass_int_ch_cluster = sorted_mass_int_ch_df.iloc[cluster_index[i]]
                mass_rep_cluster, inte_rep_cluster, ch_rep_cluster = spectra_representatives_cluster(c_index, mass_int_ch_cluster)
                # calculate precursor mass representatives 
                preMass_c = precursor_df['precursor_mass'].iloc[cluster_index[i]].values
                preInte_c = precursor_df['precursor_intensity'].iloc[cluster_index[i]].values
                precursor_mass_rep, precursor_inte_rep = precursor_mass_cluster_rep(preMass_c, preInte_c)
                precursor_ch = precursor_df['precursor_charge'].iloc[cluster_index[i][0]]
            else:
                mass_rep_cluster = sorted_mass_int_ch_df['mass'].iloc[cluster_index[i][0]]
                inte_rep_cluster = sorted_mass_int_ch_df['intensity'].iloc[cluster_index[i][0]]
                ch_rep_cluster = sorted_mass_int_ch_df['charge'].iloc[cluster_index[i][0]]
                precursor_mass_rep = precursor_df['precursor_mass'].iloc[cluster_index[i][0]]
                precursor_inte_rep = precursor_df['precursor_intensity'].iloc[cluster_index[i][0]]
                precursor_ch = precursor_df['precursor_charge'].iloc[cluster_index[i][0]]
    
            # extract 50 largest peaks for ms/ms
            sorted_idx = np.argsort(inte_rep_cluster)
            inte_rep_cluster = inte_rep_cluster[sorted_idx[-50:]]
            mass_rep_cluster = mass_rep_cluster[sorted_idx[-50:]]
            ch_rep_cluster = ch_rep_cluster[sorted_idx[-50:]]
            # sort by mass value
            sorted_mass_idx = np.argsort(mass_rep_cluster)
            mass_rep_cluster = np.sort(mass_rep_cluster)
            inte_rep_cluster = inte_rep_cluster[sorted_mass_idx]
            ch_rep_cluster = ch_rep_cluster[sorted_mass_idx]
            # ms/ms
            mass_rep_cluster_all.append(mass_rep_cluster)
            inte_rep_cluster_all.append(inte_rep_cluster)
            ch_rep_cluster_all.append(ch_rep_cluster)
            # precursor      
            precursor_mass_rep_all.append(precursor_mass_rep)
            precursor_inte_rep_all.append(precursor_inte_rep)
            precursor_charge_all.append(precursor_ch)
    
        return mass_rep_cluster_all, inte_rep_cluster_all, ch_rep_cluster_all, precursor_mass_rep_all, precursor_inte_rep_all, precursor_charge_all
    
    
    def precursor_mass_cluster_rep(preMass_c, preInte_c):
        # generate precursor mass representatives in a cluster
        """
        preMass_c: precursor mass values in a cluster
        preInte_c: precursor intesity values in a cluster
        """
        preMass_diff = np.intp(preMass_c) - np.intp(preMass_c[0])
        unique, counts = np.unique(preMass_diff, return_counts = True)
        index = np.argmax(counts)
        idx = [ix for ix, x in enumerate(preMass_diff) if x==unique[index]]
        # weighted average
        preMass_bin_weighted_avg = ((preMass_c[idx]) * (preInte_c[idx])).sum() / (preInte_c[idx]).sum()
        preInte_bin_avg = (preInte_c[idx]).mean()
        return preMass_bin_weighted_avg, preInte_bin_avg
    
    
    def spectra_representatives_cluster(c_index, sorted_mass_int_ch_df):
        #  generate a representative for a cluster after clustering
        sorted_mass1 = sorted_mass_int_ch_df['mass'].iloc[0]
        sorted_mass2 = sorted_mass_int_ch_df['mass'].iloc[1]
        sorted_int1 = sorted_mass_int_ch_df['intensity'].iloc[0]
        sorted_int2 = sorted_mass_int_ch_df['intensity'].iloc[1]
        sorted_ch1 = sorted_mass_int_ch_df['charge'].iloc[0]
        sorted_ch2 = sorted_mass_int_ch_df['charge'].iloc[1]
            
        err = 0.0
        ppm = 10
        tol2 = np.maximum(err*np.ones(np.shape(sorted_mass1)), ppm*sorted_mass1*1e-6)
        th1_all = 0.0+tol2
        w_mass_all, int_sum_all, ch_all = weighted_mass_sum_inte_cal(sorted_mass1, sorted_int1, sorted_ch1, 
                                                                     sorted_mass2, sorted_int2, sorted_ch2, 
                                                                     th1_all)
    
        w_mass_idx = np.argsort(w_mass_all)
        sorted_w_mass_all = np.sort(w_mass_all)
        sorted_int_sum_all = np.array(int_sum_all)[w_mass_idx]
        sorted_ch_all = np.array(ch_all)[w_mass_idx]    
    
        for j in range(2, len(c_index)):
            tol2 = np.maximum(err*np.ones(np.shape(sorted_w_mass_all)), ppm*sorted_w_mass_all*1e-6)
            th1_all = 0.0+tol2
            sorted_mass2 = sorted_mass_int_ch_df['mass'].iloc[j]
            sorted_int2 = sorted_mass_int_ch_df['intensity'].iloc[j]
            sorted_ch2 = sorted_mass_int_ch_df['charge'].iloc[j]
            w_mass1, int1_sum, ch1 = weighted_mass_sum_inte_cal(sorted_w_mass_all, sorted_int_sum_all, sorted_ch_all, 
                                                                sorted_mass2, sorted_int2, sorted_ch2, 
                                                                th1_all)
    
    
            w_mass_idx = np.argsort(w_mass1)
            sorted_w_mass_all = np.sort(w_mass1)
            sorted_int_sum_all = np.array(int1_sum)[w_mass_idx]
            sorted_ch_all = np.array(ch1)[w_mass_idx]   
    
        return sorted_w_mass_all, sorted_int_sum_all, sorted_ch_all
    
    
    @nb.njit
    def weighted_mass_sum_inte_cal(sorted_mass1, sorted_int1, sorted_ch1, sorted_mass2, sorted_int2, sorted_ch2, th1_all):
        w_mass_all = []
        int_sum_all = []
        ch_all = []
        idx1 = []
        idx2 = []
        for i in range(len(sorted_mass1)):
            int_sum = 0
            w_mass_sum = 0
            ch_matched = []
            th1 = th1_all[i]
            for j in range(len(sorted_mass2)):
                if sorted_ch1[i] == sorted_ch2[j]:
                    if abs(sorted_mass2[j] - sorted_mass1[i]) <= th1:
                        if i not in idx1:
                            w_mass_sum = w_mass_sum + (sorted_mass1[i] * sorted_int1[i] + sorted_mass2[j] * sorted_int2[j]) 
                            int_sum = int_sum + (sorted_int1[i] + sorted_int2[j])
                            ch_matched.append(sorted_ch1[i])
                        else:
                            w_mass_sum = w_mass_sum + (sorted_mass2[j] * sorted_int2[j]) 
                            int_sum = int_sum + sorted_int2[j]  
    
                        idx1.append(i)
                        idx2.append(j)
    
            if int_sum !=0:
                w_mass_all.append(w_mass_sum / int_sum)
                int_sum_all.append(int_sum)
                ch_all.extend(ch_matched)
    
        for i in range(len(sorted_mass1)):
            if i not in idx1:
                w_mass_all.append(sorted_mass1[i])
                int_sum_all.append(sorted_int1[i])
                ch_all.append(sorted_ch1[i])
    
        for j in range(len(sorted_mass2)):
            if j not in idx2:
                w_mass_all.append(sorted_mass2[j])
                int_sum_all.append(sorted_int2[j])
                ch_all.append(sorted_ch2[j])
    
        return w_mass_all, int_sum_all, ch_all
    
    
    def cluster_representative_min_evalue(ms_df):
        lib_cluster_ids = np.unique(ms_df['cluster_id'].values)
        idx_all = []
        for i in range(len(lib_cluster_ids)):
            temp_df = ms_df[ms_df['cluster_id']==lib_cluster_ids[i]]
            if temp_df['e_value'].notnull().any():
                idx = temp_df['e_value'].idxmin() 
            else:
                idx = temp_df['num_mass'].idxmax() 
            idx_all.append(idx)
        lib_mass_inte_ch_ext = ms_df.loc[idx_all,['precursor_mass','precursor_intensity','precursor_charge',
                                                  'mass','intensity','charge']]
        # sort by mass value
        sorted_mass_v_all = []
        sorted_inte_v_all = []
        sorted_ch_v_all = []
        for i in range(len(lib_mass_inte_ch_ext)):
            mass_v = lib_mass_inte_ch_ext['mass'].iloc[i]
            inte_v = lib_mass_inte_ch_ext['intensity'].iloc[i]
            ch_v = lib_mass_inte_ch_ext['charge'].iloc[i]
            sorted_mass_idx = np.argsort(mass_v)
            sorted_mass_v = np.sort(mass_v)
            sorted_inte_v = inte_v[sorted_mass_idx]
            sorted_ch_v = ch_v[sorted_mass_idx]
            sorted_mass_v_all.append(sorted_mass_v)
            sorted_inte_v_all.append(sorted_inte_v)
            sorted_ch_v_all.append(sorted_ch_v)
        pre_mass = lib_mass_inte_ch_ext['precursor_mass'].values.tolist()
        pre_inte = lib_mass_inte_ch_ext['precursor_intensity'].values.tolist()
        pre_ch = lib_mass_inte_ch_ext['precursor_charge'].values.tolist()
        return sorted_mass_v_all, sorted_inte_v_all, sorted_ch_v_all, pre_mass, pre_inte, pre_ch

    
    def decoy_mass_gen(spectra_mass_rep):
        # amino acids information
        aa_name = ['alanine', 'arginine', 'asparagine', 'aspartic acid', 'cysteine', 'glutamic acid', 'glutamine', 'glycine', 'histidine', 'isoleucine',
                    'lysine', 'methionine', 'phenylalanine', 'proline', 'serine', 'threonine', 'tryptophan', 'tyrosine', 'valine']
        aa_mono_mass = [71.03711381, 156.10111105, 114.04292747, 115.02694307, 103.00918451, 129.04259314, 128.05857754,  57.02146373,
                   137.05891188, 113.08406401, 128.09496305, 131.04048464, 147.06841395,  97.05276387,  87.03202844,
                   101.0476785 , 186.07931298, 163.06332857,  99.06841394]
        data = {'amino acid': aa_name,
                'calculated mono mass': aa_mono_mass}
        aa_mass_df = pd.DataFrame(data)
        aa_mass_val = aa_mass_df['calculated mono mass'].values # use mono mass
        columns1 = list(map(lambda x: 'mass_' + str(x), np.arange(50)))
        decoy_mass_all = []
        for j in range(len(spectra_mass_rep)):
            random.seed(j)
            mass1 = spectra_mass_rep[columns1].iloc[j].values
            mass1 = mass1[~np.isnan(mass1)]
            pre_mass1 = spectra_mass_rep['precursor_mass'].iloc[j]
            decoy_mass = np.zeros((len(mass1),), dtype=float)
            for i in range(len(mass1)):
                aa_random_select = random.sample(aa_mass_df.index.tolist(), 2)
                aa_mass_diff = aa_mass_val[aa_random_select][0]-aa_mass_val[aa_random_select][1]
                mass1_shift = mass1[i] + aa_mass_diff
                while mass1_shift<0:
                    aa_random_select = random.sample(aa_mass_df.index.tolist(), 2)
                    aa_mass_diff = aa_mass_val[aa_random_select][0]-aa_mass_val[aa_random_select][1]
                    mass1_shift = mass1[i] + aa_mass_diff
                while mass1_shift > pre_mass1:
                    aa_random_select = random.sample(aa_mass_df.index.tolist(), 2)
                    aa_mass_diff = aa_mass_val[aa_random_select][0]-aa_mass_val[aa_random_select][1]
                    mass1_shift = mass1[i] + aa_mass_diff
                decoy_mass[i] = mass1_shift
            decoy_mass_all.append(decoy_mass)
    
        decoy_mass_df = pd.DataFrame(data=decoy_mass_all, columns=columns1)
        return decoy_mass_df
    

    def representative_table_gen(target_decoy_spectra_reps, rep_name, conn):
        # generate representative table in library
        col_drop_names = ['file_name','scan', 'title', 'retention_time', 'ms_level', 'ms_one_id', 'ms_one_scan', 'precursor_window_begin', 'precursor_window_end', 'activation','precursor_mz',
                          'precursor_feature_id', 'proteoform_id', 'protein_accession', 'protein_description', 'first_residue', 'last_residue',
                          'protein_sequence', 'proteoform', 'e_value', 'spectrum_level_q_value', 'proteoform_level_q_value', 'num_unexpected_modifications', 'unexpected_modifications', 'fixed_ptms',
                          'adjusted_precursor_mass', 'prsm_id', 'fragmentation', 'num_original_peaks','proteoform_intensity', 'feature_intensity', 'feature_score', 'feature_apex_time', 'num_protein_hits',
                          'proteoform_mass','protein_n_terminal_form', 'num_variable_ptms', 'variable_ptms', 'num_matched_peaks', 'num_matched_fragment_ions', 'special_amion_acids', 'miscore', 'cluster_id','flag']    
        drop_cols = [col for col in col_drop_names if col in target_decoy_spectra_reps.columns]
        target_decoy_spectra_reps = target_decoy_spectra_reps.drop(drop_cols, axis=1)
        target_decoy_spectra_reps.to_sql(name = rep_name + '_representatives', con = conn, index=False, if_exists='append') 
        print(rep_name + '_representatives generated for this file!\n')
    

    pbar = tqdm(total=6, desc="Processing")
    start_time = time.time()        
    # step1: preprocessing the ms data
    # get spectra data and check
    ms_df, ms_org_df = ms_incorr_remove(lib_ms_df) 
    # convert to vector, log and normalization
    sorted_mass_int_ch_df = ms_spectrum_preprocess(ms_df,'log')
    elapsed_time = time.time() - start_time
    pbar.set_description(f"Processing: ({elapsed_time:.1f}s)")
    time.sleep(1/6)
    pbar.update()
    
    # step2:clustering
    metric = 'cosine'
    precursor_tol = 2.2
    ppm = 10
    err = 0.0
    thre = 0.7
    cluster_index2 = clustering_precursor_mass_index(ms_df, precursor_tol)
    cluster_index = clustering_index_50peaks(sorted_mass_int_ch_df['mass'], sorted_mass_int_ch_df['intensity'], sorted_mass_int_ch_df['charge'], err, ppm, thre, cluster_index2, metric)
    # find each spectrum's cluster id
    cluster_id = np.zeros((len(ms_df),), dtype=int)
    for i in range(len(cluster_index)):
        temp_idx = cluster_index[i]
        cluster_id[temp_idx] = np.ones((len(temp_idx),), dtype=int) * i
    
    ms_df['cluster_id']=cluster_id    
    elapsed_time = time.time() - start_time
    pbar.set_description(f"Processing ({elapsed_time:.1f}s)")
    time.sleep(2/6)
    pbar.update()
    
    # step3:generate representatives
    # read ms_org_df with all peaks
    mass_v_df = ms_org_df.loc[:,['mass','intensity','charge']]
    # sort mass
    mass_v_df2 =  ms_spectrum_preprocess(mass_v_df,'sort')
    # generate ms representatives
    if rep_method == "average":
        mass_rep, inte_rep, ch_rep, pre_mass_rep, pre_int_rep, pre_ch = cluster_representative_50pks(cluster_index, mass_v_df2, ms_df)
    else:
        mass_rep, inte_rep, ch_rep, pre_mass_rep, pre_int_rep, pre_ch = cluster_representative_min_evalue(ms_df)
    
    elapsed_time = time.time() - start_time
    pbar.set_description(f"Processing ({elapsed_time:.1f}s)")
    time.sleep(3/6)
    pbar.update()
    # convert intensity to log scale            
    inte_rep_convert = []
    for i in range(len(inte_rep)):
        data_inte = inte_rep[i]
        data_inte = np.log2(data_inte)
        data_inte = data_inte / np.linalg.norm(data_inte)
        inte_rep_convert.append(data_inte)
    
    elapsed_time = time.time() - start_time
    pbar.set_description(f"Processing ({elapsed_time:.1f}s)")
    time.sleep(4/6)
    pbar.update()
    
    #step4: building the ms representative library
    data = {'representative_id': np.arange(len(pre_mass_rep)) + 1,
            'precursor_mass': pre_mass_rep,
            'precursor_intensity': pre_int_rep,
            'precursor_charge': pre_ch}
    ms_rep_df = pd.DataFrame(data)
    # generate a dataframe
    columns1 = list(map(lambda x: 'mass_' + str(x), np.arange(50)))
    columns2 = list(map(lambda x: 'intensity_' + str(x), np.arange(50)))
    columns3 = list(map(lambda x: 'charge_' + str(x), np.arange(50)))  
    columns4 = list(map(lambda x: 'norm_intensity_' + str(x), np.arange(50)))  

    mass_rep_df = pd.DataFrame(data=mass_rep, columns=columns1)
    inte_rep_df = pd.DataFrame(data=inte_rep, columns=columns2)
    charge_rep_df = pd.DataFrame(data=ch_rep, columns=columns3)
    # intensity for normalized log scale
    inte_rep_cvt_df = pd.DataFrame(data=inte_rep_convert, columns=columns4)
    ms_mass_inte_ch_rep = pd.concat([ms_rep_df, mass_rep_df, inte_rep_df, charge_rep_df, inte_rep_cvt_df], axis=1)
    # add proteoform, proteoform id, protein accession, etc   
    # add temporary proteoform grouped ID by precursor charge
    ms_df['Cluster ID'] = ms_df.groupby(['precursor_charge', 'proteoform_id']).ngroup()
    lib_cluster_ids = np.unique(ms_df['cluster_id'].values)
    ptm_idx_all = []
    for i in range(len(lib_cluster_ids)):
        temp_df = ms_df.loc[ms_df['cluster_id']==lib_cluster_ids[i],:].sort_values(by='e_value')
        ptm_idx = temp_df.index.values[0]
        ptm_idx_all.append(ptm_idx)
    
    col_names = ['file_id','file_name','spectrum_id','title','scan','retention_time','ms_level','ms_one_id', 'ms_one_scan','precursor_window_begin', 'precursor_window_end', 'activation', 'precursor_mz','precursor_feature_id', 
                 'proteoform_id', 'Cluster ID','protein_accession','protein_description', 'first_residue', 'last_residue', 
                 'protein_sequence', 'proteoform', 'e_value', 'spectrum_level_q_value','proteoform_level_q_value','num_unexpected_modifications', 'unexpected_modifications', 'fixed_ptms',
                 'adjusted_precursor_mass', 'prsm_id', 'fragmentation', 'num_original_peaks','proteoform_intensity', 'feature_intensity', 'feature_score', 'feature_apex_time', 'num_protein_hits',
                 'proteoform_mass','protein_n_terminal_form', 'num_variable_ptms', 'variable_ptms', 'num_matched_peaks', 'num_matched_fragment_ions', 'special_amion_acids', 'miscore', 'cluster_id']
    existing_cols = [col for col in col_names if col in ms_df.columns]
    lib_ptm = ms_df.loc[ptm_idx_all, existing_cols]                 
    
    # remove inconsistent proteoform identifications
    lib_ptm_ext = lib_ptm.dropna(subset='proteoform')
    lib_ptm_ext_cluster_id = lib_ptm_ext['Cluster ID'].unique()
    rows_to_drop = []
    for i in range(len(lib_ptm_ext_cluster_id)):
        cluster_ids = lib_ptm_ext.loc[lib_ptm_ext['Cluster ID']==lib_ptm_ext_cluster_id[i],'cluster_id'].unique()
        temp_df = lib_ptm_ext.loc[lib_ptm_ext['Cluster ID']==lib_ptm_ext_cluster_id[i],:].sort_values(by='e_value')
        if len(cluster_ids)>1:
            rows_to_drop.extend(temp_df.index.values[1:])
    
    lib_ptm_ext_filtered = lib_ptm_ext.drop(rows_to_drop)
    # update the ptm after removing inconsistent proteoform 
    col_to_replace_names = ['proteoform_id', 'Cluster ID','protein_accession','protein_description', 'first_residue', 'last_residue', 
                            'protein_sequence', 'proteoform', 'e_value', 'spectrum_level_q_value','proteoform_level_q_value','num_unexpected_modifications', 'unexpected_modifications','fixed_ptms',
                            'adjusted_precursor_mass', 'prsm_id', 'fragmentation', 'num_original_peaks','proteoform_intensity', 'feature_intensity', 'feature_score', 'feature_apex_time', 'num_protein_hits',
                            'proteoform_mass','protein_n_terminal_form', 'num_variable_ptms', 'variable_ptms', 'num_matched_peaks', 'num_matched_fragment_ions', 'special_amion_acids', 'miscore', 'cluster_id']
    
    columns_to_replace  = [col for col in col_to_replace_names if col in lib_ptm.columns] 
    lib_ptm.loc[~lib_ptm.index.isin(lib_ptm_ext_filtered.index),columns_to_replace] = np.nan
    lib_ptm = lib_ptm.drop(['Cluster ID'], axis=1)
    elapsed_time = time.time() - start_time
    pbar.set_description(f"Processing ({elapsed_time:.1f}s)")
    time.sleep(5/6)
    pbar.update()
    
    # step5: generate decoy mass library            
    # generate decoy library for mass
    decoy_mass_df = decoy_mass_gen(ms_mass_inte_ch_rep)
    decoy_mass_inte_ch_rep = pd.concat([ms_rep_df, decoy_mass_df, inte_rep_df, charge_rep_df, inte_rep_cvt_df], axis=1)
    
    ms_mass_inte_ch_rep['flag'] = np.ones(len(ms_mass_inte_ch_rep), dtype='int32')
    decoy_mass_inte_ch_rep['flag'] = np.zeros(len(decoy_mass_inte_ch_rep), dtype='int32')
    # combine
    target_decoy_ms_rep = pd.concat([ms_mass_inte_ch_rep, decoy_mass_inte_ch_rep], axis=0)
    target_decoy_ms_rep.reset_index(inplace = True, drop=True)

    # add addition information for proteoform, protein accession, seqence, residues etc information
    target_decoy_ptm_rep = pd.concat([lib_ptm, lib_ptm], axis=0)
    target_decoy_ptm_rep.reset_index(inplace = True, drop=True)
    
    target_decoy_ms_rep = pd.concat([target_decoy_ms_rep, target_decoy_ptm_rep], axis=1)
    # add columns of reviewed, reviewed_proteoform and status
    # target_decoy_ms_rep['reviewed_proteoform'] = ''
    # target_decoy_ms_rep['reviewed'] = 'False'
    # target_decoy_ms_rep['reviewed_time'] = ''
    # extract data with identified proteoforms in the library
    target_decoy_ms_rep_ID = target_decoy_ms_rep.dropna(subset='proteoform')
    target_decoy_ms_rep_ID.reset_index(inplace = True, drop=True)
    
    elapsed_time = time.time() - start_time
    pbar.set_description(f"Processing ({elapsed_time:.1f}s)")
    time.sleep(6/6)
    pbar.update()
    pbar.close() 
    
    # step6: save library
    # building a database on your local directory
    curr_path = os.getcwd()
    directory = "TopLib"
    path = os.path.join(curr_path, directory)
    try:
        os.makedirs(path, exist_ok = True)
        if data_mode=='file':
            lib_filename = 'lib_spectra_ms2.db'
            # output representatives to a db file
            wfile = os.path.join(curr_path, directory, lib_filename)   
            conn = sqlite3.connect(wfile)
            target_decoy_ms_rep_ID.to_sql(name = 'target_decoy_spectra_representatives', con = conn, index=False, if_exists='replace') 
            conn.close()
            print('toplib library has been built!')
        else:
            # generate representatives tables
            lib_file = os.path.join(curr_path, directory, "toplib.db")   
            conn = sqlite3.connect(lib_file) 
            target_spectra_reps = target_decoy_ms_rep[target_decoy_ms_rep['flag']==1]
            representative_table_gen(target_spectra_reps, rep_method, conn)  
            conn.close()
        return ms_df

    except OSError as error:
        print("\nDirectory '%s' can not be created" % directory)
        return None
            

if __name__ == "__main__":
    if len(sys.argv) <= 3:
        print("Usage: python script.py <input_msalign_file> <input_tsv_file> <representative_method>")
        sys.exit()
    else:
        msalign_file_name, tsv_file_name, rep_method, cutoff = user_command()
        # file name check
        filenames = [msalign_file_name, tsv_file_name]
        if all(os.path.isfile(f) for f in filenames):
            # get spectra data
            ms_df = get_spectra_by_file(msalign_file_name, tsv_file_name)
            # extract spectra according to ppir cutoff
            if cutoff:
                ms_df = ms_first_intensity_remove(ms_df, cutoff) 
            ms_df = ms_df.drop(['precursor_intensity_list'], axis=1)
            # representative building
            spectra_df = ms_rep_library_building(ms_df, rep_method, 'file')
        else:
            print('At least one of input files is incorrect and please check file name and path!')

        
