import sys
import os
import pandas as pd
import numpy as np
import numba as nb
from distance_calculation import dist_func_cal
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as shc 
import msalign_file as mf
from db_add_file import get_msfilename, filename_check
import random
import sqlite3 
from tqdm import tqdm
import time
    

def get_user_input():
    print("Please enter your msalign file name (with extension .msalign):")    
    input_msalignfile = sys.stdin.readline()
    msalign_file = input_msalignfile.strip()
    print("Please enter your tsv file name (with extension .tsv):")    
    input_tsvfile = sys.stdin.readline()
    tsv_file = input_tsvfile.strip()
    print("Please enter your method for generating consensus spectra: 1 = Single; 2 = Average.")    
    input_rep_method = sys.stdin.readline()
    rep_method_code = input_rep_method.strip()    
    if rep_method_code == '1':
        rep_method = 'Single'
    else:
        rep_method = 'Average' 

    user_para_dict = {'file_name': {'msalign_file': msalign_file, 'tsv_file': tsv_file},
                     'representative_method': rep_method}
    return user_para_dict


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
    # add dummy file_id column for this file
    ms_df_new['file_id'] = 1
    return ms_df_new


def ms_rep_library_building(ms_df, msalign_file_name, consensus_method, data_mode):
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
    
    
    def ms_incorr_remove(ms_df):
        # extract 50 peaks and remove incorrect ms
        ms_org_df = ms_df.copy()
        num_pks = 50
        num_spec = len(ms_df)
        ms_index = ms_df.index.values
        drop_idx = []
        for i in range(num_spec):
            spec = ms_df.iloc[i].copy()
            mass_list = np.array(spec['mass'])
            inte_list = np.array(spec['intensity'])
            charge_list = np.array(spec['charge'])
            # sort by intensity and retain top peaks
            sorted_inte_idx = np.argsort(inte_list)[::-1][0:num_pks]
            spec['mass'] = mass_list[sorted_inte_idx]
            spec['intensity'] = inte_list[sorted_inte_idx]
            spec['charge'] = charge_list[sorted_inte_idx]
            ms_df.iloc[i] = spec
            if len(mass_list[sorted_inte_idx])<=1:
                drop_idx.append(ms_index[i])
            
        # drop num mass <=1
        ms_df.drop(drop_idx, inplace = True)
        ms_org_df.drop(drop_idx, inplace = True)
        # drop invalid precursor intensity
        drop_idx2 = ms_df[ms_df['precursor_intensity']==0].index.values
        ms_df.drop(drop_idx2, inplace = True)
        ms_org_df.drop(drop_idx2, inplace = True)
    
        ms_df.reset_index(drop=True, inplace=True)
        ms_org_df.reset_index(drop=True, inplace=True)
        return ms_df, ms_org_df
        
    
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
            idx_ch=g_p1.get_group(charge_list[j]).index.values.tolist()
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
    

    pbar = tqdm(total=6, desc="Processing")
    start_time = time.time()        
    # step1: preprocessing the ms data
    # get spectra data and check
    ms_df, ms_org_df = ms_incorr_remove(ms_df)            
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
    if consensus_method == "Average":
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
    lib_ptm = ms_df.loc[ptm_idx_all,['file_id','file_name','spectrum_id','scan','retention_time','precursor_feature_id', 
                                     'proteoform_id', 'Cluster ID','protein_accession','protein_description', 'first_residue', 'last_residue', 
                                     'protein_sequence', 'proteoform', 'e_value', 'spectrum_level_q_value','proteoform_level_q_value','cluster_id']]
    # remove inconsistent proteoform identifications
    lib_ptm_ext = lib_ptm.dropna(subset='proteoform')
    lib_ptm_ext_cluster_id = lib_ptm_ext['Cluster ID'].unique()
    rows_to_drop = []
    for i in range(len(lib_ptm_ext_cluster_id)):
        cluster_ids = lib_ptm_ext.loc[lib_ptm_ext['Cluster ID']==lib_ptm_ext_cluster_id[i],'cluster_id'].unique()
        temp_df = lib_ptm_ext.loc[lib_ptm_ext['Cluster ID']==lib_ptm_ext_cluster_id[i],:].sort_values(by='e_value')
        if len(cluster_ids)>1:
            rows_to_drop.extend(temp_df.index.values[1:])
            # cluster_idx_all.append(i)
    
    lib_ptm_ext_filtered = lib_ptm_ext.drop(rows_to_drop)
    # update the ptm after removing inconsistent proteoform 
    columns_to_replace = ['proteoform_id', 'Cluster ID','protein_accession','protein_description', 'first_residue', 'last_residue', 
                          'protein_sequence', 'proteoform', 'e_value', 'spectrum_level_q_value','proteoform_level_q_value']
    lib_ptm.loc[~lib_ptm.index.isin(lib_ptm_ext_filtered.index),columns_to_replace] = np.nan
    lib_ptm = lib_ptm.drop(['Cluster ID'], axis=1)
    elapsed_time = time.time() - start_time
    pbar.set_description(f"Processing ({elapsed_time:.1f}s)")
    time.sleep(5/6)
    pbar.update()
    
    # step5: generate decoy mass library            
    # generate decoy library for mass
    # decoy_mass_df = decoy_mass_rep_gen(ms_mass_inte_ch_rep)
    decoy_mass_df = decoy_mass_gen(ms_mass_inte_ch_rep)
    decoy_mass_inte_ch_rep = pd.concat([ms_rep_df, decoy_mass_df, inte_rep_df, charge_rep_df, inte_rep_cvt_df], axis=1)
    
    ms_mass_inte_ch_rep['flag'] = np.ones(len(ms_mass_inte_ch_rep), dtype='int32')
    decoy_mass_inte_ch_rep['flag'] = np.zeros(len(decoy_mass_inte_ch_rep), dtype='int32')
    # combine
    target_decoy_ms_rep = pd.concat([ms_mass_inte_ch_rep, decoy_mass_inte_ch_rep], axis=0)
    target_decoy_ms_rep.reset_index(inplace = True, drop=True)

    # add addition information for proteoform, protein accession, seuqnce, residues etc information
    target_decoy_ptm_rep = pd.concat([lib_ptm, lib_ptm], axis=0)
    target_decoy_ptm_rep.reset_index(inplace = True, drop=True)
    
    target_decoy_ms_rep = pd.concat([target_decoy_ms_rep, target_decoy_ptm_rep], axis=1)
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
        # print("\nDirectory '%s' created successfully" % directory)
        lib_name = get_msfilename(msalign_file_name) # remove the path
        if consensus_method == 'Single':
            lib_filename = lib_name + '_single_representatives.db'
        else:
            lib_filename = lib_name + '_average_representatives.db'
        if data_mode=='file':
            # output representatives to a db file
            wfile = os.path.join(curr_path, directory, lib_filename)   
            conn = sqlite3.connect(wfile)
            target_decoy_ms_rep_ID.to_sql(name = 'target_decoy_spectra_representatives', con = conn, index=False, if_exists='replace') 
            conn.close()
        # extract target representatives
        target_spectra_reps = target_decoy_ms_rep[target_decoy_ms_rep['flag']==1]
        return ms_df, target_spectra_reps

    except OSError as error:
        print("\nDirectory '%s' can not be created" % directory)
        return None
            

if __name__ == "__main__":
    if len(sys.argv) == 1:
        user_inputs = get_user_input()
        msalign_file_name = user_inputs['file_name']['msalign_file']
        tsv_file_name = user_inputs['file_name']['tsv_file']
        consensus_method = user_inputs['representative_method']
        # file name check
        curr_path = os.getcwd()
        directory = "TopLib"
        msalign_file_name = os.path.join(curr_path, directory, msalign_file_name)  
        tsv_file_name = os.path.join(curr_path, directory, tsv_file_name)             
        filenames = [msalign_file_name, tsv_file_name]
        file_flag = filename_check(filenames)
        if file_flag == 1:
            # get spectra data
            ms_df = get_spectra_by_file(msalign_file_name, tsv_file_name)
            # representative building
            spectra_df, ms_spectra_rep = ms_rep_library_building(ms_df, msalign_file_name, consensus_method, 'file')
        else:
            print('At least one of input files is incorrect and please check file name and path!')
    else:
        print("Usage: python script.py")
        sys.exit()
        
