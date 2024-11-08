import sys
import os
import argparse
import pandas as pd
import numpy as np
import numba as nb
import msalign_file as mf
from distance_calculation import peak_50_cosine_mass_ppm
import sqlite3 
from datetime import datetime
from tqdm import tqdm



def user_command():
    parser = argparse.ArgumentParser(description="Top-down mass spectral parameters specification")
    parser.add_argument("library_filename", type=str,
                         help="A library file name (*.db)")
    parser.add_argument("msalign_filename", type=str,
                        help="An msalign file name (*.msalign)")
    parser.add_argument("-T", "--precursor_error_type", type=str, default='ppm', 
                        help="Precursor mass error tolerance type (ppm or Da). Default value: ppm")
    parser.add_argument("-E", "--precursor_error", type=float, default = 10, 
                        help="Precursor mass error tolerance. Default value: 10 ppm or 2.2 Da")
    parser.add_argument("-e", "--fragment_error", type=float, default = 10, 
                        help="Fragment mass error tolerance (in ppm). Default value: 10 ppm")
    parser.add_argument("-c", "--charge_state", type=str, default ='No', 
                        help="Charge matching requirement. Default value: Not required")
    args = parser.parse_args()
    lib_filename = args.library_filename
    msalign_filename = args.msalign_filename
    precursor_error_type = args.precursor_error_type
    precursor_error = args.precursor_error
    fragment_error = args.fragment_error
    charge_state = args.charge_state
    return lib_filename, msalign_filename, precursor_error_type, precursor_error, fragment_error, charge_state 
    

# def ms_rep_library_query(lib_filename, msalign_file_name, method_sel, pre_charge_use):
def ms_rep_library_query(lib_filename, msalign_file_name, tol_type='ppm', tol_val=None, 
                         frag_tol_val=10, pre_charge_use='no'):
    # set tol_val based on tol_type if it's not provided
    if tol_val is None:
        tol_val = 10 if tol_type == 'ppm' else 2.2                   
    
    def get_lib_representatives(lib_filename):
        # load sql library
        # Connecting to sqlite
        conn = sqlite3.connect(lib_filename) 
        # read toplib data      
        ptm_col_names = ['representative_id', 'file_id','file_name','spectrum_id','scan','precursor_mass', 'precursor_intensity', 'precursor_charge',
                         'retention_time', 'precursor_feature_id', 'proteoform_id', 'protein_accession', 'protein_description', 'first_residue', 'last_residue',
                         'protein_sequence', 'proteoform', 'e_value', 'spectrum_level_q_value', 'proteoform_level_q_value','flag']
        columns1 = list(map(lambda x: 'mass_' + str(x), np.arange(50)))
        columns2 = list(map(lambda x: 'intensity_' + str(x), np.arange(50)))
        columns3 = list(map(lambda x: 'charge_' + str(x), np.arange(50)))  
        columns4 = list(map(lambda x: 'norm_intensity_' + str(x), np.arange(50)))  
        
        query = "SELECT * FROM target_decoy_spectra_representatives"
        target_decoy_spectra_rep = pd.read_sql_query(query, con=conn)  
        target_decoy_mass_rep = pd.concat([target_decoy_spectra_rep[ptm_col_names],  target_decoy_spectra_rep[columns1]], axis=1)
        target_decoy_inte_rep  = target_decoy_spectra_rep[columns2]
        target_decoy_ch_rep  = target_decoy_spectra_rep[columns3]
        target_decoy_inte_rep_convert  = target_decoy_spectra_rep[columns4]
        # extract spectra with proteoform identifications
        target_decoy_mass_rep_ID = target_decoy_mass_rep.dropna(subset='proteoform')
        spectra_rep_index = target_decoy_mass_rep_ID.index.values.tolist()
        target_decoy_inte_rep_ID = target_decoy_inte_rep.loc[spectra_rep_index,:]
        target_decoy_inte_rep_convert_ID = target_decoy_inte_rep_convert.loc[spectra_rep_index,:]
        target_decoy_ch_rep_ID = target_decoy_ch_rep.loc[spectra_rep_index,:]
        # reindex
        target_decoy_mass_rep_ID.reset_index(inplace = True, drop=True)
        target_decoy_inte_rep_ID.reset_index(inplace = True, drop=True)
        target_decoy_inte_rep_convert_ID.reset_index(inplace = True, drop=True)
        target_decoy_ch_rep_ID.reset_index(inplace = True, drop=True)
        conn.close()
        return target_decoy_mass_rep_ID, target_decoy_inte_rep_ID, target_decoy_ch_rep_ID, target_decoy_inte_rep_convert_ID
    
    
    def ms_query_preprocessing(msalign_file_name):
        # msalign_file_name: combined msalign file
        ms_df_query, mass_df_query = mf.read_msalign(msalign_file_name)             
        # remove incorrect spectra
        drop_idx2 = ms_df_query[ms_df_query['precursor_intensity']==0].index.values
        ms_df_query.drop(drop_idx2, inplace = True)
        ms_df_query.reset_index(drop=True, inplace=True)
        drop_idx = ms_df_query[ms_df_query['num_mass']<=1].index.values
        ms_df_query.drop(drop_idx, inplace = True)
        return ms_df_query, mass_df_query
    
    
    def ms_spectrum_query_preprocess(mass_rep1):
        # ms/ms spectra preprocessing
        # extract 50 largest peaks
        # use logrithmic algorithm for intensity and normalization
        # precalcuate the error tolerance
        spec_id_list = np.unique(mass_rep1['spectrum_id'].values)
        mass1_ext = pd.DataFrame()
        for i in range(len(spec_id_list)):
            ms_ext_df = mass_rep1.loc[mass_rep1['spectrum_id']==spec_id_list[i],:].sort_values(by='intensity',ascending=False).iloc[0:50]
            mass1_ext = pd.concat([mass1_ext, ms_ext_df], axis=0)
    
        # convert to vectors
        sorted_mass_all = []
        sorted_int_all = []
        sorted_ch_all = []
        specId_list = mass1_ext['spectrum_id'].unique()
        for i in range(len(specId_list)):
            mass_v = mass1_ext.loc[mass1_ext['spectrum_id']==specId_list[i],'mass'].values.tolist()
            int_v = mass1_ext.loc[mass1_ext['spectrum_id']==specId_list[i],'intensity'].values.tolist()
            ch_v = mass1_ext.loc[mass1_ext['spectrum_id']==specId_list[i],'charge'].values.tolist()
            # sort
            sorted_mass = np.sort(mass_v)
            sorted_mass_idx = np.argsort(mass_v)
            sorted_int = np.array(int_v)[sorted_mass_idx]
            sorted_int = np.log2(sorted_int)
            sorted_int = sorted_int / np.linalg.norm(sorted_int)
            sorted_ch = np.array(ch_v)[sorted_mass_idx]
            # add tolerance
            sorted_mass_all.append(sorted_mass)
            sorted_int_all.append(sorted_int)
            sorted_ch_all.append(sorted_ch)
    
        data = {'mass': sorted_mass_all,
               'intensity': sorted_int_all,
               'charge': sorted_ch_all}
        sorted_mass_int_ch_df = pd.DataFrame(data)
        return sorted_mass_int_ch_df 
    
    
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
                th2 = 1.00235+tol2
                th3 = 1.00235-tol2
                if pre_ch_use == 'yes':
                    if pre_ch_query == pre_ch_lib[j]:
                        if abs(pre_mass_lib[j] - pre_mass_query) <= th1 or ((abs(pre_mass_lib[j] - pre_mass_query) <= th2) and (abs(pre_mass_lib[j] - pre_mass_query) >= th3)):
                            idx1.append(spectra_mass_rep_index[j])
                            if spectra_query_index[i] not in idx2:
                                idx2.append(spectra_query_index[i])
                else:
                    if abs(pre_mass_lib[j] - pre_mass_query) <= th1 or ((abs(pre_mass_lib[j] - pre_mass_query) <= th2) and (abs(pre_mass_lib[j] - pre_mass_query) >= th3)):
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
    
    
    def ms_query(target_decoy_mass_rep, target_decoy_inte_rep_convert, target_decoy_ch_rep, ms_df_query, mass_df_query, tol_type, pre_tol_val, frag_tol_val, pre_charge_use):
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
        flag_all = []
        with tqdm(total=num_spec, desc="Processing") as pbar:
            for i in range(num_spec):
                pbar.update(1)
                spectra_query = ms_df_query.loc[ms_df_query['spectrum_id']==spec_id_list[i],['spectrum_id','file_name','precursor_mass','precursor_intensity','precursor_charge']]
                spec_query_id = spectra_query['spectrum_id'].values[0]
                spec_query_filename = spectra_query['file_name'].values[0]
                spectra_query.reset_index(drop=True, inplace=True)
                mass_query = mass_df_query[mass_df_query['spectrum_id'].isin(spectra_query['spectrum_id'].values)]
                # convert to vector
                mass_query.reset_index(inplace = True, drop=True)
                mass_query_v = ms_spectrum_query_preprocess(mass_query)
                # step1: precursor match querying
                spectra_mass_rep_index = target_decoy_mass_rep.index.values
                pre_mass_lib = target_decoy_mass_rep['precursor_mass'].values
                pre_ch_lib = target_decoy_mass_rep['precursor_charge'].values
        
                spectra_query_index = spectra_query.index.values
                pre_mass_query_all = spectra_query['precursor_mass'].values
                pre_ch_query_all = spectra_query['precursor_charge'].values
                if tol_type == 'Da':
                    idx1, idx2 = precursor_mass_query_Da(pre_mass_lib, pre_ch_lib, spectra_mass_rep_index, pre_mass_query_all, pre_ch_query_all, spectra_query_index, pre_tol_val, pre_charge_use)
                else:
                    idx1, idx2 = precursor_mass_query(pre_mass_lib,pre_ch_lib, spectra_mass_rep_index, pre_mass_query_all, pre_ch_query_all,spectra_query_index, pre_tol_val, pre_charge_use)
                # step2: calculate distance
                if len(idx1):
                    idx12, res_cos = dist_cal(target_decoy_mass_rep[columns1].to_numpy(), target_decoy_inte_rep_convert[columns2].to_numpy(), target_decoy_ch_rep[columns3].to_numpy(), 
                                        mass_query_v['mass'].values[0], mass_query_v['intensity'].values[0], mass_query_v['charge'].values[0], thre, np.array(idx1), np.array(idx2), frag_tol_val)
        
                    if idx12:
                        flag=target_decoy_mass_rep['flag'].iloc[idx12[0][0]]
                        flag_all.append(flag)
                        res_cos_all.extend(res_cos)
                        idx12_all.append(idx12[0][0])
                        spectra_query_id_all.append(spec_query_id) 
                        spectra_query_filename_all.append(spec_query_filename)
    
        data = {'spectrum_id': spectra_query_id_all,
                'file_name': spectra_query_filename_all,
                'calculated_cosine_dist': res_cos_all,
                'lib_index': idx12_all,
                'flag': flag_all}
        combined_target_decoy_res = pd.DataFrame(data)
    
        combined_target_decoy_res = combined_target_decoy_res.sort_values(by='spectrum_id')
        lib_index_all = combined_target_decoy_res['lib_index'].values
        combined_target_decoy_res['representative_id'] = target_decoy_mass_rep.loc[lib_index_all,'representative_id'].values
        combined_target_decoy_res['similarity'] = 1-combined_target_decoy_res['calculated_cosine_dist'].values
        combined_target_decoy_res['file_name(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'file_name'].values
        combined_target_decoy_res['spectrum_id(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'spectrum_id'].values
        combined_target_decoy_res['scan(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'scan'].values
        combined_target_decoy_res['retention_time(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'retention_time'].values
        combined_target_decoy_res['proteoform_id(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'proteoform_id'].values
        combined_target_decoy_res['protein_accession(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'protein_accession'].values
        combined_target_decoy_res['precursor_mass(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'precursor_mass'].values
        combined_target_decoy_res['precursor_charge(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'precursor_charge'].values
        combined_target_decoy_res['precursor_feature_id(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'precursor_feature_id'].values
        combined_target_decoy_res['protein_description(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'protein_description'].values
        combined_target_decoy_res['protein_sequence(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'protein_sequence'].values
        combined_target_decoy_res['first_residue(lib)'] = target_decoy_mass_rep.loc[lib_index_all,'first_residue'].values
        combined_target_decoy_res['last_residue(lib)'] = target_decoy_mass_rep.loc[lib_index_all,'last_residue'].values
        combined_target_decoy_res['proteoform(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'proteoform'].values
        combined_target_decoy_res['e_value(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'e_value'].values
        combined_target_decoy_res['spectrum_level_q_value(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'spectrum_level_q_value'].values
        combined_target_decoy_res['proteoform_level_q_value(lib)'] = target_decoy_mass_rep.loc[lib_index_all, 'proteoform_level_q_value'].values
        return combined_target_decoy_res
    
    
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
    
    
    # load library data
    target_decoy_mass_rep, target_decoy_inte_rep ,target_decoy_ch_rep, target_decoy_inte_rep_convert = get_lib_representatives(lib_filename)
    # processing query spectrum
    ms_df_query, mass_df_query = ms_query_preprocessing(msalign_file_name)    
    start_time = datetime.now()        
    # search against library
    combined_target_decoy_res = ms_query(target_decoy_mass_rep, target_decoy_inte_rep_convert, target_decoy_ch_rep, 
                                         ms_df_query, mass_df_query, tol_type, tol_val, frag_tol_val, pre_charge_use)  
    # FDR calculation
    spec_id_1fdr = fdr_cal(combined_target_decoy_res)
    lib_search = combined_target_decoy_res.loc[(combined_target_decoy_res['spectrum_id'].isin(spec_id_1fdr)) & (combined_target_decoy_res['flag']==1),'spectrum_id'].values
    lib_ident1 = combined_target_decoy_res[combined_target_decoy_res['spectrum_id'].isin(lib_search)]
    lib_ident1 = lib_ident1.sort_values(by='spectrum_id')
    lib_ident1 = lib_ident1.drop(['lib_index', 'flag','calculated_cosine_dist'], axis=1) 
    lib_ident1.reset_index(inplace = True, drop=True)
    sim_cutoff = 0.3
    lib_ident1_sim = lib_ident1[lib_ident1['similarity']>=sim_cutoff]                     
    end_time = datetime.now()
    total_time = end_time - start_time
    print("Total runtime: ", total_time)    
   
    # report results
    curr_path = os.getcwd()
    directory = "toplib_output"
    path = os.path.join(curr_path, directory)
    try:
        os.makedirs(path, exist_ok = True)
        # print("Directory '%s' created successfully" % directory)
        w_filename = 'query_res.tsv'
        wfile = os.path.join(curr_path, directory, w_filename)
        write_result_fun(wfile, lib_ident1_sim)
        print("Query results have been stored under toplib_output folder!")
    except OSError as error:
        print("Directory '%s' can not be created" % directory)
        

if __name__ == "__main__":
    if len(sys.argv) <=2:
        print("At least three inputs, usage: python script.py <db_filename> <msalign_filename>")
        sys.exit()
    else: 
        db_filename, msalign_filename, tol_type, tol_val, frag_tol_val, charge_option = user_command()
        # file name check
        curr_path = os.getcwd()
        directory = "TopLib"
        db_filename = os.path.join(curr_path, directory, db_filename)   
        filenames = [msalign_filename, db_filename]
        if all(os.path.isfile(f) for f in filenames):
            ms_rep_library_query(db_filename, msalign_filename, tol_type, tol_val, frag_tol_val, charge_option)
        else:
            print('At least one of input files is incorrect and please check file name and path!')        

