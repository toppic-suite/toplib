import pandas as pd
import numpy as np
import sqlite3 
import random
import sys
import os
import argparse


def load_msalign_data_pred(msalign_file):
    rows = []
    dataset_id = None
    mzml_filename = None
    precursor_charge = None

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
                precursor_charge = None

                mass_all = []
                intensity_all = []
                charge_all = []
                confidence_all = []

            elif line.startswith("DATASET_ID="):
                dataset_id = line.split("=", 1)[1]

            elif line.startswith("MZML_FILE_NAME="):
                mzml_filename = line.split("=", 1)[1]

            elif line.startswith("PRECURSOR_CHARGE="):
                precursor_charge = int(line.split("=", 1)[1])

            elif line == "END IONS":
                rows.append({
                    "dataset_id": dataset_id,
                    "mzml_file_name": mzml_filename,
                    "precursor_charge": precursor_charge,
                    "mass_all": mass_all,
                    "intensity_all": intensity_all,
                    "charge_all": charge_all,
                    "confidence_all": confidence_all
                })

                # clean up explicitly
                dataset_id = None
                mzml_filename = None
                precursor_charge = None

            else:
                # Peak line
                parts = line.split()
                if len(parts) >= 4:
                    # mandatory fields
                    mass_all.append(float(parts[0]))
                    intensity_all.append(float(parts[1]))
                    charge_all.append(int(parts[2]))
                    confidence_all.append(parts[3])
                   
    msalign_df = pd.DataFrame(rows)
#     rename_map = {
#         "dataset_id": "DATASET_id",
#         "mzml_file_name": "MZML_file_name",
#         "mass_all": "mass",
#         "intensity_all": "intensity",
#         "charge_all": "charge",
#         "confidence_all": "confidence"
#     }
#     msalign_df = msalign_df.rename(columns=rename_map)
    
    return msalign_df


def gen_lib_data(tsv_filename, msalign_file):
    target_filtered3 = pd.read_csv(tsv_filename, sep='\t')
    target_filtered3 = target_filtered3.drop(columns=['flag'])
    
    ms_pred_df = load_msalign_data_pred(msalign_file)
    
    target_filtered3['mass'] = ms_pred_df['mass_all'].values
    target_filtered3['intensity'] = ms_pred_df['intensity_all'].values
    target_filtered3['charge'] = ms_pred_df['charge_all'].values
    target_filtered3 = target_filtered3.reset_index(drop=True)
    return target_filtered3


def peak_list_extract_based_precursor(ms_df):
    for i in range(len(ms_df)):
        precursor_mass = ms_df.at[i, 'precursor_mass']

        mass = np.asarray(ms_df.at[i, 'mass'])
        mask = mass <= precursor_mass

        ms_df.at[i, 'mass'] = mass[mask].tolist()
        ms_df.at[i, 'intensity'] = np.asarray(ms_df.at[i, 'intensity'])[mask].tolist()
        ms_df.at[i, 'charge'] = np.asarray(ms_df.at[i, 'charge'])[mask].tolist()
       
    return ms_df


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


def ms_spectrum_preprocess(ms_all):
    # sort, normalize
    sorted_mass_all = []
    sorted_inte_all = []
    sorted_ch_all = []
    sorted_norm_inte_all = []
    for i in range(len(ms_all)):
        # sort
        sorted_mass_idx = np.argsort(ms_all['mass'].iloc[i])
        sorted_mass = np.array(ms_all['mass'].iloc[i])[sorted_mass_idx]
        sorted_inte = np.array(ms_all['intensity'].iloc[i])[sorted_mass_idx]
        # normalization
        # pre_inte = ms_all['precursor_intensity'].iloc[i] # wrg
        # sorted_inte = sorted_inte / sorted_inte.max() # wrg
        # sorted_norm_inte = sorted_inte / np.linalg.norm(sorted_inte) # wrg
        # scale_inte = pre_inte * sorted_inte # wrg
       # log and normalization
        sorted_log_inte = np.log2(sorted_inte)
        sorted_norm_inte = sorted_log_inte / np.linalg.norm(sorted_log_inte)
        sorted_ch = np.array(ms_all['charge'].iloc[i])[sorted_mass_idx]
        
        sorted_inte_all.append(sorted_inte)
        sorted_mass_all.append(sorted_mass)
        sorted_ch_all.append(sorted_ch)
        sorted_norm_inte_all.append(sorted_norm_inte)
        
    data = {'mass': sorted_mass_all,
            'intensity': sorted_inte_all,
            'charge': sorted_ch_all,
            'norm_intensity': sorted_norm_inte_all}
    sorted_mass_inte_ch_df = pd.DataFrame(data)
    return sorted_mass_inte_ch_df


def mass_inte_ch_list_gen(sorted_mass_inte_ch_df):
    mass_rep_all = []
    inte_rep_all = []
    ch_rep_all = []
    inte_rep_convert_all = []
    for i in range(len(sorted_mass_inte_ch_df)):
        mass_rep = sorted_mass_inte_ch_df['mass'].iloc[i]
        inte_rep = sorted_mass_inte_ch_df['intensity'].iloc[i]
        ch_rep = sorted_mass_inte_ch_df['charge'].iloc[i]
        inte_rep_convert = sorted_mass_inte_ch_df['norm_intensity'].iloc[i]
        
        mass_rep_all.append(mass_rep)
        inte_rep_all.append(inte_rep)
        ch_rep_all.append(ch_rep)
        inte_rep_convert_all.append(inte_rep_convert)
    
    return mass_rep_all, inte_rep_all, ch_rep_all, inte_rep_convert_all 


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
    
    mass_matrix = spectra_mass_rep[columns1].to_numpy()
    precursor_mass = spectra_mass_rep["precursor_mass"].to_numpy()
#     n_spec, n_peaks = mass_matrix.shape  
    decoy_mass_all = []
    for j in range(len(spectra_mass_rep)):
        random.seed(j)
#         if (j % 10000)==0:
#             print(j)
            
        mass1 = mass_matrix[j]
        pre_mass1 = precursor_mass[j]  
        mass1 = mass1[~np.isnan(mass1)]
        decoy_mass = np.zeros((len(mass1),), dtype=float)

        for i in range(len(mass1)):
            tries = 0
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
                tries += 1
                if tries >= 10:
                    mass1_shift = pre_mass1
                    break
            decoy_mass[i] = mass1_shift
        decoy_mass_all.append(decoy_mass)

    # decoy_mass_all = [pad_to_50(x) for x in decoy_mass_all] 
    decoy_mass_df = pd.DataFrame(data=decoy_mass_all, columns=columns1)
    return decoy_mass_df


def td_pred_ms_lib_building(tsv_filename, msalign_file, toplib_dir):
    ms_df = gen_lib_data(tsv_filename, msalign_file)
    # filter based on precursor mass
    ms_df = peak_list_extract_based_precursor(ms_df)
    # get spectra data and 50 peaks
    ms_df, ms_org_df = ms_incorr_remove(ms_df) 
    print(f"ms_df: {len(ms_df)}")
    # convert to vector, and normalization
    sorted_mass_inte_ch_df = ms_spectrum_preprocess(ms_df)
      
    # generate a dataframe
    columns1 = list(map(lambda x: 'mass_' + str(x), np.arange(50)))
    columns2 = list(map(lambda x: 'intensity_' + str(x), np.arange(50)))
    columns3 = list(map(lambda x: 'charge_' + str(x), np.arange(50)))  
    columns4 = list(map(lambda x: 'norm_intensity_' + str(x), np.arange(50)))  
      
    mass_rep, inte_rep, ch_rep, inte_rep_convert = mass_inte_ch_list_gen(sorted_mass_inte_ch_df)
      
    mass_rep_df = pd.DataFrame(data=mass_rep, columns=columns1)
    inte_rep_df = pd.DataFrame(data=inte_rep, columns=columns2)
    charge_rep_df = pd.DataFrame(data=ch_rep, columns=columns3)
    # intensity for normalized scale
    inte_rep_cvt_df = pd.DataFrame(data=inte_rep_convert, columns=columns4)
    
    ms_rep_df = ms_df.drop(columns=['mass', 'intensity','charge'])
    ms_mass_inte_ch_rep = pd.concat([ms_rep_df, mass_rep_df, inte_rep_df, charge_rep_df, inte_rep_cvt_df], axis=1) 
    
    # build decoy 
    decoy_mass_df = decoy_mass_gen(ms_mass_inte_ch_rep)
    decoy_mass_inte_ch_rep = pd.concat([ms_rep_df, decoy_mass_df, inte_rep_df, charge_rep_df, inte_rep_cvt_df], axis=1)
    # combine decoy 
    ms_mass_inte_ch_rep['flag'] = np.ones(len(ms_mass_inte_ch_rep), dtype='int32')
    decoy_mass_inte_ch_rep['flag'] = np.zeros(len(decoy_mass_inte_ch_rep), dtype='int32')
    # combine
    target_decoy_ms_rep = pd.concat([ms_mass_inte_ch_rep, decoy_mass_inte_ch_rep], axis=0)
    target_decoy_ms_rep.reset_index(inplace = True, drop=True)
    
    # save
    os.makedirs(toplib_dir, exist_ok = True)
    
    lib_name = os.path.splitext(os.path.basename(msalign_file))[0] + '_representatives'
    lib_filename = lib_name + '.db'
    wfile = os.path.join(toplib_dir, lib_filename)
    
    conn = sqlite3.connect(wfile)
    target_decoy_ms_rep.to_sql(
        name="target_decoy_spectra_representatives",
        con=conn,
        index=False,
        if_exists="append",
        chunksize=1000,     
        method="multi"      
    )
    
    conn.close()
    print("Toplib library has been built!")
    print(f"Toplib library saved to {wfile}")
    
 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build top-down MS spectral library for TD-Pred spectra")
    parser.add_argument(
        "-m", "--msalign",
        required=True,
        help="Input msalign file"
    )

    parser.add_argument(
        "-t", "--tsv",
        required=True,
        help="Input TSV file"
    )

    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory for TopLib database"
    )

    args = parser.parse_args()
    if not os.path.isfile(args.msalign):
        print(f"Error: msalign file not found: {args.msalign}")
        sys.exit(1)

    if not os.path.isfile(args.tsv):
        print(f"Error: TSV file not found: {args.tsv}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    td_pred_ms_lib_building(
        tsv_filename=args.tsv,
        msalign_file=args.msalign,
        toplib_dir=args.outdir
    )