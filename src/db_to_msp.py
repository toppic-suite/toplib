import sys
import os
import pandas as pd
import numpy as np
import sqlite3 
from tqdm import tqdm
import re
from db_add_file import get_msfilename

def get_spectra_representatives(toplib_filepath):
    conn = sqlite3.connect(toplib_filepath)
    query = "SELECT * FROM target_decoy_spectra_representatives"
    target_decoy_spectra_rep = pd.read_sql_query(query, con=conn)
    target_spectra_rep = target_decoy_spectra_rep[target_decoy_spectra_rep['flag']==1]
    conn.close()
    return target_spectra_rep


def extract_acetyl_residue(seq):
    if seq.startswith('[Acetyl]'):
        match = re.search(r'\[Acetyl\]-([A-Z])', seq)
        if match:
            return match.group(1)  # the first amino acid after Acetyl
    return None


# extract previous/next residue
def remove_residue(proteoform):
    previous_residue_idx = proteoform.find('.')
    next_residue_idx = proteoform.rfind('.')
    previous_residue = proteoform[0:previous_residue_idx]
    next_residue = proteoform[next_residue_idx+1:]
    proteoform_remove_residue = proteoform[previous_residue_idx+1:next_residue_idx]
    return pd.Series([proteoform_remove_residue, previous_residue, next_residue], index=['proteoform_remove_residue', 'previous_residue', 'next_residue'])


def ptm_pos_val(proteoform_remove_residue, fixed_ptms):    
    ptm_list = []
    # acetyle
    if 'Acetyl' in proteoform_remove_residue:
        ptm_val = 'Acetyl'
        ptm_pos = 0 
        ptm_residue = extract_acetyl_residue(proteoform_remove_residue)
        ptm_list.append((ptm_pos, ptm_residue, ptm_val))
    if fixed_ptms:
        ptms = fixed_ptms.split(';')
        for ptm in ptms:
            ptm_val = ptm.split(':')[0]
            ptm_pos = ptm.split(':')[1].replace('[','').replace(']','')
            ptm_pos = int(ptm_pos)-1
            if 'Carbamidomethylation' in ptm_val: 
                ptm_list.append((ptm_pos, 'C', 'CAM'))
    return ptm_list


def mod_pos_val(mod_val):   
    if mod_val:
        shift_val = mod_val.split(':')[0]
        shift_pos = mod_val.split(':')[1].replace('[','').replace(']','')    
        shift_start_pos = int(shift_pos.split('-')[0])-1
        shift_termin_pos = int(shift_pos.split('-')[-1])-1
    else:
        shift_start_pos, shift_termin_pos, shift_val = None, None, None
    return shift_start_pos, shift_termin_pos, shift_val


def write_msp_format(target_spectra_rep1, output_filename):
    curr_path = os.getcwd()
    directory = "TopLib"
    wfile = os.path.join(curr_path, directory, output_filename)   
    # extract spectral data and convert to a msalign file
    columns1 = list(map(lambda x: 'mass_' + str(x), np.arange(50)))
    columns2 = list(map(lambda x: 'intensity_' + str(x), np.arange(50)))
    columns3 = list(map(lambda x: 'charge_' + str(x), np.arange(50))) 
    # writing to a msalign file
    with tqdm(total=len(target_spectra_rep1), desc="Processing") as pbar:
        with open(wfile, 'w') as as_file:
            for i in range(len(target_spectra_rep1)):
                pbar.update(1)
                mass_v = target_spectra_rep1[columns1].iloc[i].values
                mass_v = mass_v[~np.isnan(mass_v)]
                inte_v = target_spectra_rep1[columns2].iloc[i].values
                inte_v = inte_v[~np.isnan(inte_v)]
                ch_v = target_spectra_rep1[columns3].iloc[i].values
                ch_v = ch_v[~np.isnan(ch_v)]
                # modifications
                proteoform_remove_residue = target_spectra_rep1['proteoform_remove_residue'].iloc[i]
                fixed_ptms = target_spectra_rep1['fixed_ptms'].iloc[i]
                mod_val = target_spectra_rep1['unexpected_modifications'].iloc[i]    
                ptm_list = ptm_pos_val(proteoform_remove_residue, fixed_ptms)
                shift_start_pos, shift_termin_pos, shift_val = mod_pos_val(mod_val)        
                ptm_str = "".join("({},{},{})".format(p, q, r) for p, q, r in ptm_list)
                # write
                as_file.write("Name: {}/{}_{}{}\n".format(target_spectra_rep1['protein_sequence'].iloc[i], 
                                                        target_spectra_rep1['precursor_charge'].iloc[i],
                                                        len(ptm_list),
                                                        ptm_str)
                             )
                as_file.write("MW: {}\n".format(target_spectra_rep1['precursor_mass'].iloc[i]))
                as_file.write((
                               "Comment: Spec=Consensus Mods={}{} Fullname={} Charge={} Protein=\"{}\" Description=\"{}\" " 
                               "Retention time={} First residue={} Last residue={} Proteoform={} Previous residue={} Next residue={} " 
                               "Mass shift position={} {} Shift value={} "
                               "E-value={} Spectrum-level-Q-value={} Proteoform-level-Q-value={}\n").format(
                                len(ptm_list), ptm_str,
                                target_spectra_rep1['sequence_fullname'].iloc[i],
                                target_spectra_rep1['precursor_charge'].iloc[i], 
                                target_spectra_rep1['protein_accession'].iloc[i],
                                target_spectra_rep1['protein_description'].iloc[i],
                                target_spectra_rep1['retention_time'].iloc[i],
                                int(target_spectra_rep1['first_residue'].iloc[i]-1),
                                int(target_spectra_rep1['last_residue'].iloc[i]-1),
                                target_spectra_rep1['proteoform'].iloc[i],
                                target_spectra_rep1['previous_residue'].iloc[i],
                                target_spectra_rep1['next_residue'].iloc[i],
                                shift_start_pos, shift_termin_pos, shift_val,
                                target_spectra_rep1['e_value'].iloc[i],
                                target_spectra_rep1['spectrum_level_q_value'].iloc[i],
                                target_spectra_rep1['proteoform_level_q_value'].iloc[i])
                )
                as_file.write("Num peaks: {}\n".format(len(mass_v)))
                for m in range(len(mass_v)):
                    as_file.write("{:.5f}\t{:.2f}\t{:d}\n".format(mass_v[m], inte_v[m], int(ch_v[m])))
                as_file.write('\n\n')                    
        
        
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_toplib_file>")
        sys.exit()
    else:
        curr_path = os.getcwd()
        directory = "TopLib"
        lib_filename = sys.argv[1] # library file (.db)
        db_filename = os.path.join(curr_path, directory, lib_filename)             
        # output_filename = sys.argv[2] # only provide an output filename  
        temp_filename = get_msfilename(db_filename)
        output_filename = temp_filename + "_nist.msp"        
        # file name check
        if os.path.isfile(db_filename):
            # get spectra data
            target_spectra_rep = get_spectra_representatives(db_filename)
            # add a new proteoform column after removing residues 
            target_spectra_rep[['proteoform_remove_residue','previous_residue','next_residue']] = target_spectra_rep['proteoform'].apply(remove_residue)
            # add previous/next residues to protein sequence
            target_spectra_rep['sequence_fullname'] = target_spectra_rep.apply(lambda row: f"{row['previous_residue']}.{row['protein_sequence']}.{row['next_residue']}", axis=1)
            # write to NIST format
            write_msp_format(target_spectra_rep, output_filename)
        else:
            print('No such library file exists in the TopLib folder!')
