import sys
import pandas as pd
from tqdm import tqdm


def read_msalign(file_msalign):
    msalign = []
    with open(file_msalign, 'r', encoding='UTF-8') as f:
        for line in f.readlines():
            msalign.append(line.rstrip('\n'))

    index_begins = []
    for line in range(len(msalign)):
        if "BEGIN IONS" in msalign[line]:
            index_begins.append(line)

    index_ends = []
    for line in range(len(msalign)):
        if "END IONS" in msalign[line]:
            index_ends.append(line)

    spectrum_num = len(index_begins)

    # extract the fields of precursor mass, intensity, ID and Scan, etc.
    file_name_list = []
    precursor_mass_list = []
    precursor_mass_int_list = []
    precursor_intensity_list = []
    precursor_charge_list = []
    precursor_feature_id_list = []
    spec_id_list = []
    scan_num_list = []
    retention_time_list = []
    title_list = []
    
    ms_level_list = []
    ms_one_id_list = []
    ms_one_scan_list = []
    precursor_wins_begin_list = []
    precursor_wins_end_list = []
    precursor_mz_list = []
    activation_list = []
    
    # extract the data of mass and intensity from the list
    mass_val_list_tot = []
    intensity_val_list_tot = []
    charge_val_list_tot = []
    confidence_val_list_tot = []
    mass_num_list_tot = []

    mass_df_spec_list = []
    mass_df_mass_list = []
    mass_df_inte_list = []
    mass_df_charge_list = []
    mass_df_confidence_list = []

    for i in range(spectrum_num):
        l_start = index_begins[i]
        l_end = index_ends[i]
        spectrum_data = msalign[l_start:l_end + 1]

        mass_row_indexes = [row for row in range(len(spectrum_data)) if spectrum_data[row][0].isdigit()]

        # skip spectra without masses
        if len(mass_row_indexes) == 0:
            continue

        mass_list = [spectrum_data[i] for i in mass_row_indexes]

        file_name_line = [ss for ss in spectrum_data if "FILE_NAME" in ss][0]
        precursor_mz_line = [ss for ss in spectrum_data if "PRECURSOR_MZ" in ss][0]
        precursor_mass_line = [ss for ss in spectrum_data if "PRECURSOR_MASS" in ss][0]
        precursor_intensity_line = [kk for kk in spectrum_data if "PRECURSOR_INTENSITY" in kk][0]
        precursor_charge_line = [tt for tt in spectrum_data if "PRECURSOR_CHARGE" in tt][0]
        precursor_feature_id_line = [tt for tt in spectrum_data if "PRECURSOR_FEATURE_ID" in tt][0]
        spec_id_line = [tt for tt in spectrum_data if "SPECTRUM_ID=" in tt][0]
        scan_num_line = [vv for vv in spectrum_data if "SCANS=" in vv][0]
        retention_time_line = [t for t in spectrum_data if "RETENTION_TIME" in t][0]
        title_line = [t for t in spectrum_data if "TITLE" in t][0]
        
        level_line = [t for t in spectrum_data if "LEVEL" in t][0]
        ms_one_id_line = [t for t in spectrum_data if "MS_ONE_ID" in t][0]
        ms_one_scan_line = [t for t in spectrum_data if "MS_ONE_SCAN" in t][0]
        precursor_wins_begin_line = [t for t in spectrum_data if "PRECURSOR_WINDOW_BEGIN" in t][0]
        precursor_wins_end_line = [t for t in spectrum_data if "PRECURSOR_WINDOW_END" in t][0]
        activation_line = [t for t in spectrum_data if "ACTIVATION" in t][0]
        

        file_name = file_name_line.split('=')[1]
        precursor_mz_info = precursor_mz_line.split('=')[1]
        precursor_mz = float(precursor_mz_info.split(':')[0])
        precursor_mass_info = precursor_mass_line.split('=')[1]
        precursor_mass = float(precursor_mass_info.split(':')[0])
        precursor_intensity_info = precursor_intensity_line.split('=')[1]
        precursor_intensity = float(precursor_intensity_info.split(':')[0])
        precursor_charge_info = precursor_charge_line.split('=')[1]
        precursor_charge = int(precursor_charge_info.split(':')[0])
        precursor_feature_id_info = precursor_feature_id_line.split('=')[1]
        precursor_feature_id = int(precursor_feature_id_info.split(':')[0])
        spec_id = int(spec_id_line.split('=')[1])
        scan_num = int(scan_num_line.split('=')[1])
        retention_time = float(retention_time_line.split('=')[1])
        title = title_line.split('=')[1]
        
        ms_level = int(level_line.split('=')[1])
        ms_one_ids = int(ms_one_id_line.split('=')[1])
        ms_one_scans = int(ms_one_scan_line.split('=')[1])
        precursor_wins_begin_lines = float(precursor_wins_begin_line.split('=')[1])
        precursor_wins_end_lines = float(precursor_wins_end_line.split('=')[1])   
        activation_lines = activation_line.split('=')[1]   

        precursor_mass_int = round(precursor_mass)

        # skip spectrum with a zero precursor mass
        if precursor_mass == 0.0:
            continue

        file_name_list.append(file_name)
        precursor_mz_list.append(precursor_mz)
        precursor_mass_list.append(precursor_mass)
        precursor_mass_int_list.append(precursor_mass_int)
        precursor_intensity_list.append(precursor_intensity)
        precursor_charge_list.append(precursor_charge)
        precursor_feature_id_list.append(precursor_feature_id)
        spec_id_list.append(spec_id)
        scan_num_list.append(scan_num)
        retention_time_list.append(retention_time)
        title_list.append(title)
        
        
        ms_level_list.append(ms_level)
        ms_one_id_list.append(ms_one_ids)
        ms_one_scan_list.append(ms_one_scans)
        precursor_wins_begin_list.append(precursor_wins_begin_lines)
        precursor_wins_end_list.append(precursor_wins_end_lines)
        activation_list.append(activation_lines)
        
        mass_val_list = []
        intensity_val_list = []
        charge_val_list = []
        confidence_val_list = []

        for j in range(len(mass_list)):
            ms2 = mass_list[j].split('\t')
            mass_val = float(ms2[0])
            intensity_val = float(ms2[1])
            charge_val = int(ms2[2])  # charge value
            confidence_val =  float(ms2[3])
            if mass_val > precursor_mass:
                continue

            mass_val_list.append(mass_val)
            intensity_val_list.append(intensity_val)
            charge_val_list.append(charge_val)
            confidence_val_list.append(confidence_val)

            mass_df_spec_list.append(spec_id)
            mass_df_mass_list.append(mass_val)
            mass_df_inte_list.append(intensity_val)
            mass_df_charge_list.append(charge_val)
            mass_df_confidence_list.append(confidence_val)

        mass_val_list_tot.append(mass_val_list)
        intensity_val_list_tot.append(intensity_val_list)
        charge_val_list_tot.append(charge_val_list)
        confidence_val_list_tot.append(confidence_val_list)
        
        mass_num_list_tot.append(len(mass_val_list))

    data_1 = {'file_name': file_name_list,
              'spectrum_id': spec_id_list,
              'title': title_list,
              'scan': scan_num_list,
              'retention_time': retention_time_list,          
              'ms_level':  ms_level_list,
              'ms_one_id': ms_one_id_list,
              'ms_one_scan': ms_one_scan_list,
              'precursor_window_begin': precursor_wins_begin_list,       
              'precursor_window_end': precursor_wins_end_list,  
              'activation': activation_list,
              'precursor_mz': precursor_mz_list,
              'precursor_mass': precursor_mass_list,
              'precursor_mass_round': precursor_mass_int_list,
              'precursor_intensity': precursor_intensity_list,
              'precursor_charge': precursor_charge_list,
              'precursor_feature_id': precursor_feature_id_list,
              'num_mass': mass_num_list_tot,
              'mass': mass_val_list_tot,
              'intensity': intensity_val_list_tot,
              'charge': charge_val_list_tot,
              'confidence': confidence_val_list_tot
              }
    msalign_df_1 = pd.DataFrame(data_1)
    print(msalign_df_1)

    data_2 = {'spectrum_id': mass_df_spec_list,
              'mass': mass_df_mass_list,
              'intensity': mass_df_inte_list,
              'charge': mass_df_charge_list,
              'confidence': mass_df_confidence_list}
    msalign_df_2 = pd.DataFrame(data_2)
    print(msalign_df_2)

    return msalign_df_1, msalign_df_2


def write_msalign(output_filename, ms_df, flag=0):
    # writing to a msalign file
    with tqdm(total=len(ms_df), desc="Processing") as pbar:
        with open(output_filename, 'w') as as_file:
            for ss in range(len(ms_df)):
                pbar.update(1)
                as_file.write('BEGIN IONS' + "\n")
                as_file.write("FILE_NAME={}\n".format(ms_df['file_name'].iloc[ss]))
                as_file.write("SPECTRUM_ID={}\n".format(ms_df['spectrum_id'].iloc[ss]))
                as_file.write("TITLE={}\n".format(ms_df['title'].iloc[ss]))
                as_file.write("SCANS={}\n".format(ms_df['scan'].iloc[ss]))
                as_file.write("RETENTION_TIME={}\n".format(ms_df['retention_time'].iloc[ss]))         
                as_file.write("LEVEL={}\n".format(ms_df['ms_level'].iloc[ss]))
                as_file.write("MS_ONE_ID={}\n".format(ms_df['ms_one_id'].iloc[ss]))
                as_file.write("MS_ONE_SCAN={}\n".format(ms_df['ms_one_scan'].iloc[ss]))
                as_file.write("PRECURSOR_WINDOW_BEGIN={}\n".format(ms_df['precursor_window_begin'].iloc[ss]))
                as_file.write("PRECURSOR_WINDOW_END={}\n".format(ms_df['precursor_window_end'].iloc[ss]))
                as_file.write("ACTIVATION={}\n".format(ms_df['activation'].iloc[ss]))
                as_file.write("PRECURSOR_MZ={:.5f}\n".format(ms_df['precursor_mz'].iloc[ss]))         
                as_file.write("PRECURSOR_CHARGE={}\n".format(int(ms_df['precursor_charge'].iloc[ss])))
                as_file.write("PRECURSOR_MASS={:.5f}\n".format(ms_df['precursor_mass'].iloc[ss]))
                as_file.write("PRECURSOR_INTENSITY={:.2f}\n".format(ms_df['precursor_intensity'].iloc[ss]))
                as_file.write("PRECURSOR_FEATURE_ID={:d}\n".format(int(ms_df['precursor_feature_id'].iloc[ss])))            
                if flag:
                    as_file.write("PROTEIN_ACCESSION={}\n".format(ms_df['protein_accession'].iloc[ss]))
                    as_file.write("FIRST_RESIDUE={}\n".format(int(ms_df['first_residue'].iloc[ss])))
                    as_file.write("LAST_RESIDUE={}\n".format(int(ms_df['last_residue'].iloc[ss])))
                    as_file.write("PREVIOUS_RESIDUE={}\n".format(ms_df['previous_residue'].iloc[ss]))
                    as_file.write("NEXT_RESIDUE={}\n".format(ms_df['next_residue'].iloc[ss]))
                    as_file.write("SEQUENCE={}\n".format(ms_df['protein_sequence'].iloc[ss]))
                    # write additional proteoform info.
                    mod_val = ms_df['unexpected_modifications'].iloc[ss]
                    fixed_ptms = ms_df['fixed_ptms'].iloc[ss]
                    if mod_val:
                        shift_val = mod_val.split(':')[0]
                        shift_pos = mod_val.split(':')[1].replace('[','').replace(']','')    
                        shift_start_pos = shift_pos.split('-')[0]
                        shift_termin_pos = shift_pos.split('-')[-1]
                        as_file.write("ANNOTATION=SHIFT {} {} {}\n".format(shift_start_pos, shift_termin_pos, shift_val)) 
                    
                    ptm_val = []
                    ptm_pos = []
                    if fixed_ptms:
                        ptms = fixed_ptms.split(';')
                        for ptm in ptms:
                            ptm_val.append(ptm.split(':')[0])
                            ptm_pos.append(ptm.split(':')[1].replace('[','').replace(']',''))
                        ptm_data = list(zip(ptm_val, ptm_pos))
                        for k in range(len(ptm_data)):
                            as_file.write("ANNOTATION=PTM {} {} {}\n".format(ptm_data[k][1], ptm_data[k][1], ptm_data[k][0])) 
                
                    if (not mod_val) and (not fixed_ptms):
                        as_file.write('ANNOTATION=' + "" + "\n")
                
                mass_list = ms_df['mass'].iloc[ss]
                inte_list = ms_df['intensity'].iloc[ss]
                charge_list = ms_df['charge'].iloc[ss]
                confidence_list = ms_df['confidence'].iloc[ss]
                for m in range(len(mass_list)):
                    as_file.write("{:.5f}\t{:.2f}\t{}\t{:.2f}\n".format(mass_list[m], inte_list[m], charge_list[m], confidence_list[m]))
                as_file.write('END IONS' + "\n")
                as_file.write("\n")


def write_one_spectrum(as_file, ms_row):
    as_file.write('BEGIN IONS' + "\n")
    as_file.write("FILE_NAME={}\n".format(ms_row['file_name']))
    as_file.write("SPECTRUM_ID={}\n".format(ms_row['spectrum_id']))
    as_file.write("TITLE={}\n".format(ms_row['title']))
    as_file.write("SCANS={}\n".format(ms_row['scan']))
    as_file.write("RETENTION_TIME={}\n".format(ms_row['retention_time']))
    as_file.write("LEVEL={}\n".format(ms_row['ms_level']))
    as_file.write("MS_ONE_ID={}\n".format(ms_row['ms_one_id']))
    as_file.write("MS_ONE_SCAN={}\n".format(ms_row['ms_one_scan']))
    as_file.write("PRECURSOR_WINDOW_BEGIN={}\n".format(ms_row['precursor_window_begin']))
    as_file.write("PRECURSOR_WINDOW_END={}\n".format(ms_row['precursor_window_end']))
    as_file.write("ACTIVATION={}\n".format(ms_row['activation']))
    as_file.write("PRECURSOR_MZ={:.5f}\n".format(ms_row['precursor_mz']))         
    as_file.write("PRECURSOR_CHARGE={}\n".format(ms_row['precursor_charge']))
    as_file.write("PRECURSOR_MASS={:.5f}\n".format(ms_row['precursor_mass']))
    as_file.write("PRECURSOR_INTENSITY={:.2f}\n".format(ms_row['precursor_intensity']))
    as_file.write("PRECURSOR_FEATURE_ID={:d}\n".format(ms_row['precursor_feature_id']))

    mass_list = ms_row['mass']
    inte_list = ms_row['intensity']
    charge_list = ms_row['charge']
    confidence_list = ms_row['confidence']
    for m in range(len(mass_list)):
        as_file.write("{:.5f}\t{:.2f}\t{}\t{:.2f}\n".format(mass_list[m], inte_list[m], charge_list[m], confidence_list[m]))


    as_file.write('END IONS' + "\n")
    as_file.write("\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file_path> <output_file_path>")
    else:
        df1, df2 = read_msalign(sys.argv[1])
        # write_malign(sys.argv[2], df1)
        output_file = open(sys.argv[2], 'w')
        write_one_spectrum(output_file, df1.iloc[0])
        output_file.close()
