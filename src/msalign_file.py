# -*- coding: utf-8 -*-
"""
Created on Mon April 1st 2024

@author: kunli

This function runs a msalign file and generate MS/MS spectra and fragmental mass data
"""
import sys
import pandas as pd


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

    # extract the data of mass and intensity from the list
    mass_val_list_tot = []
    intensity_val_list_tot = []
    charge_val_list_tot = []
    mass_num_list_tot = []

    mass_df_spec_list = []
    mass_df_mass_list = []
    mass_df_inte_list = []
    mass_df_charge_list = []

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
        precursor_mass_line = [ss for ss in spectrum_data if "PRECURSOR_MASS" in ss][0]
        precursor_intensity_line = [kk for kk in spectrum_data if "PRECURSOR_INTENSITY" in kk][0]
        precursor_charge_line = [tt for tt in spectrum_data if "PRECURSOR_CHARGE" in tt][0]
        precursor_feature_id_line = [tt for tt in spectrum_data if "PRECURSOR_FEATURE_ID" in tt][0]
        spec_id_line = [tt for tt in spectrum_data if "SPECTRUM_ID=" in tt][0]
        scan_num_line = [vv for vv in spectrum_data if "SCANS=" in vv][0]
        retention_time_line = [t for t in spectrum_data if "RETENTION_TIME" in t][0]

        file_name = file_name_line.split('=')[1]
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

        precursor_mass_int = round(precursor_mass)

        # skip spectrum with a zero precursor mass
        if precursor_mass == 0.0:
            continue

        file_name_list.append(file_name)
        precursor_mass_list.append(precursor_mass)
        precursor_mass_int_list.append(precursor_mass_int)
        precursor_intensity_list.append(precursor_intensity)
        precursor_charge_list.append(precursor_charge)
        precursor_feature_id_list.append(precursor_feature_id)
        spec_id_list.append(spec_id)
        scan_num_list.append(scan_num)
        retention_time_list.append(retention_time)


        mass_val_list = []
        intensity_val_list = []
        charge_val_list = []

        for j in range(len(mass_list)):
            ms2 = mass_list[j].split('\t')
            mass_val = float(ms2[0])
            intensity_val = float(ms2[1])
            charge_val = int(ms2[2])  # charge value
            if mass_val > precursor_mass:
                continue

            mass_val_list.append(mass_val)
            intensity_val_list.append(intensity_val)
            charge_val_list.append(charge_val)

            mass_df_spec_list.append(spec_id)
            mass_df_mass_list.append(mass_val)
            mass_df_inte_list.append(intensity_val)
            mass_df_charge_list.append(charge_val)

        mass_val_list_tot.append(mass_val_list)
        intensity_val_list_tot.append(intensity_val_list)
        charge_val_list_tot.append(charge_val_list)
        mass_num_list_tot.append(len(mass_val_list))

    data_1 = {'file_name': file_name_list,
              'spectrum_id': spec_id_list,
              'scan': scan_num_list,
              'retention_time': retention_time_list,
              'precursor_mass': precursor_mass_list,
              'precursor_mass_round': precursor_mass_int_list,
              'precursor_intensity': precursor_intensity_list,
              'precursor_charge': precursor_charge_list,
              'precursor_feature_id': precursor_feature_id_list,
              'num_mass': mass_num_list_tot,
              'mass': mass_val_list_tot,
              'intensity': intensity_val_list_tot,
              'charge': charge_val_list_tot
              }
    msalign_df_1 = pd.DataFrame(data_1)
    print(msalign_df_1)

    data_2 = {'spectrum_id': mass_df_spec_list,
              'mass': mass_df_mass_list,
              'intensity': mass_df_inte_list,
              'charge': mass_df_charge_list}
    msalign_df_2 = pd.DataFrame(data_2)
    print(msalign_df_2)

    return msalign_df_1, msalign_df_2


def write_malign(output_filename, ms_df):
    with open(output_filename, 'w') as as_file:
        for ss in range(len(ms_df)):
            as_file.write('BEGIN IONS' + "\n")
            as_file.write("FILE_NAME={}\n".format(ms_df['file_name'].iloc[ss]))
            as_file.write("SPECTRUM_ID={}\n".format(ms_df['spectrum_id'].iloc[ss]))
            as_file.write("SCANS={}\n".format(ms_df['scan'].iloc[ss]))
            as_file.write("RETENTION_TIME={}\n".format(ms_df['retention_time'].iloc[ss]))
            as_file.write("PRECURSOR_CHARGE={}\n".format(ms_df['precursor_charge'].iloc[ss]))
            as_file.write("PRECURSOR_MASS={:.5f}\n".format(ms_df['precursor_mass'].iloc[ss]))
            as_file.write("PRECURSOR_INTENSITY={:.2f}\n".format(ms_df['precursor_intensity'].iloc[ss]))
            as_file.write("PRECURSOR_FEATURE_ID={:d}\n".format(ms_df['precursor_feature_id'].iloc[ss]))

            mass_list = ms_df['mass'].iloc[ss]
            inte_list = ms_df['intensity'].iloc[ss]
            charge_list = ms_df['charge'].iloc[ss]
            for m in range(len(mass_list)):
                as_file.write("{:.5f}\t{:.2f}\t{}\n".format(mass_list[m], inte_list[m], charge_list[m]))
            as_file.write('END IONS' + "\n")
            as_file.write("\n")


def write_one_spectrum(as_file, ms_row):
    as_file.write('BEGIN IONS' + "\n")
    as_file.write("FILE_NAME={}\n".format(ms_row['file_name']))
    as_file.write("SPECTRUM_ID={}\n".format(ms_row['spectrum_id']))
    as_file.write("SCANS={}\n".format(ms_row['scan']))
    as_file.write("RETENTION_TIME={}\n".format(ms_row['retention_time']))
    as_file.write("PRECURSOR_CHARGE={}\n".format(ms_row['precursor_charge']))
    as_file.write("PRECURSOR_MASS={:.5f}\n".format(ms_row['precursor_mass']))
    as_file.write("PRECURSOR_INTENSITY={:.2f}\n".format(ms_row['precursor_intensity']))
    as_file.write("PRECURSOR_FEATURE_ID={:d}\n".format(ms_row['precursor_feature_id']))

    mass_list = ms_row['mass']
    inte_list = ms_row['intensity']
    charge_list = ms_row['charge']
    for m in range(len(mass_list)):
        as_file.write("{:.5f}\t{:.2f}\t{}\n".format(mass_list[m], inte_list[m], charge_list[m]))
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
