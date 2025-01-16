from pyteomics import mzml, mgf
import pandas as pd
import numpy as np
import msalign_file as mf
import os
import sys


def mzML_ms_extract(filename):
    # extract valid mzml data
    spectrum = []
    with mzml.read(filename) as spectra:
        for mzml_in in spectra:
            spectrum.append(mzml_in)

    spectrum2 = []
    for ii in range(len(spectrum)):
        if 'precursorList' in list(spectrum[ii].keys()) and 'peak intensity' in list(spectrum[ii]['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0].keys()):
            s1={'params': {'title': spectrum[ii]['spectrum title'].split(',')[0],
              'scans': int(spectrum[ii]['spectrum title'].split(',')[1][spectrum[ii]['spectrum title'].split(',')[1].find('scan='):-1][5:]),
              'rtinseconds': float(spectrum[ii]['scanList']['scan'][0]['scan start time']) * 60,
              'pepmass_mz': float(spectrum[ii]['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']), 
              'intensity': float(str(spectrum[ii]['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['peak intensity'])),
              'charge': spectrum[ii]['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']},
              'm/z array': spectrum[ii]['m/z array'],
              'intensity array': spectrum[ii]['intensity array'] 
            }
            spectrum2.append(s1)
    return spectrum2

    
def matchLib_mzML(msalign_filename, mzml_filepath):
    # read ms lib data
    ms_df, mass_df = mf.read_msalign(msalign_filename)
        
    ms_filenames = sorted(ms_df['file_name'].unique())
    mzml_combine = pd.DataFrame()        
    for ff in range(len(ms_filenames)):
        ms_df_file = ms_df[ms_df['file_name']==ms_filenames[ff]]
        filename = ms_filenames[ff].split('/')[-1]
        # process mzml file
        normalized_path = os.path.normpath(mzml_filepath)
        mzml_filename = os.path.join(normalized_path, filename)
        spectrum2 = mzML_ms_extract(mzml_filename)
        ms_df2 = pd.DataFrame(data = spectrum2)
        
        # format conversion
        ms_df2['m/z array'] = ms_df2['m/z array'].apply(lambda x: np.round(np.array(x), 8).tolist())
        ms_df2['intensity array'] = ms_df2['intensity array'].apply(lambda x: np.round(np.array(x), 8).tolist())
        expanded_df = pd.json_normalize(ms_df2['params'])
        result_df = pd.concat([ms_df2.drop(columns=['params']), expanded_df], axis=1)
        
        # find those spectra stored in the library and retain the same order as ms_df
        mzml_ext = ms_df_file[['scan']].merge(result_df, left_on='scan', right_on='scans', how='left')
        mzml_ext['filename'] = filename
        mzml_ext.drop(columns=['scan'], axis=1, inplace=True)
        # add original index
        mzml_ext['original_index'] = ms_df_file.index
        # combine
        mzml_combine = pd.concat([mzml_combine, mzml_ext], axis=0)
    
    mzml_combine.sort_values(by='original_index', inplace=True)
    mzml_combine.drop(columns=['original_index'], axis=1, inplace=True)
    # reset index
    mzml_combine.reset_index(inplace=True, drop=True)
    return mzml_combine


def mgf_write(mzml_combine, mgf_filename):
    # extract mzml data and convert to MGF format
    mzml_spectrum = []
    for i in range(len(mzml_combine)):
        s1={'params': {'title': mzml_combine['title'].iloc[i],
          'pepmass': (mzml_combine['pepmass_mz'].iloc[i], ),
          'scans': mzml_combine['scans'].iloc[i],
          'rtinseconds': mzml_combine['rtinseconds'].iloc[i],
          'charge': mzml_combine['charge'].iloc[i]},
          'm/z array': mzml_combine['m/z array'].iloc[i],
          'intensity array': mzml_combine['intensity array'].iloc[i]
        }
        mzml_spectrum.append(s1)
        # write it to a mgf file
    with open(mgf_filename, 'w') as out:
        for spec in mzml_spectrum:
            mgf.write([spec], out)
            
         

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_msalign_file> <input_mzML_path> <input_mgf_wfile>")
        sys.exit()
    else:
        msalign_filename = sys.argv[1]
        mzml_filepath = sys.argv[2] # input the mzML file path only (e.g.,D:\sw480_experiment_F)
        mgf_filename = sys.argv[3] 
        mzml_combine = matchLib_mzML(msalign_filename, mzml_filepath)
        mgf_write(mzml_combine, mgf_filename)
        
    
