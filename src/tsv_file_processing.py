import os
import sys
import remove_comment_lines as rcl
import remove_duplicated as rd
import remove_inconsistent_prsm as ri_prsm
import filter_prsm_based_on_proteoform as filter_prsm

def prsm_file_processing(prsm_file_path, proteoform_file_path):
    
    # remove comments for prsm
    prefix, ext = os.path.splitext(prsm_file_path)
    prsm_file_remove_comment = prefix + "_remove_comment.tsv"
    rcl.remove_first_29_lines(prsm_file_path, prsm_file_remove_comment)
    
    # remove comments for proteoform
    prefix, ext = os.path.splitext(proteoform_file_path)
    form_file_remove_comment = prefix + "_remove_comment.tsv"
    rcl.remove_first_29_lines(proteoform_file_path, form_file_remove_comment)
    
    # remove duplicate proteoform
    df_form_remove_dup = rd.process_tsv(form_file_remove_comment)
    # remove inconsistent prsm
    df_inconsistent = ri_prsm.process_tsv(prsm_file_remove_comment)
    # filter based on proteoform 
    prefix, ext = os.path.splitext(prsm_file_path)
    tsv_prsm_filtered = prefix + "_filtered.tsv"
    df_filtered = filter_prsm.process_tsv(df_inconsistent, df_form_remove_dup, tsv_prsm_filtered)    
    

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_prsm_file_path> <input_proteoform_file_path>")
    else:
        df = prsm_file_processing(sys.argv[1], sys.argv[2])
        