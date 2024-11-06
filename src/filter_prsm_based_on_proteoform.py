

def process_tsv(df, form_df, output_file_path):
    proteoform_id_set = set(form_df['Proteoform ID'])
    df_filtered = df[df['Proteoform ID'].isin(proteoform_id_set)]
    df_filtered = df_filtered[df_filtered['Proteoform-level Q-value']<=0.01]    
    # Write the result back to a new TSV file
    df_filtered.to_csv(output_file_path, sep='\t', index=False)
    return df_filtered


