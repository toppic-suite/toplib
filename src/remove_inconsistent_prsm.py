import pandas as pd


def process_tsv(input_file_path):
    # Load the TSV file into a DataFrame
    df = pd.read_csv(input_file_path, sep='\t')
    # Sort the DataFrame by 'Proteoform ID' and then by 'E-value' in increasing order
    df_sorted = df.sort_values(by=['Proteoform ID', 'E-value'])
    df_sorted = df_sorted.reset_index(drop=True)

    # Initialize a list to keep track of rows to drop
    rows_to_drop = set()

    # Iterate through the DataFrame
    for i in range(len(df_sorted)):
        if i % 1000 == 0:
            print("processing", i)
        if i > 0 and df_sorted.iloc[i]['Proteoform ID'] == df_sorted.iloc[i-1]['Proteoform ID']:
            continue
        for j in range(i + 1, len(df_sorted)):
            if df_sorted.iloc[i]['Proteoform ID'] == df_sorted.iloc[j]['Proteoform ID']:
                if df_sorted.iloc[i]['Protein accession'] != df_sorted.iloc[j]['Protein accession']:
                    rows_to_drop.add(j)
                    print((i,j))
            if df_sorted.iloc[i]['Proteoform ID'] != df_sorted.iloc[j]['Proteoform ID']:
                break

    # Drop the identified rows
    df_filtered = df_sorted.drop(rows_to_drop)
    return df_filtered
