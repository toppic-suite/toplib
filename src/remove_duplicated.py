import pandas as pd


def process_tsv(input_file_path):
    # Load the TSV file into a DataFrame
    df = pd.read_csv(input_file_path, sep='\t')
    # Sort the DataFrame by 'Protein accession' and then by 'E-value' in increasing order
    df_sorted = df.sort_values(by=['Protein accession', 'E-value'])
    df_sorted = df_sorted.reset_index(drop=True)
    # Initialize a list to keep track of rows to drop
    rows_to_drop = set()

    # Iterate through the DataFrame
    for i in range(len(df_sorted)):
        if i % 100 == 0:
            print(i)
        if i in rows_to_drop:
            continue
        for j in range(i + 1, len(df_sorted)):
            if df_sorted.iloc[i]['Protein accession'] == df_sorted.iloc[j]['Protein accession']:
                if abs(float(df_sorted.iloc[i]['Precursor mass']) - float(df_sorted.iloc[j]['Precursor mass'])) < 2.2:
                    rows_to_drop.add(j)
            if df_sorted.iloc[i]['Protein accession'] != df_sorted.iloc[j]['Protein accession']:
                break

    # Drop the identified rows
    df_filtered = df_sorted.drop(rows_to_drop)
    return df_filtered

