import pandas as pd
import numpy as np
import sys
import os
import argparse

def load_msalign_data_pred(msalign_file):
    rows = []
    dataset_id = None
    mzml_filename = None
    msalign_filename = None
    precursor_charge = None
    activation = None
    collision_energy = None
    database_sequence = None
    instrument = None

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
                msalign_filename = None
                precursor_charge = None
                activation = None
                collision_energy = None
                database_sequence = None
                instrument = None

                mass_all = []
                intensity_all = []
                charge_all = []
                confidence_all = []

            elif line.startswith("DATASET_ID="):
                dataset_id = line.split("=", 1)[1]

            elif line.startswith("MZML_FILE_NAME="):
                mzml_filename = line.split("=", 1)[1]

            elif line.startswith("MSALIGN_FILE_NAME="):
                msalign_filename = line.split("=", 1)[1]

            elif line.startswith("PRECURSOR_CHARGE="):
                charge_str = line.split("=", 1)[1]
                precursor_charge = int(charge_str.split(":", 1)[0])

            elif line.startswith("DATABASE_SEQUENCE="):
                database_sequence = line.split("=", 1)[1]
                    
            elif line.startswith("COLLISION_ENERGY="):
                collision_energy = float(line.split("=", 1)[1])
                
            elif line.startswith("INSTRUMENT="):
                instrument = line.split("=", 1)[1]
                
            elif line.startswith("ACTIVATION="):
                activation = line.split("=", 1)[1]

            elif line == "END IONS":
                rows.append({
                    "dataset_id": dataset_id,
                    "mzml_file_name": mzml_filename,
                    "msalign_file_name": msalign_filename,
                    "precursor_charge": precursor_charge,
                    "activation": activation,
                    "instrument": instrument,
                    "collision_energy": collision_energy,
                    "database_sequence": database_sequence,
                    "mass_all": mass_all,
                    "intensity_all": intensity_all,
                    "charge_all": charge_all,
                    "confidence_all": confidence_all
                })

                # clean up explicitly
                dataset_id = None
                mzml_filename = None
                precursor_charge = None
                msalign_filename = None
                activation = None
                collision_energy = None
                database_sequence = None
                instrument = None

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


def write_msalign(output_msalign, df):
    with open(output_msalign, "w") as fout:

        for _, row in df.iterrows():
            fout.write("BEGIN IONS\n")
            fout.write(f"DATASET_ID={row['dataset_id']}\n")
            fout.write(f"MZML_FILE_NAME={row['mzml_file_name']}\n")
            fout.write(f"MSALIGN_FILE_NAME={row['msalign_file_name']}\n")
            fout.write(f"PRECURSOR_CHARGE={row['precursor_charge']}\n")
            fout.write(f"ACTIVATION={row['activation']}\n")
            fout.write(f"INSTRUMENT={row['instrument']}\n")
            fout.write(f"COLLISION_ENERGY={row['collision_energy']}\n")
            fout.write(f"DATABASE_SEQUENCE={row['database_sequence']}\n")

            for mass, inten, charge, conf in zip(
                row["mass_all"],
                row["absolute_intensity"],   # use scaled intensity
                row["charge_all"],
                row["confidence_all"]
            ):
                fout.write(
                    f"{mass:.5f}\t{inten:.2f}\t{charge}\t{conf}\n"
                )

            fout.write("END IONS\n\n")

    print(f"Finished writing {len(df):,} spectra.")


def scale_predicted_msalign(input_msalign, metadata_tsv, output_msalign):
    # tsv_file = r"F:\v1.0.3\merged_info_add_project_id_v2\combined\msv000091610\toplib\rebuilt\td_pred\lib_without_mass_shift_info.tsv"
    tsv_df = pd.read_csv(metadata_tsv, sep='\t')
    # msalign_file = r"F:\v1.0.3\merged_info_add_project_id_v2\combined\msv000091610\toplib\rebuilt\td_pred\OneDrive_1_8-4-2026\lib_without_mass_shift_meta_predict_v0_convert.msalign"
    td_ms2_df = load_msalign_data_pred(input_msalign)

    td_abs_inte_all = []
    for i in range(len(td_ms2_df)):
        pre_inte = tsv_df['precursor_intensity'].iloc[i]
        td_abs_inte = pre_inte * np.array(td_ms2_df['intensity_all'].iloc[i])
        td_abs_inte_all.append(td_abs_inte)

    td_ms2_df['absolute_intensity'] = td_abs_inte_all

    # write it to a msalign file second column in absolute_intensity
    write_msalign(output_msalign, td_ms2_df)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Scale TD-Pred relative fragment intensities using precursor intensities from a TSV file."
        )
    )

    parser.add_argument(
        "--msalign",
        required=True,
        help="Input TD-Pred predicted MSAlign file."
    )
    parser.add_argument(
        "--tsv",
        required=True,
        help="Input TSV containing precursor_intensity."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output MSAlign file with scaled fragment intensities."
    )

    args = parser.parse_args()

    if not os.path.isfile(args.msalign):
        print(f"Error: MSAlign file not found: {args.msalign}")
        sys.exit(1)

    if not os.path.isfile(args.tsv):
        print(f"Error: TSV file not found: {args.tsv}")
        sys.exit(1)

    scale_predicted_msalign(
        input_msalign=args.msalign,
        metadata_tsv=args.tsv,
        output_msalign=args.output
    )


if __name__ == "__main__":
    main()