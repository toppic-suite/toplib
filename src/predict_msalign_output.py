import re
import sys

def fix_predicted_msalign(input_file, output_file):
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        seq_len = None

        for line in fin:
            line_strip = line.strip()

            # capture sequence length
            if line_strip.startswith("DATABASE_SEQUENCE="):
                seq = line_strip.split("=", 1)[1].strip()
                seq_len = len(seq)
                fout.write(line)
                continue

            # keep metadata
            if not line_strip or not line_strip[0].isdigit():
                fout.write(line)
                continue

            parts = line_strip.split()
            if len(parts) < 4:
                continue

            mass, intensity, charge, ion_full = parts[:4]

            # split ion type and position
            match = re.match(r"([a-zA-Z_]+)(\d+)", ion_full)
            if match:
                ion_type = match.group(1)
                pos = int(match.group(2))
            else:
                ion_type = ion_full
                pos = 0

            # compute N-terminal position
            if ion_type.startswith("b"):# or ion_type.startswith("c"):
                nterm_pos = pos
            elif ion_type.startswith("y"):# or ion_type.startswith("z"):
                if seq_len is not None:
                    nterm_pos = seq_len - pos
                else:
                    nterm_pos = 0
            else:
                nterm_pos = 0

            # dummy confidence
            confidence = 1

            # final line (match real format)
            new_line = f"{mass}\t{intensity}\t{charge}\t{confidence}\t{ion_type}{pos}\t{pos}\t{nterm_pos}\t0\t0\t0"

            fout.write(new_line + "\n")


if __name__ == "__main__":
    msalign_file = sys.argv[1]
    output_file = sys.argv[2]
    fix_predicted_msalign(msalign_file, output_file)
