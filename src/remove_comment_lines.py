import sys

def remove_first_29_lines(input_file_path, output_file_path):
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    # Remove the first 29 lines
    lines = lines[29:]

    with open(output_file_path, 'w') as file:
        file.writelines(lines)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file_path> <output_file_path>")
    else:
        remove_first_29_lines(sys.argv[1], sys.argv[2])
        
    