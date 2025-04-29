# load the library for reticulate first

def extract_sequences(input_fasta, output_fasta, ids_file):
    # Read IDs from the text file
    with open(ids_file) as f:
        ids = set(line.strip() for line in f)

    # Initialize variables
    sequences = {}
    current_id = None

    # Read the input FASTA file
    with open(input_fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # If the line starts with ">", it indicates a new sequence
                current_id = line[1:]  # Extract the ID
                if current_id in ids:
                    sequences[current_id] = []
            elif current_id in ids:
                # If we are currently processing a sequence and the ID is in the IDs list, append the line to the sequence
                sequences[current_id].append(line)

    # Write the extracted sequences to a new FASTA file
    with open(output_fasta, 'w') as f:
        for seq_id, seq_data in sequences.items():
            f.write(f'>{seq_id}\n')
            f.write(''.join(seq_data))
            f.write('\n')


# Example usage:
input_fasta = '../SNR-u2023011/analysis/snrIII/Potrx01-protein.fa'
output_fasta = '/mnt/picea/home/schoudhary/shruti/SNRIII/data/output.fasta'
ids_file = '../SNR-u2023011/analysis/snrIII/queryManoj.txt'

extract_sequences(input_fasta, output_fasta, ids_file)

# For Anna
import gzip

import gzip

def extract_sequences(input_fasta, output_fasta, ids_file, missing_ids_file):
    # Read base IDs from the text file
    with open(ids_file, 'r') as f:
        ids = set(line.strip() for line in f)

    # Initialize variables
    sequences = {}
    found_ids = set()
    current_id = None

    # Determine if the input file is gzip compressed
    is_gzipped = input_fasta.endswith('.gz')

    # Open the input FASTA file with appropriate function
    open_func = gzip.open if is_gzipped else open

    # Use 'rb' mode for .gz files to handle binary content and decode as needed
    with open_func(input_fasta, 'rt' if not is_gzipped else 'rb') as f:
        for line in f:
            if is_gzipped:
                line = line.decode('utf-8').strip()  # Decode binary data to text
            else:
                line = line.strip()

            if line.startswith('>'):
                # Extract the base ID by removing any suffix after the first period
                current_id = line[1:].split()[0].split('.')[0]
                if current_id in ids:
                    print(f"Found ID: {current_id}")  # Debugging: Print matching IDs
                    sequences[current_id] = []
                    found_ids.add(current_id)
            elif current_id in ids:
                # If we are currently processing a sequence and the base ID is in the IDs list, append the line to the sequence
                sequences[current_id].append(line)

    # Check for missing IDs
    missing_ids = ids - found_ids
    if missing_ids:
        print(f"Missing IDs: {missing_ids}")
        # Write missing IDs to the missing_ids_file
        with open(missing_ids_file, 'w') as f:
            for missing_id in missing_ids:
                f.write(missing_id + '\n')

    # Open the output FASTA file (supports both regular and .gz files)
    open_func = gzip.open if output_fasta.endswith('.gz') else open

    with open_func(output_fasta, 'wt') as f:  # 'wt' mode for text writing
        for seq_id, seq_data in sequences.items():
            f.write(f'>{seq_id}\n')
            f.write(''.join(seq_data))
            f.write('\n')

    print(f"Extraction complete. Number of sequences written: {len(sequences)}")
    print(f"Number of missing IDs: {len(missing_ids)}")

# File paths
input_fasta = '/mnt/picea/home/schoudhary/shruti/PRN2-FACS/reference/fasta/Potra02_transcripts.fasta.gz'
output_fasta = '/mnt/picea/home/schoudhary/shruti/SNRIII/data/output.fasta'
ids_file ='/mnt/picea/home/schoudhary/shruti/SNRIII/data/queryAnna.txt'
missing_ids_file = '/mnt/picea/home/schoudhary/shruti/SNRIII/data/missing_ids.txt'

# Extract sequences and capture missing IDs
extract_sequences(input_fasta, output_fasta, ids_file, missing_ids_file)
