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
output_fasta = '../SNR-u2023011/analysis/snrIII/output.fasta'
ids_file = '../SNR-u2023011/analysis/snrIII/queryManoj.txt'

extract_sequences(input_fasta, output_fasta, ids_file)
