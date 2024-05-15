import os
import argparse
from Bio import SeqIO
import pandas as pd

def calculate_sequence_lengths(fasta_folder):
    all_sequences = []
    for file_name in os.listdir(fasta_folder):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(fasta_folder, file_name)
            sequences = calculate_sequence_length(file_path, file_name)
            all_sequences.extend(sequences)
    return all_sequences

def calculate_sequence_length(fasta_file, file_name):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_id = file_name
        sequence_length = len(record.seq)
        sequences.append((sequence_id, sequence_length))
    return sequences

def create_dataframe(sequences):
    df = pd.DataFrame(sequences, columns=['sequence_id', 'sequence_length'])
    return df

def main():
    parser = argparse.ArgumentParser(description="Calculate sequence lengths from FASTA files in a folder and store them in a DataFrame.")
    parser.add_argument("input_folder", help="Input folder containing FASTA files")
    parser.add_argument("output_file", help="Output CSV file path")
    args = parser.parse_args()

    input_folder = args.input_folder
    output_file = args.output_file

    sequences = calculate_sequence_lengths(input_folder)
    df = create_dataframe(sequences)
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    main()
