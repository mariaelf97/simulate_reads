import itertools
import pandas as pd

ISOLATES = [i for i in open("/home/maryam/mnt/tb_seqs/seq_simulation/amplicon_read_simulations/amplicon_alignment_wgsim/primer_v2/mixed_samples_diff_cov/isolates.txt").read().split('\n') if len(i) >0]
ABUNDANCES = [i for i in open("/home/maryam/mnt/tb_seqs/seq_simulation/amplicon_read_simulations/amplicon_alignment_wgsim/primer_v2/mixed_samples_diff_cov/abundances.txt").read().split('\n') if len(i) >0]

isolate_combinations = list(itertools.combinations(ISOLATES, 2))
abundance_combinations = [("25000","75000"),("250000","750000")]

for isolate_combination in isolate_combinations:
    for abundance_combination in abundance_combinations:
        file1_path = "/home/maryam/mnt/tb_seqs/seq_simulation/amplicon_read_simulations/amplicon_alignment_wgsim/primer_v2/" + isolate_combination[0] + "/" + abundance_combination[0] + "/reads_1.fastq"
        file2_path = "/home/maryam/mnt/tb_seqs/seq_simulation/amplicon_read_simulations/amplicon_alignment_wgsim/primer_v2/" + isolate_combination[1] + "/" + abundance_combination[1] + "/reads_1.fastq"
        output_path = "/home/maryam/mnt/tb_seqs/seq_simulation/amplicon_read_simulations/amplicon_alignment_wgsim/primer_v2/mixed_samples_diff_cov/" + isolate_combination[0] + "_" + abundance_combination[0] + "_" + isolate_combination[1] + "_" + abundance_combination[1] + ".fastq"
        with open(file1_path, 'r') as file1:
            lines_file1 = file1.readlines()
        with open(file2_path, 'r') as file2:
            lines_file2 = file2.readlines()
        combined_lines = lines_file1 + lines_file2
        with open(output_path, 'w') as output_file:
            output_file.writelines(combined_lines)























