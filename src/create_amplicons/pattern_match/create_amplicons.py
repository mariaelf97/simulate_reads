import math
import os
from os.path import basename
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from regex import regex
import argparse


def find_closest_match(pattern,reference_seq):
    """function to find a string allowing up to 3 mismatches"""
    for m in regex.finditer(r"\L<primer_string>{s<=2,i<1,d<1}",
                                reference_seq, primer_string=[pattern]):
        return m.start()


def make_amplicon(left_primer_loc,right_primer_loc,
                  primer_seq_y, reference):
    """function to create amplicons based on the string match location"""
    if math.isnan(left_primer_loc) or math.isnan(right_primer_loc):
        amplicon = ""
        return amplicon
    else:
        amplicon = reference[int(left_primer_loc - 1): int(right_primer_loc + len(primer_seq_y))]
        return amplicon

def main():
    parser = argparse.ArgumentParser(description="Create amplicons for a genome using a primer set.")
    parser.add_argument("--genome_path", "-g", help="Path to the genome of interest.")
    parser.add_argument("--output", "-o", help="Folder where the output will go")
    parser.add_argument("--primers_file", "-p", help="Path to primer bed file")

    args = parser.parse_args()
    genome_filename_short = ".".join(basename(args.genome_path).split(".")[:-1])
    reference = next(SeqIO.parse(args.genome_path, "fasta"))
    reference_seq = str(reference.seq)
    col_names = ["ref","start", "end", "left_right", "chr","strand", "primer_seq"]
    primer_bed = pd.read_csv(args.primers_file, sep= "\t", names=col_names)
    primer_bed["amplicon_number"] = primer_bed["left_right"].str.split('_').str[1]
    df = pd.merge(
        primer_bed.loc[primer_bed["left_right"].str.contains("LEFT")],
        primer_bed.loc[primer_bed["left_right"].str.contains("RIGHT")],
        on=["amplicon_number"]
    )
    # select needed columns
    merged_df = df[["amplicon_number","primer_seq_x","primer_seq_y"]]
    merged_df["comp_rev"] = merged_df.apply(lambda row: Seq(row['primer_seq_y']).reverse_complement(), axis=1)
    merged_df["left_primer_loc"] = merged_df.apply(lambda row: reference_seq.find(row["primer_seq_x"]), axis=1)
    merged_df["right_primer_loc"] = merged_df.apply(lambda row: reference_seq.find(str(row["comp_rev"])), axis=1)
    merged_df["amplicon_length"] = merged_df.apply(lambda row: row["right_primer_loc"] - row["left_primer_loc"], axis=1)
    failed_amplicons = pd.concat([merged_df[merged_df['amplicon_length']>2000],
                                 merged_df[merged_df['amplicon_length'] == 0],
                                  merged_df[merged_df['amplicon_length'] < 0]])
    failed_amplicons["alt_left_primer_loc"] = failed_amplicons.apply(lambda row: find_closest_match(row["primer_seq_x"],reference_seq), axis=1)
    failed_amplicons["alt_right_primer_loc"] = failed_amplicons.apply(lambda row: find_closest_match(str(row["comp_rev"]),reference_seq,), axis=1)
    merged_df = pd.merge(
        merged_df,
        failed_amplicons[["amplicon_number","alt_left_primer_loc","alt_right_primer_loc"]],
        on=["amplicon_number"],
        how="outer"
    )
    merged_df['left_primer_loc'] = merged_df.apply(
        lambda row: row['alt_left_primer_loc'] if row['left_primer_loc'] == -1 else row['left_primer_loc'],
        axis=1)
    merged_df['right_primer_loc'] = merged_df.apply(
        lambda row: row['alt_right_primer_loc'] if row['right_primer_loc'] == -1 else row['right_primer_loc'],
        axis=1)

    merged_df["amplicon_length"] = merged_df["right_primer_loc"] - merged_df["left_primer_loc"] + len(merged_df["primer_seq_y"])
    merged_df["amplicon_sequence"] = merged_df.apply(lambda row: make_amplicon(row["left_primer_loc"],
                                                                          row["right_primer_loc"],
                                                                          row["primer_seq_y"],reference_seq), axis=1)
    for row in merged_df.itertuples():
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        with open(f"{args.output}/amplicon_{row.amplicon_number}" + ".fasta", "w") as f:
            f.write(f">{reference.id}_amplicon_{row.amplicon_number}" + "\n")
            f.write(row.amplicon_sequence + "\n\n")


if __name__ == "__main__":
    main()