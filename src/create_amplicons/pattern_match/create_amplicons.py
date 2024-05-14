import os
from os.path import basename
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from regex import regex


def find_closest_match(pattern):
    for m in regex.finditer(r"\L<primer_string>{s<=2}",
                            reference_seq, primer_string=[pattern]):
        return m.start()


def make_amplicon(left_primer_loc,right_primer_loc, alt_left_primer_loc, alt_right_primer_loc,primer_seq_y, reference):
    if left_primer_loc == "-1" or right_primer_loc == "-1":
        amplicon = reference[alt_left_primer_loc - 1: alt_right_primer_loc + len(primer_seq_y)]
    else:
        amplicon = reference[left_primer_loc - 1: right_primer_loc + len(primer_seq_y)]
    return len(amplicon)


# def write_amplicon(df, reference, genome_filename_short, output):
#     reference_seq = str(reference.seq)
#     for r in df.itertuples():
#         amplicon_number = r.amplicon_number
#         if r.left_primer_loc == "-1" or r.right_primer_loc == "-1":
#             amplicon = reference_seq[r.alt_left_primer_loc - 1: r.alt_right_primer_loc + len(r.primer_seq_y)]
#         else:
#             amplicon = reference_seq[r.left_primer_loc - 1: r.right_primer_loc + len(r.primer_seq_y)]
#         if not os.path.exists(output):
#             os.makedirs(output)
#         with open(f"{output}/{genome_filename_short}_amplicon_{amplicon_number}" + ".fasta", "w") as f:
#             f.write(f">{reference.id}_amplicon_{amplicon_number}" + "\n")
#             f.write(amplicon + "\n\n")


if __name__ == "__main__":
    import argparse

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
    failed_amplicons["alt_left_primer_loc"] = failed_amplicons.apply(lambda row: find_closest_match(row["primer_seq_x"]), axis=1)
    failed_amplicons["alt_right_primer_loc"] = failed_amplicons.apply(lambda row: find_closest_match(str(row["comp_rev"])), axis=1)
    merged_df = pd.merge(
        merged_df,
        failed_amplicons[["amplicon_number","alt_left_primer_loc","alt_right_primer_loc"]],
        on=["amplicon_number"],
        how="outer"
    )
    merged_df["amplicon_length"] = merged_df.apply(lambda row: make_amplicon(row["left_primer_loc"],
                                                                          row["right_primer_loc"],
                                                                          row["alt_left_primer_loc"],
                                                                          row["alt_right_primer_loc"],
                                                                          row["primer_seq_y"],
                                                                          reference_seq), axis=1)
    merged_df.to_csv(args.output, index=False)
    #write_amplicon(merged_df, reference, genome_filename_short, args.output)