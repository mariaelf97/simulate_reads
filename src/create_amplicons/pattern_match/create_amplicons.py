import itertools
from itertools import permutations
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
    matches = [m.start() for m in regex.finditer(r"\L<primer_string>{s<1}",
                                reference_seq, primer_string=[pattern])]
    return matches


def make_amplicon(left_primer_loc,right_primer_loc,
                  primer_seq_y, reference):
    """function to create amplicons based on the string match location"""
    if math.isnan(left_primer_loc) or math.isnan(right_primer_loc):
        amplicon = ""
        return amplicon
    else:
        amplicon = reference[int(left_primer_loc - 1): int(right_primer_loc + len(primer_seq_y))]
        return amplicon


def evaluate_matches(left_primer_coordinates, right_primer_coordinates):
    """function to evaluate which coordinates found for each primer makes a valid amplicon"""
    if len(left_primer_coordinates) !=0 and len(right_primer_coordinates)!=0:
        valid_combinations = []
        combinations = list(itertools.product(left_primer_coordinates, right_primer_coordinates))
        for combination in combinations:
            amplicon_length = combination[1] - combination[0]
            if 0 < amplicon_length <= 2000:
                valid_combinations.append(combination)
            else:
                pass
        return valid_combinations
    else:
        return []


def main():
    parser = argparse.ArgumentParser(description="Create amplicons for a genome using a primer set.")
    parser.add_argument("--genome_path", "-g", help="Path to the genome of interest.")
    parser.add_argument("--output", "-o", help="Folder where the output will go")
    parser.add_argument("--primers_file", "-p", help="Path to primer bed file")

    args = parser.parse_args()
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
    merged_df["left_primer_loc"] = merged_df.apply(lambda row: find_closest_match(row["primer_seq_x"],reference_seq), axis=1)
    merged_df["right_primer_loc"] = merged_df.apply(lambda row: find_closest_match(str(row["comp_rev"]),reference_seq), axis=1)
    merged_df["valid_combinations"] = ""
    d = pd.DataFrame()
    for i in range(len(merged_df)):
        merged_df["valid_combinations"][i] = evaluate_matches(merged_df["left_primer_loc"][i], merged_df["right_primer_loc"][i])
        for j in range(len(merged_df["valid_combinations"][i])):
            primer_start = merged_df["valid_combinations"][i][j][0]
            primer_end = merged_df["valid_combinations"][i][j][1]
            amplicon_number = merged_df["amplicon_number"][i]
            amplicon_length = primer_end - primer_start
            temp = pd.DataFrame(
                {
                    "amplicon_number": [amplicon_number],
                    'primer_start': [primer_start],
                    'primer_end': [primer_end],
                    'amplicon_length': [amplicon_length]
                }
            )

            d = pd.concat([d, temp])
    all_amplicons = pd.merge(d, merged_df[["amplicon_number","primer_seq_y"]], how='outer', sort=False, on='amplicon_number')
    all_amplicons = all_amplicons.fillna(0)


    #all_amplicons.to_csv(args.output, index=False)
    all_amplicons["amplicon_sequence"] = all_amplicons.apply(lambda row: make_amplicon(row["primer_start"],
                                                                          row["primer_end"],
                                                                          row["primer_seq_y"],reference_seq), axis=1)
    for row in all_amplicons.itertuples():
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        with open(f"{args.output}/amplicon_{row.amplicon_number}_" + str(row.Index) + ".fasta", "w") as f:
            f.write(f">{reference.id}_amplicon_{row.amplicon_number}" + "\n")
            f.write(row.amplicon_sequence + "\n\n")


if __name__ == "__main__":
    main()