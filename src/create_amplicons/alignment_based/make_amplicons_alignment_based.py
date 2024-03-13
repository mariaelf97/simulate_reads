## Code credit to SWAMPy : https://github.com/goldman-gp-ebi/SWAMPy
import subprocess
import os
from os.path import join, basename
from io import StringIO
from Bio import SeqIO
import pandas as pd
import logging


def build_index(genome_path, output):
    # create bowtie2 index for the reference
    directory="indices"
    path = join(output, directory)
    if not os.path.exists(path):
        os.makedirs(path)
    subprocess.run([
        "bowtie2-build", genome_path, path, "-o",path], capture_output=True)


def align_primers(genome_filename_short,output, primers_files, verbose):
    # run bowtie2 aligner
    directory="indices"
    path = os.path.join(output, directory)
    alignment = subprocess.run(
        ["bowtie2",
         "-x", path,
         "-U", primers_files], capture_output=True)

    alignment = StringIO(alignment.stdout.decode("UTF-8"))

    # read alignment data as a dataframe
    df = pd.read_csv(alignment, sep="\t", skiprows=[0, 1, 2], header=None, names=[i for i in range(20)])
    df = pd.DataFrame(df[[0, 2, 3, 9]])
    df = df.rename(columns={0: "name", 2: "ref", 3: "start", 9: "seq"})

    # split the column "name" to extract useful data
    df["seq_len"] = df["seq"].apply(len)
    df["amplicon_number"] = df["name"].apply(lambda x: int(x.split("_")[1]))
    df["handedness"] = df["name"].apply(lambda x: x.split("_")[2])
    df["is_alt"] = df["name"].apply(lambda x: len(x.split("_")) > 3)

    # remove rows where the alignment mismatched
    drop_rows = []
    for r in df.itertuples():
        if r.ref == "*":
            if verbose:
                logging.info(f"Dropping amplicon {r.amplicon_number}, couldn't find a match for the primer {r.seq}")

            drop_rows.append(r.Index)

    df.drop(drop_rows, inplace=True)

    # inner join the dataframe with itself, to get the pairs of primers and their start/end positions
    df = pd.merge(
        df.loc[df["handedness"] == "LEFT"],
        df.loc[df["handedness"] == "RIGHT"],
        on=["amplicon_number", "is_alt"]
    )

    # rename the columns to more understandable names
    df = pd.DataFrame(
        df[["ref_x", "amplicon_number", "is_alt",
            "start_x", "seq_x", "seq_len_x",
            "start_y", "seq_y", "seq_len_y"]]
    )

    df = df.rename(columns={
        "ref_x": "ref",
        "start_x": "left",
        "seq_x": "left_primer",
        "seq_len_x": "left_primer_length",
        "start_y": "right",
        "seq_y": "right_primer",
        "seq_len_y": "right_primer_length"})

    df["amplicon_filepath"] = genome_filename_short + "_amplicon_" + df["amplicon_number"].map(str) + df["is_alt"].map(
        lambda x: "_alt" if x else "") + ".fasta"

    if verbose:
        logging.info("First 5 rows: ")
        logging.info(df.head())

    return df


def write_amplicon(df, reference, genome_filename_short, output, verbose=False):
    for r in df.itertuples():
        amplicon_number = r.amplicon_number
        alt = r.is_alt
        reference_string = str(reference.seq)
        amplicon = reference_string[r.left - 1: r.right + r.right_primer_length - 1]

        if verbose:
            logging.info(f">{reference.id}_amplicon_{amplicon_number}" + ("_alt" if alt else ""))
            logging.info("length: " + str(r.right - r.left))
            logging.info(amplicon + "\n")
            logging.info(r.left_primer + "-" * (r.right - r.left - r.left_primer_length) + r.right_primer + "\n")
        directory = "amplicons"
        path = join(output, directory)
        if not os.path.exists(path):
            os.makedirs(path)
        with open(f"{path}/{genome_filename_short}_amplicon_{amplicon_number}" + (
        "_alt" if alt else "") + ".fasta", "w") as f:

            f.write(f">{reference.id}_amplicon_{amplicon_number}" + ("_alt" if alt else "") + "\n")
            f.write(amplicon + "\n\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create amplicons for a genome using a primer set.")
    parser.add_argument("--genome_path", "-g", help="Path to the genome of interest.")
    parser.add_argument("--output", "-o", help="Folder where the output will go")
    parser.add_argument("--primers_file", "-p", help="Path to fastq file of primers. Default ARTIC V1 primers.")
    parser.add_argument("--verbose", help="Verbose mode.")

    args = parser.parse_args()
    genome_filename_short = ".".join(basename(args.genome_path).split(".")[:-1])
    reference = SeqIO.read(args.genome_path, format="fasta")

    build_index(args.genome_path, args.output)
    df = align_primers(genome_filename_short, args.output, args.primers_file, verbose=args.verbose)
    write_amplicon(df, reference, genome_filename_short, args.output, verbose=args.verbose)