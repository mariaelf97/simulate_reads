import argparse

from Bio import SeqIO
from Bio.Seq import Seq


def extract_amplicons(genome_file, forward_reverse_primer_file, output_file, primer_name):
    sequences = []
    with open(forward_reverse_primer_file, "r") as fastq_file:
        for record in SeqIO.parse(fastq_file, "fastq"):
            sequence = str(record.seq)
            sequences.append(sequence)
    forward_primer = sequences[0]
    reverse_primer = Seq(sequences[1])
    complementary_reverse = reverse_primer.reverse_complement()
    complementary_reverse = str(complementary_reverse)
    for seq_record in SeqIO.parse(genome_file, "fasta"):
        forward_loc = seq_record.seq.find(forward_primer)
        reverse_loc = seq_record.seq.find(complementary_reverse)
        amplicon = seq_record.seq[forward_loc:reverse_loc + len(complementary_reverse)]
        with open(output_file, "w") as output_handle:
            output_handle.write(f">{primer_name}\n{amplicon}\n")


def main():
    parser = argparse.ArgumentParser(description="Create a FASTA file.")
    parser.add_argument("--genome", "-g", help="whole genome fasta file ", required=True)
    parser.add_argument("--primer_file", "-f", help="forward reverse primer files in fastq format", required=True)
    parser.add_argument("--output", "-o", help="Output file path", required=True)
    parser.add_argument("--primer", "-p", help="primer name to be used when writing amplicon file", required=True)

    args = parser.parse_args()
    extract_amplicons(args.genome, args.primer_file, args.output, args.primer)


if __name__ == "__main__":
    main()