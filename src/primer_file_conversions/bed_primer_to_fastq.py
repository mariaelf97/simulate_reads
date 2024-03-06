import argparse


def bed_to_fastq(bed_file, output_fastq):
    with open(bed_file, 'r') as bed, open(output_fastq, 'w') as fastq:
        for line in bed:
            parts = line.strip().split('\t')
            primer_name = parts[3]
            primer_name = primer_name.split("_")
            corrected_primer_name = primer_name[0] + "_" + primer_name[1] + "_" + primer_name[2]
            primer_seq = parts[6]
            # Write to fastq file
            fastq.write(f"@{corrected_primer_name}\n{primer_seq}\n+\n{'I' * len(primer_seq)}\n")


def main():
    parser = argparse.ArgumentParser(description='change a bed file to a fastq format file')
    parser.add_argument('-o', '--output', help='Output file path', required=True)
    parser.add_argument('-p', '--primers_bed', help='bed formatted primer file', required=True)
    args = parser.parse_args()

    bed_to_fastq(args.primers_bed, args.output)


if __name__ == "__main__":
    main()
