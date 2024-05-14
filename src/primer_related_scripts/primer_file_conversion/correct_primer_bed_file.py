import argparse

import pandas as pd


def correct_bed_file(primer_file, output):
    column_names = ['chrom', 'start', 'end', 'primer_name', "na", "strand","primer_seq"]
    bed_primer = pd.read_csv(primer_file, sep= "\t",
                             names=column_names, header=None)
    bed_primer['primer_name'] = bed_primer['primer_name'].str.split('_').str[:-1].str.join('_')
    bed_primer["chrom"] = 1
    bed_primer.to_csv(output, header=None, index=False, sep="\t")


def main():
    parser = argparse.ArgumentParser(description='fix bed file primer name field')
    parser.add_argument('-o', '--output', help='Output file path', required=True)
    parser.add_argument('-p', '--primers_bed', help='bed formatted primer file', required=True)
    args = parser.parse_args()

    correct_bed_file(args.primers_bed, args.output)


if __name__ == "__main__":
    main()
