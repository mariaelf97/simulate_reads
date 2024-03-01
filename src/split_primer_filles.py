import argparse

from Bio import SeqIO


def split_primer_files(fastq_file,output):
    with open(fastq_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fastq"))

        for i in range(0, len(records), 2):
            record1 = records[i]
            primer_name = record1.id
            output_path = output + primer_name.split("_")[1]
            with open(output_path, "w") as output_handle:
                SeqIO.write(record1, output_handle, "fastq")
                if i + 1 < len(records):
                    record2 = records[i + 1]
                    SeqIO.write(record2, output_handle, "fastq")


def main():
    parser = argparse.ArgumentParser(description='split primers into primer pair files')
    parser.add_argument('-i', '--input', help='fastq file containing all primers', required=True)
    parser.add_argument('-o', '--output', help='output path for the location to save primer files', required=True)
    args = parser.parse_args()
    split_primer_files(args.input)


if __name__ == "__main__":
    main()

