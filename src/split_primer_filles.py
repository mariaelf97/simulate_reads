from Bio import SeqIO


def main():
    fastq_file = "/home/maryam/mnt/tb_seqs/primers/primer.fastq"
    with open(fastq_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fastq"))

        for i in range(0, len(records), 2):
            record1 = records[i]
            primer_name = record1.id
            output = "/home/maryam/mnt/tb_seqs/primers" + primer_name.split("_")[1]
            with open(output, "w") as output_handle:
                SeqIO.write(record1, output_handle, "fastq")
                if i + 1 < len(records):
                    record2 = records[i + 1]
                    SeqIO.write(record2, output_handle, "fastq")


if __name__ == "__main__":
    main()

