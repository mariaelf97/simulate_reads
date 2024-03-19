# Amplicon read simulation for waste water

# Steps to simulate reads

1. Convert a bed-formatted primer file into fastq format.
For this purpose, user can use the following command after cloning the repo``python src/primer_file_conversions/bed_primer_to_fastq.py --p primer.bed --output output.fastq``
2. There are two methods you may use to create amplicons. The first method involves only using a simple string search to find the primers and generate the amplicon. The second method however, uses an alignment approach developed by Boulton et al.
- Alignment approach: ``python src/create_amplicons/alignment_based/make_amplicons_alignment_based.py -g genome.fasta -p primer.fastq -o output_directory --verbose VERBOSE``
- String search approach: This method involves passing pairs of primers in a fastq format to the script and using a snakefile to generate amplicons for all primer pairs. First generate primer pair files using ```python src/primer_file_conversions/split_primer_files.py -i primer.fastq -o output_folder```. From there, you can create a list of primer pair files created using``ls folder > primer_names.txt``. This will be used in the snakefile to generate amplicons for each primer pair. ``cd src/create_amplicons/find_pattern_based/`` and run the snakefile using ``snakemake --cores 10 --use-conda`` to generate amplicons. Please remember to change the file paths in the snakefile accordingly 
3. Merge amplicons
- Amplicons can be merged into a single multi-fasta file using ``cat amplicon_folder/*.fasta >> amplicon_library.fasta``
4. Simulate reads using amplicon library
- reads can be simulated using different read simulators. For instance, WGSIM can simulate illumina reads using ``wgsim -1 400 -2 400 -r 0.001 -R 0.001 -X 0 -e 0.001 -N 10000 amplicon_library.fasta reads1.fastq reads2.fastq``. If you need to generate reads for multiple isolates, you may use the snakefile located in``src/wgsim-amplicon-based``. Please remember to change the file paths in the snakefile accordingly.