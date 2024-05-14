# Amplicon read simulation for waste water

# Steps to simulate reads

There are two methods for amplicon extraction using bed-formatted primer file.
1. Alignment approach:
This method use bowtie2 as an aligner to map the primers to
a provided whole genome fasta file. The default parameters for
an end to end alignment are used. 
`src/create_amplicons/alignment_based`

2. Basic string match approach:
This approach uses a strict pattern matching approach allowing up to 3 mismatches. 
This method does not take indels into account.
`src/create_amplicons/pattern_match`

Once the amplicons are created, wgsim is used as a short-read 
simulator to simulate amplicon reads. 
`src/variant_calling`
`src/read_simulation`
