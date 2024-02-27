#!/bin/bash
# simulate reads using grinder
# -rf reference genome to simulate reads from
# -rs random seed to make sure the result is reproducible
# -cf coverage fold
# -rd  Desired shotgun or amplicon read length distribution specified as:
# average length, distribution ('uniform' or 'normal') and standard deviation
# -fr Use DNA amplicon sequencing using a forward and reverse PCR primer
# sequence provided in a FASTA file. The reference sequences and their
# reverse complement will be searched for PCR primer matches.
# -cb copy bias (use if the reference sequences are whole genome)
# -lb length bias sample reference sequences proportionally to their length
# (use if the reference sequences are whole genome)
# -md to model Illumina errors
# -fq provide output in fastq format
# -ql generate basic quality scores for the simulated reads
# Introduce sequencing errors in the reads under the form of homopolymeric stretches with Balzer model
for FILE in $(cat /home/mahmadi/tb_seqs/primers/primer_names.txt)
do
grinder -rf /home/mahmadi/tb_seqs/seq_simulation/all_seqs.fasta -rs 10 -cf 300 -rd 300 uniform 10 -fr /home/mahmadi/tb_seqs/primers/primer_pairs_fastq/$FILE -cb 1 -lb 0 -md poly4 3e-3 3.3e-8 -fq 1 -ql 30 20 -od /home/mahmadi/tb_seqs/seq_simulation/grinder/$FILE -hd Balzer
done
