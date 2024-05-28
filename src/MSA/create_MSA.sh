#!/bin/bash
#Long assembly to reference mapping
# (-k19 -w10 -U50,500 --rmq -r100k -g10k -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50).
# Up to 20% sequence divergence.
minimap2 -a -x asm20 ~/tb_seqs/H37Rv/H37Rv.fasta all_isolates.fasta -o alinged.sam
gofasta sam toMultiAlign -s alinged.sam -o aligned.fasta