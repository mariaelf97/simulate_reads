import os
import matplotlib.pyplot as plt

ampDir = 'test/amplicons/'
lengths0 = {}
for file in os.listdir(ampDir):    
    infile = open(os.path.join(ampDir,file))
    k=0
    for line in infile:
        if line[0] == '>' or (k>=2):
            pass
        else:
            print(file, len(line))
            lengths0[file] = len(line)
        k+=1
    infile.close()

# read in H37Rv reference, make sure alignments make sense.  
infile = open('ref/H37Rv.fasta')
refseq = ''
k=0
for line in infile:
    if line[0] == '>' or (k>=2 and len(line)<=1):
        pass
    else:
        refseq+=line.strip()
    k+=1
infile.close()

import pandas as pd 
from Bio.Seq import Seq
## read in bed file, make sure that is producing what we expect. 
bed_df = pd.read_csv('primers/primer_V1.bed',header=None,sep='\t')
### check to make sure that the sequence is the same in the reference.
for i in range(0,bed_df.shape[0]):
    start = bed_df.iloc[i,1]
    end = bed_df.iloc[i,2]
    if 'LEFT' in bed_df.iloc[i,3]:
        if refseq[start:end] != bed_df.iloc[i,6]:
            print('error with amplicon')
            print(refseq[start:end])
            print(bed_df.iloc[i,6])
    else:
        seq0 = Seq(bed_df.iloc[i,6]).reverse_complement()
        if refseq[start:end] != seq0:
            print('error with amplicon, reverse')
            print(refseq[start:end])
            print(seq0)

# things seem to make sense, so why is amplicon 68 so long?
amp68_L = bed_df[bed_df.iloc[:,3].str.contains('68') & bed_df.iloc[:,3].str.contains('LEFT')]
amp68_R = bed_df[bed_df.iloc[:,3].str.contains('68') & bed_df.iloc[:,3].str.contains('RIGHT')]

import re
primerL = amp68_L.iloc[0,6]
primerR = str(Seq(amp68_R.iloc[0,6]).reverse_complement())
#check for repeats 
print('68 Left',[m.start() for m in re.finditer(primerL,refseq)])
print('68 Right',[m.start() for m in re.finditer(primerR,refseq)])


#how does amp 44 look? 
amp44_L = bed_df[bed_df.iloc[:,3].str.contains('44') & bed_df.iloc[:,3].str.contains('LEFT')]
amp44_R = bed_df[bed_df.iloc[:,3].str.contains('44') & bed_df.iloc[:,3].str.contains('RIGHT')]

primerL = amp44_L.iloc[0,6]
primerR = str(Seq(amp44_R.iloc[0,6]).reverse_complement())
#check for repeats 
print('44 Left',[m.start() for m in re.finditer(primerL,refseq)])
print('44 Right',[m.start() for m in re.finditer(primerR,refseq)])
