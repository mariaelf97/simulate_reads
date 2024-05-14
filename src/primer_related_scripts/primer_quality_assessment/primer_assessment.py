import os

import pandas as pd
from Bio import SeqIO

primer_versions = ["V2"]
amplicon_names = pd.read_csv("/home/mahmadi/tb_seqs/primers/primer_names.txt", names= ["amplicon_name"])
isolates = pd.read_csv("/home/mahmadi/tb_seqs/assemblies/long_read_assemblies/pacbio-RSII/isolates.txt", names= ["isolate_name"])

df = pd.DataFrame(columns=["primer_version","amplicon_name","amplicon_length"])
for isolate in isolates["isolate_name"]:
    for primer in primer_versions:
        for i in range(0,164):
            path = ("/home/mahmadi/tb_seqs/seq_simulation/amplicons/hybrid/primer_V1/GCA_014900175.1/amplicons/" +
                    isolate + "/amplicons/" + isolate + "_amplicon_" + str(i) + ".fasta")
            if os.path.exists(path):
                with open(path, "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        new_row = {"isolate" : [isolate], "primer_version": [primer],"amplicon_name": "amplicon_" +  str(i),"amplicon_length": [len(record.seq)]}
                        df = pd.concat([df, pd.DataFrame(new_row)], ignore_index=True)
            else:
                new_row = {"isolate": [isolate], "primer_version": [primer], "amplicon_name": "amplicon_" + str(i),
                           "amplicon_length": "NA"}
                df = pd.concat([df, pd.DataFrame(new_row)], ignore_index=True)

df.to_csv("/home/mahmadi/tb_seqs/primers/primer_stats_hybrid.csv", index= False)
