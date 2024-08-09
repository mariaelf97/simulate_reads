import itertools
import pandas as pd
import subprocess

import itertools
import pandas as pd
import subprocess
import random


ISOLATES = [i for i in open("/home/mahmadi/measles_seqs/assemblies/unique_isolates.txt").read().split('\n') if len(i) >0]

# Generate 4 random proportion combinations
isolate_combinations = list(itertools.combinations(ISOLATES, 2))
primer_bed_file = "/home/mahmadi/measles_seqs/primers/primer_v3_400.bed"
for isolate_combination in isolate_combinations:
    proportion1 = round(random.uniform(0, 1), 2)
    proportion2 = round(1.0 - proportion1, 2)
    file1_path = "/home/mahmadi/measles_seqs/assemblies/" + isolate_combination[0] + "/" + isolate_combination[0] + ".fna"
    file2_path = "/home/mahmadi/measles_seqs/assemblies/" + isolate_combination[1] + "/" + isolate_combination[1] + ".fna"
    output_path = "/home/mahmadi/measles_seqs/seq_simulation/combinations/" + isolate_combination[0] + "_" + str(proportion1) + "_" +  isolate_combination[1] + "_" + str(proportion2)
    command = [
            "python", "/home/mahmadi/git_repos/amplicon_sequencing_simulator/workflow/scripts/amplicon_simulator_wrapper.py",
            "-s",f"{isolate_combination[0]},{isolate_combination[1]}" ,
            "-sp", f"{file1_path},{file2_path}",
            "-pr", f"{proportion1},{proportion2}",
            "-p", primer_bed_file,
            "-n", "100000",
            "-o", output_path
        ]
        # Execute the command
    subprocess.run(command, check=True)



