# required packages
import yaml
import json
import pandas as pd
from urllib import request
import requests
import logging
from datetime import datetime
import os


### borrowing the sort_lineages function from the cov-lineages repo
def sort_lineages(lin_list):
    splitted = [i.split(".") for i in lin_list]
    numeric = []
    for i in splitted:
        lin = [i[0]]
        for j in i[1:]:
            lin.append(int(j))
        numeric.append(lin)
    sorted_list = sorted(numeric)
    stringed = []
    for i in sorted_list:
        lin = [i[0]]
        for j in i[1:]:
            lin.append(str(j))
        stringed.append(lin)
    finished_list = ['.'.join(i) for i in stringed]
    return finished_list


lineages = pd.read_csv('freschi.tsv',sep='\t')['#lineage'].to_list()

lineages = [l.replace('*','') for l in lineages]

lineage_info = {}
for j,lineage in enumerate(lineages):
    children_names = [lin for lin in lineages if lin.startswith(lineage)]
    print(children_names)
    if len(lineage.split('.'))>1:
        splits = lineage.split('.')
        lineage_info[lineage] = {'children':children_names,'parent':'.'.join(splits[0:(len(splits)-1)])}
    else:
        lineage_info[lineage] = {'children':children_names}

# ### now go back and sort the lineage lists
# for lin in lineage_info.keys():
#     lineage_info[lin]['children'] = sort_lineages(lineage_info[lin]['children'])


# also write a copy of the fixed yaml file
# borrowed from cov-lineages/grinch
with open("lineages.yml", "w") as lineage_file:
    for lineage in lineage_info.keys():
        lineage_file.write("- name: " + lineage + "\n")
        # lineage_file.write("  alias: " + lineage_info[lineage]['alias'] + "\n")
        lineage_file.write("  children:\n")

        for child in lineage_info[lineage]['children']:
            lineage_file.write("      - " + child + "\n")
        if 'parent' in lineage_info[lineage].keys():
            if lineage_info[lineage]['parent'] is not None:
                lineage_file.write("  parent: " + lineage_info[lineage]['parent'] + "\n")