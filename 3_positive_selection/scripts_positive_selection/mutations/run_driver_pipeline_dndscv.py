
import os,json
import pandas as pd
import numpy as np


# Config
output_path=config["o"]
input_path=config["i"]

with open(input_path,'r') as f:
    d_info=json.load(f) # this is the .json file created in the "create_groups_for_positive_selection.ipynb"

ttypes=list(d_info.keys())
for ttype in ttypes:
    # create list of patients
    file_input = f"{output_path}" + f"/dndscv/{ttype}.json"
    with open(file_input, 'w') as f:
        json.dump(d_info[ttype], f)
# Rules
rule all:
    input:
        expand(f"{output_path}"+"/dndscv/{ttype}.dndscv.results.tsv.gz",ttype=ttypes),

rule create_dataset_input_mutations: # create input datasets
    input:
        input_path=f"{output_path}" + "/dndscv/{ttype}.json"
    output:
        input_dndscv=f"{output_path}"+"/dndscv/{ttype}.dndscv.input.tsv.gz" # temp
    params: cluster_memory = "32G"
    run:
        shell('set +eu '
        ' && PS1=dummy '
        ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && python get_variants_samples.py --patients_file {input.input_path} --ttype_name {wildcards.ttype} --output_file {output.input_dndscv}')

rule run_dndscv:
    input:
        input_dndscv=f"{output_path}"+"/dndscv/{ttype}.dndscv.input.tsv.gz",
    output:
        ouput_dndscv=f"{output_path}"+"/dndscv/{ttype}.dndscv.results.tsv.gz",
    params: cluster_memory = "128G"
    run:
        basepath=os.path.dirname(output.ouput_dndscv)
        shell('set +eu '
            ' && PS1=dummy '
            ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate dndscv && Rscript run_dndscv.R {input.input_dndscv} {basepath} {wildcards.ttype}')

