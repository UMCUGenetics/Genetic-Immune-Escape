import pandas as pd
import os

# Config
output_path=config["o"]
donors_file =config["i"]

with open(donors_file) as f:
    df = pd.read_csv(donors_file,sep="\t",names=["sample"])
samples = list(df["sample"].values)

# Rules
rule all:
    input:
        expand(f"{output_path}" + "/{sample}.immune.infiltration.tsv",sample=samples)

rule get_markers:
    input:
        rna_file="/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/isofox/{sample}/{sample}.isf.transcript_data.csv" # this is the path of isofox output, change it to your location
    output:
        output_sample=f"{output_path}" + "/{sample}.immune.infiltration.tsv"
    shell:
        'set +eu '
        ' && PS1=dummy '
        ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && ' \
        'python get_scores_infiltration.py --rna_data {input.rna_file} ' \
        '--column_name AdjTPM --sample_id {wildcards.sample} --output_file {output.output_sample}'
