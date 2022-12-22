import os,json
import pandas as pd
import gzip

# Config
output_path=config["o"]
donors_file =config["i"]

with open(donors_file) as f:
    df = pd.read_csv(donors_file,sep="\t",names=["sample"])

samples = list(df["sample"].values)

for sample in samples: 
    if not(os.path.exists(f"/hpc/cuppen/shared_resources/HMF_data/DR-104-update4/somatics/{sample}/purple/{sample}.purple.somatic.vcf.gz")): # This is the path of the Hartwig purple output for each sample
        samples.remove(sample)
        print (sample, " Removed")
# Rules
rule all:
    input:
        expand(f"{output_path}" + "{sample}/{sample}.vaf_info.tsv.gz", sample=samples),

rule parsing_input:
    input:
        input_file = "/hpc/cuppen/shared_resources/HMF_data/DR-104-update4/somatics/{sample}/purple/{sample}.purple.somatic.vcf.gz",
    output:
        parsed_muts = f"{output_path}" + "{sample}/{sample}.vaf_info.tsv.gz"
    params: cluster_memory = "16GB"
    shell:
        'set +eu '
        ' && PS1=dummy '
        ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && ' \
        'mkdir -p {output_path}/{wildcards.sample}; ' \
        'python process_mutations_purple.py --file_input {input.input_file} --output_file {output.parsed_muts}'
