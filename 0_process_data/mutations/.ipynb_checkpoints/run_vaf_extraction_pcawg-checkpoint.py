import os,json
import pandas as pd
import gzip

# Config
output_path=config["o"]
donors_file =config["i"]

with open(donors_file) as f:
    donors = f.read().splitlines()

# Rules
rule all:
    input:
        expand(f"{output_path}" + "{sample}/{sample}.vaf_info.tsv.gz", sample=donors),

rule parsing_input:
    input:
        input_file = "/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{sample}-from-jar/purplesoft3.3/{sample}T.purple.somatic.postprocessed.vcf.gz",
    output:
        parsed_muts = f"{output_path}" + "{sample}/{sample}.vaf_info.tsv.gz"
    params: cluster_memory = "24GB"
    shell:
        'set +eu '
        ' && PS1=dummy '
        ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && ' \
        'mkdir -p {output_path}/{wildcards.sample}; ' \
        'python process_mutations_purple.py --file_input {input.input_file} --output_file {output.parsed_muts}'
