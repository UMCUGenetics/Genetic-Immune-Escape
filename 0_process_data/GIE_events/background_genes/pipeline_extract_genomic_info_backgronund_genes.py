import os,shutil
import glob
import pandas as pd
# how to call it
#snakemake --profile slurm --snakefile pipeline_extract_genomic_info_backgronund_genes.py --config   o=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/hla_events/ i=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/whitelisted_samples.tsv d=hmf --drop-metadata  -np


# Config
dataset =config["d"]
donors_file =config["i"]
output_path=config["o"]

# read samples to be analyzed
df = pd.read_csv(donors_file,sep="\t")
samples = list(df["sample_id"].values)
print ("number of samples", len(samples))

# Rules
rule all:
    input:
        expand(f"{output_path}" + "{sample}/report_background_genes/"+"genomic_alterations.tsv.gz",sample=samples)

rule get_summary:
    input:
        input_data="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/{dataset}/hla_events/{sample}/report_escape/report_status_immune_escape.tsv.gz",
    output:
        output_somatic_summary = "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/{dataset}/hla_events/{sample}/report_background_genes/"+"genomic_alterations.tsv.gz"
    params: cluster_memory = "8G"
    run:
        if dataset == "hmf":
            somatic_vcf= f"/hpc/cuppen/shared_resources/HMF_data/DR-104-update4/somatics/{wildcards.sample}/purple/{wildcards.sample}.purple.somatic.vcf.gz"
            somatic_cnv = f"/hpc/cuppen/shared_resources/HMF_data/DR-104-update4/somatics/{wildcards.sample}/purple/{wildcards.sample}.purple.cnv.gene.tsv"
            sample_info = f"/hpc/cuppen/shared_resources/HMF_data/DR-104-update4/somatics/{wildcards.sample}/purple/{wildcards.sample}.purple.purity.tsv"
        else:
            somatic_vcf= f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/purplesoft3.3/{wildcards.sample}T.purple.somatic.vcf.gz"
            somatic_cnv = f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/purplesoft3.3/{wildcards.sample}T.purple.cnv.gene.tsv"
            sample_info = f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/purplesoft3.3/{wildcards.sample}T.purple.purity.tsv"

        shell('set +eu '
        ' && PS1=dummy '
        ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && ' \
        'python /home/cog/fmartinez/scripts/paper_immuno_scripts/process_data/random_genes/extract_data_random_genes.py --somatic_vcf {somatic_vcf} ' \
        '--somatic_cnv {somatic_cnv} --purity_info {sample_info} --output_file {output.output_somatic_summary}')
        






