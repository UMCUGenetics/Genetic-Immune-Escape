
import os,shutil
import glob

import pandas as pd

# Config
output_path=config["o"]
donors_file =config["i"] # Identifiers of samples to be processed

selected_samples = [ os.path.basename(ff).split(".")[0] for ff in glob.glob("/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/hla/*.hla.bam") if ff.split(".")[0][-1] !="R"]

print ("number of samples", len(selected_samples),selected_samples[0])
with open(donors_file) as f:
    df = pd.read_csv(donors_file,sep="\t")
    df=df[df["tumor_sample"].isin(selected_samples)]
normal_samples = list(df["normal_sample"].values)
samples = list(df["tumor_sample"].values)
sets = list(df["setName"].values)
d_keys = dict(zip(samples,normal_samples))
d_sets = dict(zip(samples,sets))


# Rules
rule all:
    input:
        expand(f"{output_path}" + "{sample}/{sample}_tumor.lilac.csv.gz",sample=samples),
        expand(f"{output_path}" + "{sample}/{sample}_normal.lilac.csv.gz",sample=samples),


rule select_hla_locus:
    input:
        tumor_bam = "/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/hla/{sample}.hla.bam",
    output:
        normal_bam=temp(f"{output_path}"+"{sample}/normal.bam"), # temp
        normal_bai=temp(f"{output_path}"+"{sample}/normal.bam.bai"),
        tumor_bam =temp(f"{output_path}"+"{sample}/tumor.bam"), # temp
        tumor_bai =temp(f"{output_path}"+"{sample}/tumor.bam.bai")
    run:
        normal_id=d_keys[wildcards.sample]
        normal_bam=input.tumor_bam.replace(wildcards.sample,normal_id)
        shell(f'module load sambamcram/samtools/1.7; mkdir -p {output_path}/{wildcards.sample}; cp {normal_bam} {output.normal_bam}; samtools index {output.normal_bam};  cp {input.tumor_bam} {output.tumor_bam}; samtools index {output.tumor_bam} ')

rule run_lilac_tumor:
    input:
        tumor_bam = temp(f"{output_path}" + "{sample}/tumor.bam"),
        tumor_bai =temp(f"{output_path}"+"{sample}/tumor.bam.bai"),
    output:
        output_lilac = f"{output_path}" + "{sample}/{sample}_tumor.lilac.csv"
    threads: 2
    params: cluster_memory = "32G"
    run:
        out = os.path.dirname(input.tumor_bam)
        path_base = os.path.join(output_path,f"{wildcards.sample}", "lilac_output")
        #shell('java -jar /hpc/local/CentOS7/cog/software/lilac/v1/lilac_v1.0_beta.jar -sample {wildcards.sample} -ref_genome /hpc/cuppen/shared_resources/genomes/GRCh37/Sequence/genome.fa -resource_dir /hpc/local/CentOS7/cog/software/lilac/v1/resources/ '
        #      '-reference_bam {input.normal_bam} -tumor_bam {input.tumor_bam} -output_dir {out}')
        sample_name=wildcards.sample+"_tumor"
        shell(f'java -jar /hpc/local/CentOS7/cog/software/lilac/v1/lilac_v1.0_beta.jar -sample {sample_name} -ref_genome /hpc/cuppen/shared_resources/genomes/GRCh37/Sequence/genome.fa -resource_dir /hpc/local/CentOS7/cog/software/lilac/v1/resources/ '
              '-reference_bam {input.tumor_bam} -output_dir {out}')
        #shell("gzip -f {out}/*.csv")
        #shell("gzip -f {out}/*.txt")

rule run_lilac_normal:
    input:
        normal_bam = temp(f"{output_path}" + "{sample}/normal.bam"),
        tumor_bai =temp(f"{output_path}"+"{sample}/normal.bam.bai"),
    output:
        output_lilac = f"{output_path}" + "{sample}/{sample}_normal.lilac.csv"
    threads: 2
    params: cluster_memory = "16G"
    run:
        out = os.path.dirname(input.normal_bam)
        path_base = os.path.join(output_path,f"{wildcards.sample}", "lilac_output")
        sample_name=wildcards.sample+"_normal"
        shell(f'java -jar /hpc/local/CentOS7/cog/software/lilac/v1/lilac_v1.0_beta.jar -sample {sample_name} -ref_genome /hpc/cuppen/shared_resources/genomes/GRCh37/Sequence/genome.fa -resource_dir /hpc/local/CentOS7/cog/software/lilac/v1/resources/ '
              '-reference_bam {input.normal_bam} -output_dir {out}')

rule compress_data:
    input:
        output_lilac_tumor = f"{output_path}" + "{sample}/{sample}_normal.lilac.csv",
        output_lilac_normal = f"{output_path}" + "{sample}/{sample}_tumor.lilac.csv"
    output:
        output_lilac_tumor = f"{output_path}" + "{sample}/{sample}_normal.lilac.csv.gz",
        output_lilac_normal = f"{output_path}" + "{sample}/{sample}_tumor.lilac.csv.gz"
    run:
        out = os.path.dirname(input.output_lilac_tumor)
        shell(f"gzip -f {out}/*.csv")
        shell(f"gzip -f {out}/*.txt")