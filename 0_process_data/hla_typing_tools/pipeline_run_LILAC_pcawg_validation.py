
import os,shutil
import glob,json
import pandas as pd


# Config
output_path=config["o"]
donors_file =config["i"] # Identifiers of samples to be processed, could be json or txt file

if ".txt" in donors_file:
    with open(donors_file) as f:
        donors = f.read().splitlines()
else:
    with open(donors_file) as f:
        donors=json.load(f)

# Rules
rule all:
    input:
        expand(f"{output_path}" + "{sample}/{sample}_tumor.lilac.csv.gz",sample=donors),
        expand(f"{output_path}" + "{sample}/{sample}_normal.lilac.csv.gz",sample=donors),



rule select_hla_locus:
    input:
        tumor_bam=ancient("/hpc/cuppen/shared_resources/PCAWG/pipeline5/HLA/{sample}T_6:29854528-32726735.bam"),
        normal_bam = ancient("/hpc/cuppen/shared_resources/PCAWG/pipeline5/HLA/{sample}R_6:29854528-32726735.bam"),
    output:
        normal_bam=temp(f"{output_path}"+"{sample}/normal.bam"), # temp
        normal_bai=temp(f"{output_path}"+"{sample}/normal.bam.bai"),
        tumor_bam =temp(f"{output_path}"+"{sample}/tumor.bam"), # temp
        tumor_bai =temp(f"{output_path}"+"{sample}/tumor.bam.bai")

    shell:
        'module load sambamcram/samtools/1.7; mkdir -p {output_path}/{wildcards.sample}; cp {input.normal_bam} {output.normal_bam}; samtools index {output.normal_bam};  cp {input.tumor_bam} {output.tumor_bam}; samtools index {output.tumor_bam} '


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