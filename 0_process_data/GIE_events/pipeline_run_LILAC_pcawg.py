
import os,shutil
import glob

import pandas as pd


# Config
output_path=config["o"]
donors_file =config["i"]

with open(donors_file) as f:
    donors = f.read().splitlines()


# Rules
rule all:
    input:
        expand(f"{output_path}" + "{sample}/{sample}.lilac.csv.gz",sample=donors),




rule select_hla_locus:
    input:
        tumor_bam=ancient("/hpc/cuppen/shared_resources/PCAWG/pipeline5/HLA/{sample}T_6:29854528-32726735.bam"), # path with the sliced tumor PCAWG .bams
        normal_bam = ancient("/hpc/cuppen/shared_resources/PCAWG/pipeline5/HLA/{sample}R_6:29854528-32726735.bam"), # path with the sliced tumor PCAWG .bams
    output:
        normal_bam=temp(f"{output_path}"+"{sample}/normal.bam"), # temp
        normal_bai=temp(f"{output_path}"+"{sample}/normal.bam.bai"),
        tumor_bam =temp(f"{output_path}"+"{sample}/tumor.bam"), # temp
        tumor_bai =temp(f"{output_path}"+"{sample}/tumor.bam.bai")

    shell:
        'module load sambamcram/samtools/1.7; mkdir -p {output_path}/{wildcards.sample}; cp {input.normal_bam} {output.normal_bam}; samtools index {output.normal_bam};  cp {input.tumor_bam} {output.tumor_bam}; samtools index {output.tumor_bam} '

rule run_lilac:
    input:
        tumor_bam = temp(f"{output_path}" + "{sample}/tumor.bam"),
        tumor_bai =temp(f"{output_path}"+"{sample}/tumor.bam.bai"),
        normal_bam = temp(f"{output_path}" + "{sample}/normal.bam"),
        normal_bai=temp(f"{output_path}"+"{sample}/normal.bam.bai"),
    output:
        output_lilac = f"{output_path}" + "{sample}/{sample}.lilac.csv.gz"
    threads: 2
    params: cluster_memory = "32G"
    run:
        out = os.path.dirname(input.tumor_bam)
        somatic_vcf= f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/purplesoft3.3/{wildcards.sample}T.purple.somatic.vcf.gz" # path to purple output
        somatic_cnv = f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/purplesoft3.3/{wildcards.sample}T.purple.cnv.gene.tsv" # path to purple output
        sample_info = f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/purplesoft3.3/{wildcards.sample}T.purple.purity.tsv" # path to purple output
        path_base = os.path.join(output_path,f"{wildcards.sample}", "lilac_output")
        shell('java -jar /hpc/local/CentOS7/cog/software/lilac/v1/lilac_v1.0.jar -sample {wildcards.sample} -ref_genome /hpc/cuppen/shared_resources/genomes/GRCh37/Sequence/genome.fa -resource_dir /hpc/local/CentOS7/cog/software/lilac/v1/resources/ '
              '-reference_bam {input.normal_bam} -tumor_bam {input.tumor_bam} -somatic_variants_file {somatic_vcf} -gene_copy_number_file {somatic_cnv} -output_dir {out}')
        shell("gzip -f {out}/*.csv")
        shell("gzip -f {out}/*.txt")
