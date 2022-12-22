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
        expand(f"{output_path}" + "{donor}/typing_tumor/winners.hla.txt",donor=donors),
        expand(f"{output_path}" + "{donor}/typing_normal/winners.hla.txt",donor=donors),
        expand(f"{output_path}" + "{donor}/typing_normal_xhla/report-{donor}-hla.json",donor=donors),
        expand(f"{output_path}" + "{donor}/typing_tumor_xhla/report-{donor}-hla.json", donor=donors),
        


rule select_hla_locus:
    input:
        tumor_bam=ancient("/hpc/cuppen/shared_resources/PCAWG/pipeline5/HLA/{sample}T_6:29854528-32726735.bam"), # path with the sliced tumor PCAWG .bams
        normal_bam = ancient("/hpc/cuppen/shared_resources/PCAWG/pipeline5/HLA/{sample}R_6:29854528-32726735.bam"), # path with the sliced normal PCAWG .bams
    output:
        normal_bam=f"{output_path}"+"{sample}/normal.bam", # temp
        normal_bai=f"{output_path}"+"{sample}/normal.bam.bai",
        tumor_bam =f"{output_path}"+"{sample}/tumor.bam", # temp
        tumor_bai =f"{output_path}"+"{sample}/tumor.bam.bai"

    shell:
        'module load sambamcram/samtools/1.7; mkdir -p {output_path}/{wildcards.sample}; cp {input.normal_bam} {output.normal_bam}; samtools index {output.normal_bam};  cp {input.tumor_bam} {output.tumor_bam}; samtools index {output.tumor_bam} '

rule typing_tumor: 
    input:
        tumor_bam = f"{output_path}" + "{sample}/tumor.bam"
    output:
        typing_tumor_file = f"{output_path}" + "{sample}/typing_tumor/winners.hla.txt"
    threads: 2
    params: cluster_memory = "8G"
    run:
        path_base = os.path.join(output_path,f"{wildcards.sample}", "typing_tumor")
        print ("Working on... ",  path_base)
        if not (os.path.exists(path_base)):
            os.mkdir(path_base)
        if not (os.path.exists(path_base + "/tmp")):
            os.mkdir(path_base + "/tmp")

        input_path = os.path.dirname(input.tumor_bam)
        tumor_name_bam = os.path.basename(input.tumor_bam)
        if not (os.path.exists(path_base + "/tmp")):
            os.mkdir(path_base + "/tmp")
        shell('singularity exec -C -B {path_base}:/output/ -B {path_base}/tmp:/tmp/ -B {input_path}:/input  /hpc/local/CentOS7/cog/software/polysolver/polysolver4_new.sif /home/polysolver/scripts/shell_call_hla_type /input/{tumor_name_bam} Unknown 1 hg19 STDFQ 0 /output/') # We used the singulairty container of polysolver, the binary could also be used

rule typing_normal:
    input:
        normal_bam = ancient(f"{output_path}" + "{sample}/normal.bam")
    output:
        typing_normal_file = f"{output_path}" + "{sample}/typing_normal/winners.hla.txt"
    threads: 2
    params: cluster_memory = "8G"
    run:
        path_base = os.path.join(output_path,f"{wildcards.sample}","typing_normal")
        print ("Working on... ",  path_base)
        if not (os.path.exists(path_base)):
            os.mkdir(path_base)
        if not (os.path.exists(path_base+ "/tmp")):
            os.mkdir(path_base + "/tmp")

        input_path = os.path.dirname(input.normal_bam)
        normal_name_bam = os.path.basename(input.normal_bam)
        shell('singularity exec -C -B {path_base}:/output/ -B {path_base}/tmp:/tmp/ -B {input_path}:/input  /hpc/local/CentOS7/cog/software/polysolver/polysolver4_new.sif /home/polysolver/scripts/shell_call_hla_type /input/{normal_name_bam} Unknown 1 hg19 STDFQ 0 /output/') # We used the singulairty container of polysolver, the binary could also be used
        shell('rm -rf {path_base}/tmp')

rule typing_normal_xHLA:
    input:
        normal_bam = ancient(f"{output_path}" + "{sample}/normal.bam")
    output:
        outpath = directory(f"{output_path}" + "{sample}/typing_normal_xhla/"),
        typing_normal_file = f"{output_path}" + "{sample}/typing_normal_xhla/report-{sample}-hla.json"
    threads: 4
    params: cluster_memory = "32G"
    run:
        if not (os.path.exists(output.outpath)):
            os.mkdir(output.outpath)
        shell('run_hla_typing_xHLA.sh {input.normal_bam} {output.outpath} {wildcards.sample}')
        shell('rm {output.outpath}/*.bam')




rule typing_tumor_xHLA:
    input:
        tumor_bam = ancient(f"{output_path}" + "{sample}/tumor.bam")
    output:
        outpath = directory(f"{output_path}" + "{sample}/typing_tumor_xhla/"),
        typing_tumor_file = f"{output_path}" + "{sample}/typing_tumor_xhla/report-{sample}-hla.json"
    threads: 4
    params: cluster_memory = "32G"

    run:
        if not(os.path.exists(f"{output.outpath}")):
            os.mkdir(f"{output.outpath}")

        shell('run_hla_typing_xHLA.sh {input.tumor_bam} {output.outpath} {wildcards.sample}')
        shell('rm {output.outpath}/*.bam')


