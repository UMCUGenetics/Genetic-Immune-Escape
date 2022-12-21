import os,shutil
import glob
import pandas as pd


# Config
output_path=config["o"]
donors_file =config["i"]

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
        expand(f"{output_path}" + "{sample}/typing_tumor/winners.hla.txt",sample=samples),
        expand(f"{output_path}" + "{sample}/typing_normal/winners.hla.txt",sample=samples),
        expand(f"{output_path}" + "{sample}/typing_normal_xhla/report-{sample}-hla.json",sample=samples),
        expand(f"{output_path}" + "{sample}/typing_tumor_xhla/report-{sample}-hla.json", sample=samples),
 

rule select_hla_locus: # prepare the input bams
    input:
        tumor_bam="/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/hla/{sample}.hla.bam",
    output:
        normal_bam=f"{output_path}"+"{sample}/normal.bam", # temp
        normal_bai=f"{output_path}"+"{sample}/normal.bam.bai",
        tumor_bam =f"{output_path}"+"{sample}/tumor.bam", # temp
        tumor_bai =f"{output_path}"+"{sample}/tumor.bam.bai"

    run:
        normal_id=d_keys[wildcards.sample]
        normal_bam=input.tumor_bam.replace(wildcards.sample,normal_id)
        shell('module load sambamcram/samtools/1.7; mkdir -p {output_path}/{wildcards.sample}; cp {normal_bam} {output.normal_bam}; samtools index {output.normal_bam};  cp {input.tumor_bam} {output.tumor_bam}; samtools index {output.tumor_bam} ')

rule typing_tumor: # with Polysolver
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
        shell('singularity exec -C -B {path_base}:/output/ -B {path_base}/tmp:/tmp/ -B {input_path}:/input  /hpc/local/CentOS7/cog/software/polysolver/polysolver4_new.sif /home/polysolver/scripts/shell_call_hla_type /input/{tumor_name_bam} Unknown 1 hg19 STDFQ 0 /output/')

rule typing_normal:  # with Polysolver
    input:
        normal_bam = f"{output_path}" + "{sample}/normal.bam"
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
        shell('singularity exec -C -B {path_base}:/output/ -B {path_base}/tmp:/tmp/ -B {input_path}:/input  /hpc/local/CentOS7/cog/software/polysolver/polysolver4_new.sif /home/polysolver/scripts/shell_call_hla_type /input/{normal_name_bam} Unknown 1 hg19 STDFQ 0 /output/')
        shell('rm -rf {path_base}/tmp')

rule typing_normal_xHLA: # xHLA
    input:
        normal_bam = f"{output_path}" + "{sample}/normal.bam"
    threads: 4
    params: cluster_memory = "16G"
    output:
        outpath = directory(f"{output_path}" + "{sample}/typing_normal_xhla/"),
        typing_normal_file = f"{output_path}" + "{sample}/typing_normal_xhla/report-{sample}-hla.json"
    run:
        if not (os.path.exists(output.outpath)):
            os.mkdir(output.outpath)
        shell('cd {output.outpath}; /home/cog/fmartinez/scripts/organoids_CRC/run_hla_typing_xHLA.sh {input.normal_bam} {output.outpath} {wildcards.sample}')
        shell('rm {output.outpath}/*.bam')




rule typing_tumor_xHLA: # xHLA
    input:
        tumor_bam = f"{output_path}" + "{sample}/tumor.bam"
    threads: 4
    params: cluster_memory = "16G"
    output:
        outpath = directory(f"{output_path}" + "{sample}/typing_tumor_xhla/"),
        typing_tumor_file = f"{output_path}" + "{sample}/typing_tumor_xhla/report-{sample}-hla.json"
    run:
        if not(os.path.exists(f"{output.outpath}")):
            os.mkdir(f"{output.outpath}")

        shell('cd {output.outpath}; /home/cog/fmartinez/scripts/organoids_CRC/run_hla_typing_xHLA.sh {input.tumor_bam} {output.outpath} {wildcards.sample}')
        shell('rm {output.outpath}/*.bam')







