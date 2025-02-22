import os,shutil
import glob

import pandas as pd

# Aux function
def prepare_purity_file(sample,path_purity,output_file):
    df_p = pd.read_csv(path_purity,sep="\t")
    purity = float(df_p["purity"])
    ploidy=float(df_p["ploidy"])
    with open(output_file,'w') as f:
        f.write("Ploidy\ttumorPurity\ttumorPloidy\n")
        f.write(f"tumor_{sample}\t{2}\t{purity}\t{ploidy}\n")
    return

def read_hlas_xhla(path):

    with open(path, 'r') as f:
        d = json.load(f)["hla"]
    list_normal_total = [get_4digit_xhla(hla) for hla in d["alleles"]]
    list_normal_I = list([hla for hla in list_normal_total if hla[4] != "D"]) # only class I at this moment
    return list_normal_I

def get_4digit_xhla(hla):
    v = "HLA-" + hla.replace("*", "")
    if v[-1] == "N": # this allele is raising an error...
        return v[0:-1] # remove the "N" of the final...
    return v

def create_divergence_file(path_xhla,output_path,sample):
    hlas = read_hlas_xhla(path_xhla)
    l,s_total = [], set()
    for type_hla in ["HLA-A","HLA-B","HLA-C"]:
        l.append([type_hla])
        for h in hlas:
            if type_hla in h:
                t=h.split("-")[1].replace(":","")
                l[-1].append(t)
                s_total.add(t)
    l.append(["germline_diversity"]+list(s_total))
    df=pd.DataFrame(l,columns=["IDs"]+["Allele"+str(f) for i in range(len(l[-1])-1)])
    df.to_csv(os.path.join(output_path,f"{sample}.tsv"),sep="\t",index=False)


# Config
output_path=config["o"]
donors_file =config["i"]

with open(donors_file) as f:
    donors = f.read().splitlines()
print (donors[0])
donors_s = []
for donor in donors:
    if os.path.exists(f"/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/lilac/{donor}/{donor}.lilac.csv.gz") and (os.path.exists(f"{output_path}" + f"{donor}/loh/{donor}.5.DNA.HLAlossPrediction_CI.xls")):
        donors_s.append(donor)


print (len(donors_s))
# Rules
rule all:
    input:
        expand(f"{output_path}" + "{donor}/typing_tumor/winners.hla.txt",donor=donors),
        expand(f"{output_path}" + "{donor}/typing_normal/winners.hla.txt",donor=donors),
        expand(f"{output_path}" + "{donor}/typing_tumor/{donor}.mutect.filtered.annotated",donor=donors),
        expand(f"{output_path}" + "{donor}/loh/{donor}.5.DNA.HLAlossPrediction_CI.xls",donor=donors),
        expand(f"{output_path}" + "{donor}/typing_normal_xhla/report-{donor}-hla.json",donor=donors),
        expand(f"{output_path}" + "{donor}/typing_tumor_xhla/report-{donor}-hla.json", donor=donors),
        expand(f"{output_path}" + "{donor}/summary/"+"somatic_mutations_APP.tsv.gz",donor=donors_s),
        expand(f"{output_path}" + "{donor}/divergence/{donor}.IndividualDivergence.txt",donor=donors),
        expand(f"{output_path}" + "{donor}/report_escape/"+"report_status_immune_escape.tsv.gz",donor=donors_s),



rule select_hla_locus:
    input:
        tumor_bam=ancient("/hpc/cuppen/shared_resources/PCAWG/pipeline5/HLA/{sample}T_6:29854528-32726735.bam"),
        normal_bam = ancient("/hpc/cuppen/shared_resources/PCAWG/pipeline5/HLA/{sample}R_6:29854528-32726735.bam"),
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
        shell('singularity exec -C -B {path_base}:/output/ -B {path_base}/tmp:/tmp/ -B {input_path}:/input  /hpc/local/CentOS7/cog/software/polysolver/polysolver4_new.sif /home/polysolver/scripts/shell_call_hla_type /input/{tumor_name_bam} Unknown 1 hg19 STDFQ 0 /output/')

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
        shell('singularity exec -C -B {path_base}:/output/ -B {path_base}/tmp:/tmp/ -B {input_path}:/input  /hpc/local/CentOS7/cog/software/polysolver/polysolver4_new.sif /home/polysolver/scripts/shell_call_hla_type /input/{normal_name_bam} Unknown 1 hg19 STDFQ 0 /output/')
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




rule call_hla_mutations:
    input:
        normal_bam = ancient(f"{output_path}" + "{sample}/normal.bam"),
        tumor_bam = ancient(f"{output_path}" + "{sample}/tumor.bam"),
        winners= f"{output_path}" + "{sample}/typing_normal/winners.hla.txt",
        winners_tumor = f"{output_path}" + "{sample}/typing_tumor/winners.hla.txt"
    output:
        mutect_output = f"{output_path}" + "{sample}/typing_tumor/{sample}.mutect.filtered.annotated",
        strelka_indels = f"{output_path}" + "{sample}/typing_tumor/{sample}.strelka_indels.filtered.annotated"
    threads: 2
    params: cluster_memory = "12G"
    run:
        input_path = os.path.dirname(input.tumor_bam)
        output_path_n = os.path.dirname(output.mutect_output)
        shell("mkdir -p {output_path_n}/tmp")
        shell ("singularity exec -C -B {output_path_n}:/output/ -B {output_path_n}/tmp:/tmp/  -B {input_path}:/input /hpc/local/CentOS7/cog/software/polysolver/polysolver4_new.sif bash /home/polysolver/scripts/shell_call_hla_mutations_from_type /input/normal.bam  /input/tumor.bam /input/typing_normal/winners.hla.txt hg19 STDFQ /output/ {wildcards.sample}")
        shell ("cd {output_path}/{wildcards.sample}/; tar -czvf typing_tumor.tar.gz typing_tumor/")
        shell ("mv {output_path}/{wildcards.sample}/typing_tumor.tar.gz {output_path_n}")
        shell ("echo 'singularity exec -C -B {output_path_n}:/output/ -B {output_path_n}/tmp:/tmp/ /hpc/local/CentOS7/cog/software/polysolver/polysolver4_new.sif bash /home/polysolver/scripts/shell_annotate_hla_mutations {wildcards.sample} /output/typing_tumor.tar.gz /output'")
        shell ("singularity exec -C -B {output_path_n}:/output/ -B {output_path_n}/tmp:/tmp/ /hpc/local/CentOS7/cog/software/polysolver/polysolver4_new.sif bash /home/polysolver/scripts/shell_annotate_hla_mutations {wildcards.sample} /output/typing_tumor.tar.gz /output")
        shell ("rm -rf {output_path}/{wildcards.sample}/typing_tumor/tmp")  # clean the data
        shell ("rm -rf {output_path}/{wildcards.sample}/typing_tumor/*.tar.gz")  # clean the data


rule run_lohhla:
    input:
        normal_bam = ancient(f"{output_path}" + "{sample}/normal.bam"),
        tumor_bam = ancient(f"{output_path}" + "{sample}/tumor.bam"),
        winners= f"{output_path}" + "{sample}/typing_normal/winners.hla.txt",
        purity= "/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{sample}-from-jar/purple53/{sample}T.purple.purity.tsv"
    output:
        loh=directory(output_path+"{sample}/loh"),
        loh_file=f"{output_path}" + "{sample}/loh/{sample}.5.DNA.HLAlossPrediction_CI.xls"
    threads: 2
    params: cluster_memory = "8G"
    run:
        if not (os.path.exists(output.loh)):
            os.mkdir(output.loh)
        purity_file = os.path.join(output.loh,"purity_info.txt") # this should change in real tumor samples from patients
        prepare_purity_file(wildcards.sample,input.purity,purity_file)
        shell("run_lohhla.sh {input.normal_bam} {input.tumor_bam} {output.loh} {wildcards.sample} {purity_file} {input.winners}")


rule get_summary:
    input:
        loh_hla= f"{output_path}" + "{sample}/loh/{sample}.5.DNA.HLAlossPrediction_CI.xls",
    params: cluster_memory = "8G"
    output:
        output_somatic_dir = directory(f"{output_path}" + "{sample}/summary/"),
        output_somatic_summary = f"{output_path}" + "{sample}/summary/"+"somatic_mutations_APP.tsv.gz"


    run:
        somatic_hla_mutations = f"{output_path}" + f"{wildcards.sample}/typing_tumor/{wildcards.sample}.mutect.filtered.annotated",
        somatic_hla_indels= f"{output_path}" + f"{wildcards.sample}/typing_tumor/{wildcards.sample}.strelka_indels.filtered.annotated",
        gs = glob.glob(f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/*/germline_caller/*.germline.vcf.gz")
        file_germline = ""
        if len(gs)>0:
            file_germline=f"--germline_vcf {gs[0]}"
        print ("file germline", f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/*/germline_caller/*.germline.vcf.gz")
        # files input
        somatic_vcf= f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/purplesoft3.3/{wildcards.sample}T.purple.somatic.vcf.gz"
        somatic_cnv = f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/purplesoft3.3/{wildcards.sample}T.purple.cnv.gene.tsv"
        sample_info = f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/purplesoft3.3/{wildcards.sample}T.purple.purity.tsv"

        shell('set +eu '
        ' && PS1=dummy '
        ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && ' \
        'python analyze_events_APP.py --somatic_vcf {somatic_vcf} ' \
        '--somatic_cnv {somatic_cnv} --purity_info {sample_info} --somatic_hla_mutations {somatic_hla_mutations} ' \
        '--somatic_hla_indels {somatic_hla_indels} --loh_hla {input.loh_hla} --outpath {output.output_somatic_dir} {file_germline} ')
        gs = glob.glob(f"rm -r {output_path}/{wildcards.sample}/*.bam")
        if len(gs) > 0:
            shell(f"rm -r {output_path}/{wildcards.sample}/*.bam")
            shell(f"rm -r {output_path}/{wildcards.sample}/*.bai")

rule compute_divergence:
    input:
        typing_normal= f"{output_path}" + "{sample}/typing_normal_xhla/report-{sample}-hla.json"

    output:
        divergence_info = f"{output_path}" + "{sample}/divergence/{sample}.IndividualDivergence.txt"
    run:
        output_path = os.path.dirname(output.divergence_info)
        if not (os.path.exists(output_path)):
            os.mkdir(output_path)
        create_divergence_file(input.typing_normal,output_path,wildcards.sample)
        path_script = "./hla_divergence/"
        shell(f"perl {path_script}CalculateIndividualDivergence.pl -d {path_script}AAdistMatrix_Grantham.cnv  -f {path_script}HLA_ClassI_CWDonly.fas -i {os.path.join(output_path,wildcards.sample+'.tsv')} -o {output_path} -s {wildcards.sample} ")

rule prepare_report_sample:
    input:
        input_somatic_muts = f"{output_path}" + "{sample}/summary/"+"somatic_mutations_APP.tsv.gz",
        input_dir = f"{output_path}" + "{sample}/summary/",
    output:
        report_output = f"{output_path}" + "{sample}/report_escape/"+"report_status_immune_escape.tsv.gz"
    run:
        timing_data = f"{output_path}" + f"/../timing/{wildcards.sample}/"+f"{wildcards.sample}.mutationaltiming.tsv.gz" # this is only for internal usage, not part of the manuscript
        fusion_data = f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{wildcards.sample}-from-jar/linxsoft1.17/{wildcards.sample}T.linx.fusion.tsv"
        if os.path.exists(fusion_data):
            fusion_data = "--path_fusions "+fusion_data
        else:
            fusion_data = ""
        divergence_info = f"{output_path}" + f"{wildcards.sample}/divergence/{wildcards.sample}.IndividualDivergence.txt",
        shell('set +eu '
        ' && PS1=dummy '
        ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && ' \
        'python create_report_summary.py --path_diversity {divergence_info} --input_dir {input.input_dir} {fusion_data} --timing {timing_data} --sample {wildcards.sample} --output_file {output.report_output} --path_lilac /hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/lilac/{wildcards.sample}/')

