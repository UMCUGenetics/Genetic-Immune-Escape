import pandas as pd
import os


# Config
output_path=config["o"]
dataset=config["d"]
input_path=config["i"] # this the json file created by the "create_groups_for_positive_selection.ipynb"

with open(input_path,'r') as f:
    d_info=json.load(f)
ttypes=list(d_info.keys())
for ttype in ttypes:
    # create list of patients
    file_input = f"{output_path}" + f"/{ttype}.json"
    with open(file_input, 'w') as f:
        json.dump(d_info[ttype], f)


names =["ignore_loh_focal", "ignore_loh_nonfocal","ignore_loh_hfocal", # LOH CNVs
        "ignore_deepdel_focal", "ignore_deepdel_nonfocal","ignore_deepdel_hfocal",  # Deletions CNVs
        "ignore_amp_nonfocal","ignore_amp_focal","ignore_amp_hfocal"] # Amplifications CNVs

# Rules
rule all:
    input:
        expand(f"{output_path}" + "/{conf}/{ttype}.tsv.gz",ttype=ttypes,conf=names)

rule run_pos_selection_LOH:
    input:
        input_path = ancient(f"{output_path}" + "/{ttype}.json")
    threads: 8
    params: cluster_memory = "16G"
    output:
        output=f"{output_path}" + "/{conf}/{ttype}.tsv.gz"
    run:
        wgd_status,type_analysis,focal=wildcards.conf.split("_")
        flag=""
        p=f"{output_path}" + f"/{wildcards.conf}/"
        if not(os.path.exists(p)):
            os.mkdir(p)

        shell('set +eu '
              ' && PS1=dummy '
              ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && ' \
              'python cnv/positive_selection_CNV_shuffle.py --samples_file {input.input_path} ' \
              '--tumor_type "{wildcards.ttype}" --type_analysis {type_analysis} {flag} --focal {focal}  --output_file {output.output} --dataset {dataset} ')

