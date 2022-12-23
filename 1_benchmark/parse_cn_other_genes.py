from tqdm import tqdm
import pandas as pd
import os
df_genes = pd.read_csv("../external_data/selected_regions.tsv",sep="\t")
selected_genes = df_genes[~(df_genes["gene"].str.contains("HLA-",regex=False))]["gene"].values
df_meta = pd.read_csv("../metadata/dataset_metadata_supp_table3.tsv",sep="\t") # Supp Table 3

for dataset in ["hmf","pcawg"]:
    l=[]
    if dataset=="pcawg":
        donors=df_meta[(df_meta["cohort"]=="PCAWG")&(df_meta["is_selected"]==True)]["sample_id"].values
    else:
        donors=df_meta[(df_meta["cohort"]=="Hartwig")&(df_meta["is_selected"]==True)]["sample_id"].values
    for s in donors:
        if dataset == "pcawg":
            filein=f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{s}-from-jar/purple53/{s}T.purple.cnv.gene.tsv"
        else:
            filein = f"/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/{st}/purple/{s}.purple.cnv.gene.tsv"

        if not (os.path.exists(filein)):
            continue
        dh = pd.read_csv(filein,sep="\t")
        dh=dh[dh["gene"].isin(selected_genes)][["gene","minMinorAlleleCopyNumber","minCopyNumber"]]
        dh["sampleId"] = s
        l.append(dh)

    p = pd.concat(l)
    p.to_csv(f"data/cn_data_random_genes_{dataset}.tsv",sep="\t",index=False,compression="gzip")
                                                  
