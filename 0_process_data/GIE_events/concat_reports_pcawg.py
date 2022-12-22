import glob
import pandas as pd
from tqdm import tqdm
l=[]
path_pcawg_output="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/hla_events/"
path_global_output="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/immune_escape/"
for report in tqdm(glob.glob(f"{path_pcawg_output}/*/report_escape/report_status_immune_escape.tsv.gz")): 
    df = pd.read_csv(report,sep="\t")
    l.append(df)
df_total = pd.concat(l)
df_total.to_csv(f"{path_global_output}/pcawg_reports_lilac.tsv.gz",sep="\t",index=False,compression="gzip") 
