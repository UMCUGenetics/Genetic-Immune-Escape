import glob
import pandas as pd
from tqdm import tqdm
l=[]
for report in tqdm(glob.glob("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/infiltration/*.immune.infiltration.tsv")):
    df = pd.read_csv(report,sep="\t")
    l.append(df)
df_total = pd.concat(l)
df_total.to_csv("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/infiltration/hmf_infiltration.tsv.gz",sep="\t",index=False,compression="gzip")



