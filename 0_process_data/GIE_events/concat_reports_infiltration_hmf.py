import glob
import pandas as pd
from tqdm import tqdm
l=[]
path_hmf_output="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/infiltration/"
path_infiltration_output="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/infiltration/"
for report in tqdm(glob.glob(f"{path_hmf_output}*.immune.infiltration.tsv")):
    df = pd.read_csv(report,sep="\t")
    l.append(df)
df_total = pd.concat(l)
df_total.to_csv(f"{path_infiltration_output}hmf_infiltration.tsv.gz",sep="\t",index=False,compression="gzip")



