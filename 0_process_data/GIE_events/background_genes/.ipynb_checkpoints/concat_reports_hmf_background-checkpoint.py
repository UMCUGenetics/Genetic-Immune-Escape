import glob
import pandas as pd
from tqdm import tqdm
l=[]
output_hmf="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/hla_events/"
ouptut_global="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/immune_escape/"
for report in tqdm(glob.glob(f"{output_hmf}/*/report_background_genes/genomic_alterations.tsv.gz")):
    df = pd.read_csv(report,sep="\t")
    l.append(df)
df_total = pd.concat(l)
df_total.to_csv(f"{ouptut_global}/hmf_reports_background_genes.tsv.gz",sep="\t",index=False,compression="gzip")
