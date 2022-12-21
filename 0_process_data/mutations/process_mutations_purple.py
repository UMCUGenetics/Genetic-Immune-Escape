import os, sys
import pandas as pd
import allel
import numpy as np
import click
from tqdm import tqdm
import re
import io
import gzip


def read_af_ann(t):
    try:
        purple_af=float(re.search("PURPLE_AF=([0-9\.]+);",t).group(1))
        purple_vcn=float(re.search("PURPLE_VCN=([0-9\.]+);",t).group(1))
        purple_cn=float(re.search("PURPLE_CN=([0-9\.]+);",t).group(1))
    except:
        purple_af,purple_vcn,purple_cn=np.nan,np.nan,np.nan
    try:
        ann=t.split(";")[0].split("=")[1]
    except IndexError:
        ann = ""
    
    if "|" in ann:
        csq,gene,impact=ann.split("|")[1],ann.split("|")[3],ann.split("|")[2]
    else:
        csq,gene,impact = "","",""
    return pd.Series([purple_af,purple_vcn,purple_cn,csq,gene,impact])



def vcf_reader(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={
            '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str
        }, sep='\t', low_memory=False,
    ).rename(columns={'#CHROM': 'CHROM'})


@click.command()
@click.option('--file_input',
              type=click.Path(exists=True),
              help="Input path of somatic cnv",
              required=True)
@click.option('--output_file',
              type=click.Path(),
              help="Output file",
              required=True)

def run(file_input, output_file):

    sample = os.path.basename(file_input).split(".")[0]
    df_muts = vcf_reader(file_input)
    df_muts = df_muts[
                ["CHROM", "POS", "REF", "ALT","FILTER",sample,"INFO"]].rename(columns={sample:"sample_info"})
    df_muts = df_muts[df_muts["FILTER"]!="PON"]
    df_muts[["PURPLE_AF","PURPLE_VCN","PURPLE_CN","csq","gene","impact"]] = df_muts.apply(lambda row: read_af_ann(row["INFO"]),axis=1)
    df_muts["AF_raw"] = df_muts.apply(lambda row: float(row["sample_info"].split(":")[2]),axis=1)
    df_muts.to_csv(output_file,sep="\t",index=False,compression="gzip")



if __name__ == '__main__':
    run()