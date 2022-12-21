import pandas as pd
import os
import numpy as np
import click
from tqdm import tqdm

# Extracted from Danaher et al. Journal for ImmunoTherapy of Cancer (2017) 5:18
dict_info={
# immune markers
"bcells":["BLK","CD19","MS4A1","TNFRSF17","FCRL2","KIAA0125","PNOC","SPIB","TCL1A"],
"cd45":["PTPRC"],
"cd8":["CD8A","CD8B"],
"cytox":["CTSW","GNLY","GZMA","GZMB","GZMH","KLRB1","KLRD1","KLRK1","PRF1","NKG7"],
"dc":["CCL13","CD209","HSD11B1"],
"exhausted_cd8":["CD244","EOMES","LAG3","PTGER4"],
"macroph":["CD163","CD68","CD84","MS4A4A"],
"mast_cells":["MS4A2","TPSAB1","CPA3","HDC","TPSB2"],
"neutrophil":["CSF3R","S100A12","CEACAM3","FCAR","FCGR3B","FPR1","SIGLEC5"],
"nk_cd56n":["IL21R","KIR2DL3","KIR3DL1","KIR3DL2"],
"nk":["NCR1","XCL2","XCL1"],
"tcells":["CD3D","CD3E","CD3G","CD6","SH2D1A","TRAT1"],
"th1":["TBX21"],
"treg":["FOXP3"],
"infiltration_davoli" : ["CD247","CD2","CD3E","GZMH","NKG7","PRF1","GZMK"], # davoli et al Science 2017
"cd4_davoli": ["IGFBP4","ITM2A","AMIGO2","TRAT1","CD40LG","ICOS"], # davoli et al Science 2017
"cd8_davoli" : ["CD3E","GZMK","CXCR3","BCL11B","IL7R","KLRG1"],
"t_cell_grasso": ["CCL2", "CCL3", "CCL4", "CXCL9", "CXCL10", "CD8A", "HLA-DOB", "HLA-DMB", "HLA-DOA", "GZMK",
                             "ICOS", "IRF1"], #Grasso et al. (DOI: 10.1158/2159-8290.CD-17-132)
# Other genes of interest
"hla-a":["HLA-A"],
"B2M":["B2M"],
"NLRC5":["NLRC5"],
"hla-e":["HLA-E"],
"hla-b":["HLA-B"],
"hla-c":["HLA-C"],
"PDL1":["CD274"],
"CTLA4":["CTLA4"],
"endopeptidases" : ["ERAP1","ERAP2"],
"proteasome" : ["PSME1","PSME2","PSME3","PSMB5","PSMB6","PSMB7"],
"immunoproteasome" : ["PSMB8","PSMB9","PSMB10"],
"ifn-gamma":["IFNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11", "IDO1", "PRF1", "GZMA", "HLA-DRA"]
           }

def get_scores_infiltration(df_rna,column):
    df_rna["log2_rna"] = np.log2(df_rna[column]+0.1)
    l = []
    for key in dict_info:
        list_genes = dict_info[key]
        score=np.nanmean(df_rna[df_rna["GeneName"].isin(list_genes)].groupby(["Transcript stable ID","GeneName"],as_index=False).agg({"log2_rna":np.nanmax})["log2_rna"].values)
        l.append([key, score])

    df = pd.DataFrame(l, columns=["marker", "value"]).set_index("marker").T

    return df

@click.command()
@click.option('--rna_data',
              type=click.Path(exists=True),
              help="Input path of dataframe with all samples",
              required=True)
@click.option('--ensembl_mapping',
              type=click.Path(exists=True),
              help="Mapping of ensembl ids to hugo symbol",
              required=True)
@click.option('--metadata_pcawg',
              type=click.Path(exists=True),
              help="Metadata pcawg",
              required=True)
@click.option('--output_file',
              type=click.Path(),
              help="Output file",
              required=True)
def run(rna_data, ensembl_mapping, metadata_pcawg, output_file):
    # read the data
    df = pd.read_csv(rna_data,sep="\t")
    df["Transcript stable ID"] = df.apply(lambda row: row["Transcript"].split(".")[0],axis=1)
    df.drop(columns=["Transcript"],inplace=True)
    df=df.melt(id_vars=["Transcript stable ID"],var_name="aliquot_id",value_name="tpm")
    # get donor id

    df_metadata = pd.read_csv(metadata_pcawg,sep="\t")
    df=df.merge(df_metadata[["aliquot_id", "icgc_donor_id"]].drop_duplicates())
    # get mapping ensembl
    df_mapping = pd.read_csv(ensembl_mapping,sep="\t")[["Transcript stable ID","HGNC symbol"]].drop_duplicates().rename(columns={"HGNC symbol":"GeneName"})
    df = df_mapping.merge(df)
    print (df.head())
    # for sample iterate
    df.set_index("icgc_donor_id",inplace=True)
    l=[]
    for sample in tqdm(set(df.index)):
        df_tmp = df.loc[sample]
        res = get_scores_infiltration(df_tmp, "tpm")
        res["sampleId"] = sample
        l.append(res)
    df_total = pd.concat(l)
    # save it
    df_total.to_csv(output_file,sep="\t",index=False,compression="gzip")



if __name__ == '__main__':
    run()