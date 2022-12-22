import os, sys
import pandas as pd
import allel
import numpy as np
import click

input_path="../../external_data/background_genes.tsv"
genes_interest=list(pd.read_csv(input_path,sep="\t")["gene"].values)

valid_csqns=["structural_interaction_variant",
"chromosome_number_variation",
"exon_loss_variant",
"frameshift_variant",
"stop_gained",
"stop_lost",
"start_lost",
"splice_acceptor_variant",
"splice_donor_variant",
"rare_amino_acid_variant",
"missense_variant",
"disruptive_inframe_insertion",
"conservative_inframe_insertion",
"disruptive_inframe_deletion",
"conservative_inframe_deletion",
"5_prime_UTR_truncation+exon_loss_variant",
"3_prime_UTR_truncation+exon_loss",
"splice_branch_variant",
"splice_region_variant",
"stop_retained_variant",
"initiator_codon_variant"]


def is_valid(bases):
    v = ["A", "C", "T", "G"]
    correct = True
    for b in bases:
        if not (b in v):
            return False
    return correct


def count_mnvs(ref, alt):
    count = 0

    for r, a in zip(ref, alt):
        if len(r) == len(a) and len(r) > 0 and is_valid(r) and is_valid(a) and len(a) > 1:
            count += 1
    return count

def valid_csqn(values):
    o=[]
    for v in values:
        valid =False
        for c in valid_csqns:
            if c in v:
                o.append(True)
                valid=True
                break
        if not(valid):
            o.append(False)
    return o

def is_syn(values):
    o = []
    for v in values:
        if v.startswith("synonymous_variant"):
            o.append(True)
        else:
            o.append(False)
    return o



def load_somatic_mutations(path_somatic):

    df_muts = allel.vcf_to_dataframe(path_somatic,fields="*")
    print (df_muts.columns.values)
    df_muts = df_muts[df_muts["FILTER_PASS"]][
        ["CHROM", "POS", "REF", "ALT_1", "ANN", "BIALLELIC", "PON_COUNT", "PURPLE_AF",
         "PURPLE_MACN", "PURPLE_GERMLINE","PURPLE_VCN", "TNC", "SUBCL", "HOTSPOT","NEAR_HOTSPOT"]]
    df_muts = df_muts[~pd.isnull(df_muts["ANN"])]
    if df_muts.shape[0] >0:
        df_muts[["csq", "gene", "impact"]] = df_muts.apply(
            lambda row: pd.Series([row["ANN"].split("|")[1], row["ANN"].split("|")[3], row["ANN"].split("|")[2]]), axis=1)
        df_muts = df_muts[
            (df_muts["gene"].isin(genes_interest)) & ((df_muts["impact"].isin(["MODERATE", "HIGH"])) & (valid_csqn(df_muts["csq"])))]
        
    return df_muts


def check_cnvs(dels,candidates):
    result = []
    if dels.shape[0] >0:
        gs = list(dels["gene"].unique())
        for c in candidates:
            if c in gs:
                result.append(1)
            else:
                result.append(0)
    else:
        result = list(np.zeros(len(candidates)))
    return result

def check_ploidy(df_loh,column_cns=["integer_minor","integer_major"]):
    results = []
    for m in genes_interest:
        x=df_loh[df_loh["gene"]==m]
        results.append(",".join([str(x[column_cns[0]].values[0]),str(x[column_cns[1]].values[0])]))
    return results

def is_loh_purple(row,column_cn="minMinorAlleleCopyNumber"):
    if not (column_cn in row):
        column_cn = "minMinorAllelePloidy"
    ratio_diff=  row["minMajorAlleleCopyNumber"] - row["minMinorAlleleCopyNumber"]
    return (row["minMinorAlleleCopyNumber"]<0.3) & (row["integer_major"]>=1) & (ratio_diff>0.4)

def get_chr_arm(row,cytobands):
    end = row["end"]
    start = row["start"]
    c = row["chromosome"]

    p_lenght = cytobands.loc[(c, "p")]["length"]
    q_lenght = cytobands.loc[(c, "q")]["length"]
    if start < p_lenght and end < p_lenght:
        return p_lenght
    if start < p_lenght and end > p_lenght:
        return np.nanmin([p_lenght, q_lenght])
    elif start > p_lenght:
        return q_lenght

    return np.nanmin([p_lenght, q_lenght])

def load_cytoband():
    cytobands=pd.read_csv("cytoband_sizes.tsv",sep="\t")
    cytobands["chr"] = cytobands.apply(lambda row: row["chrom"][3:],axis=1)
    cytobands.set_index(["chr","arm"],inplace=True)
    return cytobands

def is_focal(row):
    l = row["len"]
    if l < row["chr_arm_lenght"]*0.75: # lower than 75% of arm lenght
        return True
    return False

def is_hfocal(row):
    l = row["len"]
    if l < 3*10**6: # lower than 3MBs
        return True
    return False

def is_nonfocal(row):
    l =row["len"]
    if l > row["chr_arm_lenght"]*0.75: # greater than 75% of arm lenght
        return True
    return False

def check_mutations(df_somatic_muts,list_genes):
    result = []
    if df_somatic_muts.shape[0]==0:
        return ""
    if "PURPLE_VCN" in df_somatic_muts.columns.values: # hartwig
        column_vcn, column_macn = "PURPLE_VCN", "PURPLE_MACN"
    else: # pcawg
        column_vcn, column_macn = "PURPLE_CN", "PURPLE_MAP"
    for gene in list_genes:
        g = df_somatic_muts[(df_somatic_muts["gene"] == gene)][["gene", "csq", "BIALLELIC", "PURPLE_AF",column_vcn,column_macn, "SUBCL"]].drop_duplicates()
        if g.shape[0]>0:
            for gene,csq,bi,af,cn,macn,subcl_likelihood in g.values:
                allelic_status="biallelic" if bi else "monoallelic"
                if np.isfinite(subcl_likelihood):
                    clonality_purple="clonal" if subcl_likelihood < 0.80 else "subclonal"
                else:
                    clonality_purple = "clonal"
                result.append(";".join([gene,csq,allelic_status,str(clonality_purple),"{:.2f}".format(af),"{:.2f}".format(cn),"{:.2f}".format(macn)]))
    return "___".join(result)

@click.command()
@click.option('--somatic_vcf',
              type=click.Path(exists=True),
              help="Input path of vcf with somatic mutations",
              required=True)
@click.option('--somatic_cnv',
              type=click.Path(exists=True),
              help="Input path of somatic cnv",
              required=True)
@click.option('--purity_info',
              type=click.Path(exists=True),
              help="Purity and ploidy of sample",
              required=True)
@click.option('--output_file',
              type=click.Path(),
              help="Output path (directory)",
              required=True)


def run(somatic_vcf, somatic_cnv, purity_info, output_file):

    # load lenght of cytoband
    cytobands=load_cytoband()

    # create output
    outpath=os.path.dirname(output_file)
    if not(os.path.exists(outpath)):
        os.mkdir(outpath)
    # process results
    sample= os.path.basename(somatic_vcf).split(".")[0]
    list_results = [sample]
    list_names = ["sample_id"]
    # Read sample ploidy and info
    df_purity = pd.read_csv(purity_info,  sep="\t")
    df_purity=df_purity[["purity", "ploidy", "diploidProportion", "gender", "wholeGenomeDuplication", "msStatus", "msIndelsPerMb","tml", "tmlStatus","tmbPerMb", "tmbStatus", "svTumorMutationalBurden"]]
    sample_ploidy = float(df_purity["ploidy"])
    # Read somatic cnv
    df_cnv = pd.read_csv(somatic_cnv,sep="\t")
    df_cnv = df_cnv[df_cnv["gene"].isin(genes_interest)][["gene", "chromosome","start","end","minCopyNumber", "maxCopyNumber", "minMinorAlleleCopyNumber","minRegionStart","minRegionEnd"]]
    df_cnv["integer_minor"] = df_cnv.apply(lambda row: round(row["minMinorAlleleCopyNumber"]),axis=1)
    df_cnv["minMajorAlleleCopyNumber"] = df_cnv.apply(lambda row: row["minCopyNumber"] - row["minMinorAlleleCopyNumber"], axis=1)
    df_cnv["integer_major"] = df_cnv.apply(lambda row: round(row["minMajorAlleleCopyNumber"]), axis=1)
    df_cnv["is_loh"] = df_cnv.apply(lambda row: is_loh_purple(row),axis=1)
    df_loh=df_cnv[df_cnv["is_loh"]==True]
    if df_loh.shape[0] >0:
        df_loh["len"]=df_loh.apply(lambda row: row["minRegionEnd"] - row["minRegionStart"],axis=1)
        df_loh["chr_arm_lenght"] = df_loh.apply(lambda row: get_chr_arm(row,cytobands), axis=1)
        df_loh["is_focal"] = df_loh.apply(lambda row: is_focal(row),axis=1) # < = 75% of chr lenght
        df_loh["is_hfocal"] = df_loh.apply(lambda row: is_hfocal(row),axis=1) # < = 75% of chr lenght
        df_loh["is_nonfocal"] = df_loh.apply(lambda row: is_nonfocal(row),axis=1) 

    ## Determine AMP and DEL events based on HMF driver calling rules
    threshold = float(3 * sample_ploidy)
    df_amp = df_cnv[(df_cnv["minCopyNumber"] > threshold)]
    # Determine deepdels
    df_dels = df_cnv[(df_cnv["minCopyNumber"] < 0.5)]
    # Read somatic mutations
    df_muts = load_somatic_mutations(somatic_vcf)
    df_muts = df_cnv[["gene","minCopyNumber","maxCopyNumber","minMinorAlleleCopyNumber"]].merge(df_muts,how="right")
    df_muts["POS"] = df_muts["POS"].astype(int)
    df_muts["CHROM"] = df_muts["CHROM"].astype(str)
    muts_genes = check_mutations(df_muts, genes_interest)
    list_results += [muts_genes]
    list_names += ["muts_genes"]
    

    # deletions
    dels_genes = check_cnvs(df_dels, genes_interest)
    list_results += dels_genes
    list_names += ["del_" + m for m in genes_interest]
    # amplifications
    amps_genes = check_cnvs(df_amp, genes_interest)
    list_results += amps_genes
    list_names += ["amp_" + m for m in genes_interest]

    # loh
    amps_genes = check_cnvs(df_loh, genes_interest)
    list_results += amps_genes
    list_names += ["loh_" + m for m in genes_interest]

    # annotate ploidys
    loh_minor_major_p = check_ploidy(df_cnv)
    list_results += loh_minor_major_p
    list_names += ["ploidy_minor_major_" + m for m in genes_interest]
    
    # annotate focality of loh events
    if df_loh.shape[0] >0:
        loh_focal=",".join(list(df_loh[df_loh["is_focal"]]["gene"].values))
        loh_hfocal=",".join(list(df_loh[df_loh["is_hfocal"]]["gene"].values))
        loh_non_focal=",".join(list(df_loh[df_loh["is_nonfocal"]]["gene"].values))
    else:
        loh_focal,loh_non_focal,loh_hfocal="","",""
    list_results+=[loh_focal,loh_non_focal,loh_hfocal]
    list_names+=["loh_focal","loh_nonfocal","loh_hfocal"]
    
    # save data
    df_total = pd.DataFrame([list_results],columns=list_names)
    df_total.to_csv(output_file,sep="\t",index=False)
    
    # Somatic mutations

if __name__ == '__main__':
    run()