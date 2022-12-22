import pandas as pd
import os
import numpy as np
import click
import itertools


# hla locus
mhc_I = ["HLA-A","HLA-B","HLA-C"]
mhc_I_other = ["HLA-E","HLA-H","HLA-G","HLA-F"]
mhc_II = ["HLA-DPA1","HLA-DPB1","HLA-DQA1", "HLA-DQA2","HLA-DQB1", "HLA-DQB2", 	"HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRA"] # HLA-D*
# APP related proteins
transport_mhc = ["TAP1","TAP2","TAPBP","TAPBPL"]
scaffold_mhc = ["B2M","CANX","CALR","PDIA3"]
endopeptidases = ["ERAP1","ERAP2"]
proteasome = ["PSME1","PSME2","PSME3","PSMB5","PSMB6","PSMB7"]
immunoproteasome = ["PSMB8","PSMB9","PSMB10"]
# chaperones global
chaperones = ["HSPA1A","HSPA1B","HSP90AA1","STUB1","HSPBP1"]
# nmd factors
nmd = ["UPF1","UPF2","SMG1"]
#interferon
interferon=["JAK1","JAK2","STAT1","IRF1","IRF2","APLNR","IFNGR1", "IFNGR2"]
#TFs HLA locus
tfs=["NLRC5","RFX5","CIITA"]
# elongation factors
elongation = ["EIF5B","EIF4E","EIF4G1","EIF4A1","EIF5","EIF1","EIF3A"]
stop_readthrough = ["ETF1","GSPT1"]
# breaks
inhibitors = ["CTLA4","CD274"]

# CD58
cd58_nk = ["CD58"]

# Epigenetic regulators
epigenetic_regulators = ["SETDB1"]

coding_csq = ["frameshift_variant","stop_lost", "inframe_deletion", "inframe_insertion","stop_gain","start_lost","stop_lost","missense_variant","structural_interaction_variant",""]


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

def check_cnvs_condensed(dels,candidates):
    result = []
    if dels.shape[0] >0:
        gs = list(dels["gene"].unique())
        for c in candidates:
            if c in gs:
                result.append(c)
    return ",".join(result)

def check_loh(df_loh,column_cn="minMinorAlleleCopyNumber"):
    results = []
    for m in mhc_I:
        if not(column_cn in df_loh.columns.values):
            column_cn="minMinorAllelePloidy"
        if not ("PVal_unique" in df_loh.columns.values):
            results.append("-")
            continue
        lost=df_loh[(df_loh["gene"]==m)&(df_loh["integer_minor"]==0)&(df_loh["PVal_unique"]<0.01)]
        if lost.shape[0] >0:
            results.append(lost["LossAllele"].values[0])
        else:
            results.append("-")
    return results


def check_ploidy(df_loh,column_cn="integer_minor"):
    results = []
    for m in mhc_I:
        results.append(int(df_loh[df_loh["gene"]==m][column_cn].values[0]))

    return results



def is_loh_purple(df_loh,column_cn="minMinorAlleleCopyNumber"):
    if not (column_cn in df_loh.columns.values):
        column_cn = "minMinorAllelePloidy"
    df_loh["ratio_diff"] = df_loh.apply(lambda row: row["minMajorAlleleCopyNumber"] - row["minMinorAlleleCopyNumber"],axis=1 )
    return df_loh[(df_loh["integer_minor"]==0)&(df_loh["minMajorAlleleCopyNumber"]>0.5)&(df_loh["ratio_diff"]>0.4)].shape[0] >0

def is_imbalance_purple(df_loh,column_cn="minMinorAlleleCopyNumber"):
    if not (column_cn in df_loh.columns.values):
        column_cn = "minMinorAllelePloidy"
    df_loh["ratio"] = df_loh.apply(lambda row: (row["integer_minor"]) / (row["integer_major"]+10**-9),axis=1)
    df_loh["ratio_diff"] = df_loh.apply(lambda row: row["minMajorAlleleCopyNumber"] - row["minMinorAlleleCopyNumber"],
                                        axis=1)
    return df_loh[(df_loh["integer_minor"]>0)&(df_loh["ratio"]<.9)&(df_loh["ratio_diff"]>0.4)].shape[0] >0

def is_loh_lohhla(df_loh):
    if not ("PVal_unique" in df_loh.columns.values):
        return False
    return df_loh[df_loh["PVal_unique"]<0.01].shape[0] > 0

def check_loh_imbalance(df_loh,column_cn="minMinorAlleleCopyNumber"):
    results = []
    for m in mhc_I:
        if not("PVal_unique" in df_loh.columns.values):
            results.append("-")
            continue
        if not(column_cn in df_loh.columns.values):
            column_cn="minMinorAllelePloidy"
        lost=df_loh[(df_loh["gene"]==m)&(df_loh["integer_minor"]>0)&(df_loh["PVal_unique"]<0.01)] # then there is an imbalance
        if lost.shape[0] >0:
            results.append(lost["LossAllele"].values[0])
        else:
            results.append("-")
    return results

def check_mutations(df_somatic_muts,list_genes):
    result = []
    if df_somatic_muts.shape[0]==0:
        return ""
    if "PURPLE_VCN" in df_somatic_muts.columns.values: # hmf
        column_vcn, column_macn = "PURPLE_VCN", "PURPLE_MACN"
    else: # pcawg
        column_vcn, column_macn = "PURPLE_CN", "PURPLE_MAP"
    for gene in list_genes:
        if "timing_class" in df_somatic_muts.columns.values:
            g=df_somatic_muts[(df_somatic_muts["gene"]==gene)][["gene","csq","BIALLELIC","PURPLE_AF",column_vcn,column_macn,"germline_variant","SUBCL","timing_class"]].drop_duplicates()
        else:
            g = df_somatic_muts[(df_somatic_muts["gene"] == gene)][["gene", "csq", "BIALLELIC", "PURPLE_AF",column_vcn,column_macn, "germline_variant", "SUBCL"]].drop_duplicates()
            g["timing_class"] = ""

        if g.shape[0]>0:
            for gene,csq,bi,af,cn,macn,germline,subcl_likelihood,timing_class in  g.values:
                allelic_status="biallelic" if bi else "monoallelic"
                germline_status="germline"+germline

                if np.isfinite(subcl_likelihood):
                    clonality_purple="clonal" if subcl_likelihood < 0.80 else "subclonal"
                else:
                    clonality_purple = "clonal"
                clonality_timing=timing_class
                result.append(";".join([gene,csq,allelic_status,germline_status,str(clonality_purple),str(clonality_timing),"{:.2f}".format(af),"{:.2f}".format(cn),"{:.2f}".format(macn)]))
    return "___".join(result)

def valid_csq(csqns):
    valid = []
    for csq in csqns:
        found=False
        for option in coding_csq:
            if option in csq:
                found=True
        valid.append(found)
    return valid

def get_category(row):
    if ";" in str(row["info_germline"]):
        if row["lof"] == True:
            return "_lof"
        else:
            return "_aachange"
    else:
        return "_wt"

def get_loh_lilac(df):
    v=False
    for x in ["A","B","C"]:
        v  = v or  ((np.abs(df[f"{x}1_CN"].values[0] - df[f"{x}2_CN"].values[0]) > 0.4) and (df[f"{x}2_CN"].values[0] < 0.3 or df[f"{x}1_CN"].values[0] < 0.3))
    return v

def get_imbalance_lilac(df):
    v=False
    for x in ["A","B","C"]:
        v  = v or  (((df[f"{x}1_CN"].values[0] / df[f"{x}2_CN"].values[0]) >= 2.0) and (df[f"{x}2_CN"].values[0] > 0.7)) or (((df[f"{x}2_CN"].values[0] / df[f"{x}1_CN"].values[0]) >= 2.0) and (df[f"{x}1_CN"].values[0] > 0.7))
    return v

def get_mut_lilac(df):
    v=False
    for x in ["A1_MUT","A2_MUT","B1_MUT","B2_MUT","C1_MUT","C2_MUT"]:
        v  = v or  (df[x].values[0] > 0.0)
    return v


def load_data(path_base,timing_file,sample,path_lilac):
    # load somatic mutations
    somatic_muts = os.path.join(path_base, "somatic_mutations_APP.tsv.gz")
    df_somatic_muts = pd.read_csv(somatic_muts, sep="\t")
    df_somatic_muts["POS"] = df_somatic_muts["POS"].astype(int)
    df_somatic_muts["CHROM"] = df_somatic_muts["CHROM"].astype(str)
    # load somatic mutations with timing
    try:
        df_timing = pd.read_csv(timing_file,sep="\t")
        df_timing["POS"] = df_timing["POS"].astype(int)
        df_timing["CHROM"] = df_timing["CHROM"].astype(str)
        df_somatic_muts = df_somatic_muts.merge(df_timing.rename(columns={"ALT": "ALT_1"}), how="left")
    except:
        pass # ignore it


    # load germline mutations in genes with somatic variants
    somatic_muts_genes = df_somatic_muts["gene"].unique()
    file_germline = os.path.join(path_base, "germline_mutations_APP.tsv.gz")

    if len(somatic_muts_genes) > 0 and os.path.exists(file_germline):
        df_g = pd.read_csv(file_germline,sep="\t")
        df_g = df_g[(df_g["gene"].isin(somatic_muts_genes))&(valid_csq(df_g["csq"]))]
        df_somatic_muts["germline_variant"] = "_wt"
        if df_g.shape[0] > 0:
            # if there are mutations overlapping...
            df_somatic_muts["transcript"] = df_somatic_muts.apply(lambda row: row["ANN"].split("|")[6], axis=1)
            df_g["transcript"] = df_g.apply(lambda row: row["ANN"].split("|")[6] if not (":" in row["ANN"].split("|")[6]) else row["ANN"].split("|")[6].split(":")[2], axis=1)

            df_g["info_germline"] = df_g.apply(lambda row:  row["gene"] + ";" + row["ID"] + ";" + row["ANN"],axis=1)
            df_g["lof"] = df_g.apply(lambda row: row["csq"] in ["stop_gain","start_lost","frameshift_variant","stop_lost"],axis=1)
            df_somatic_muts=df_somatic_muts.merge(df_g[["transcript","gene","info_germline","lof"]].drop_duplicates(),how="left")
        if "info_germline" in df_somatic_muts.columns.values:
            df_somatic_muts["germline_variant"] = df_somatic_muts.apply(lambda row:  get_category(row),axis=1 )
        else:
            df_somatic_muts["germline_variant"] = "_wt"
    else:
        df_somatic_muts["germline_variant"] = "_wt"

    dels = os.path.join(path_base, "DEL_genes.tsv.gz")
    df_dels = pd.read_csv(dels, sep="\t")
    amps = os.path.join(path_base, "AMP_genes.tsv.gz")
    df_amps = pd.read_csv(amps, sep="\t")
    loh = os.path.join(path_base, "LOHA_HLA.tsv.gz")
    df_loh = pd.read_csv(loh, sep="\t").drop_duplicates()
    # set integer allele ploidy
    df_loh["integer_minor"] = df_loh.apply(lambda row: round(row["minMinorAlleleCopyNumber"]),axis=1)
    df_loh["minMajorAlleleCopyNumber"] = df_loh.apply(lambda row: row["minCopyNumber"] - row["minMinorAlleleCopyNumber"], axis=1)
    df_loh["integer_major"] = df_loh.apply(lambda row: round(row["minMajorAlleleCopyNumber"]), axis=1)
    purity = os.path.join(path_base, "purity.tsv.gz")
    df_p = pd.read_csv(purity, sep="\t")
    # HLA
    muts_hla = os.path.join(path_base, "somatic_mutations_HLA.tsv.gz")
    try:
        df_muts_hla = pd.read_csv(muts_hla, sep="\t")
    except FileNotFoundError:
        df_muts_hla = pd.DataFrame([])
    try:
        indels_hla = os.path.join(path_base, "somatic_indels_HLA.tsv.gz")
        df_indels_hla = pd.read_csv(indels_hla, sep="\t").rename(columns={"individual":"sampleID","CHROM":"o"})
        df_indels_hla.rename(columns={"o":"individual"},inplace=True)
    except FileNotFoundError:
        df_indels_hla = pd.DataFrame([])
    # path lilac
    df_lilac = pd.read_csv(os.path.join(path_lilac,sample+".lilac.csv.gz"))
    df_lilac_qc = pd.read_csv(os.path.join(path_lilac,sample+".lilac.qc.csv.gz"))

    df_lilac["Allele_number"] = ["A1", "A2", "B1", "B2", "C1", "C2"]
    df_lilac["Allele_number_loh"] = list(
        [x + y for x, y in itertools.product(["A1", "A2", "B1", "B2", "C1", "C2"], ["_CN"])])
    df_lilac["nonsynonymous_mut"] = df_lilac["SomaticMissense"] + df_lilac["SomaticNonsenseOrFrameshift"] + df_lilac[
        "SomaticSplice"] + df_lilac["SomaticInframeIndel"]
    df_lilac["Allele_number_mut"] = list(
        [x + y for x, y in itertools.product(["A1", "A2", "B1", "B2", "C1", "C2"], ["_MUT"])])
    df_lilac["sample_id"] = sample



    return df_somatic_muts, df_dels, df_amps, df_loh, df_p, df_muts_hla, df_indels_hla, df_lilac, df_lilac_qc

def check_mutations_hla(df_muts,df_indels):
    r = []
    if df_muts.shape[0] == 0 and df_indels.shape[0] ==0:
        return ["","",""]
    for hla in mhc_I:
        text = ""
        if "gene" in df_muts.columns.values and hla in df_muts["gene"].unique():
            df_muts["synonymous"] = df_muts.apply(lambda row: "synonymous" if row["protein_change"][2] == row["protein_change"][-1] else "non-synonymous",axis=1 )
            text+=";".join(list(df_muts[(df_muts["gene"]==hla)][["gene","individual","protein_change","synonymous"]].values[0]))
        if df_indels.shape[0] >0:
            if hla in df_indels["gene"].unique():
                df_indels["synonymous"] = df_indels.apply(
                    lambda row: "synonymous" if row["protein_change"][2] == row["protein_change"][
                        -1] else "non-synonymous", axis=1)
                if text != "":
                    text+="___"

                text+=";".join(list(df_indels[(df_indels["gene"]==hla)][["gene","individual","protein_change","synonymous"]].values[0]))
        r.append(text)
    return r

def get_count_timing(df_total,timing):

    categories_timing = ["clonal [early]", "clonal [late]", "clonal [NA]", "subclonal"]
    if os.path.exists(timing):
        try:
            df_timing = pd.read_csv(timing,sep="\t")
        except:
            for cat in categories_timing:
                df_total[cat] = 0
            df_total["timing_muts"]=0
            return df_total

        if df_timing.shape[0] == 0:
            for cat in categories_timing:
                df_total[cat] = 0
            df_total["timing_muts"] = 0
            return df_total
        for cat in categories_timing:
            df_total[cat] = df_timing[df_timing["timing_class"]==cat].shape[0]
        df_total["timing_muts"] = df_timing.shape[0]
    else:
        for cat in categories_timing:
            df_total[cat] = 0
        df_total["timing_muts"] = 0

    return df_total



def get_diversity(path_diversity):
    if not(os.path.exists(path_diversity)):
        return np.nan,np.nan,np.nan,np.nan
    df = pd.read_csv(path_diversity,sep="\t",usecols=["ID","AlleleNumber","Divergence_Average","Divergence_Sum","AlleleList"])
    print (df.head())
    total_n,total_avg,total_sum = df[df["ID"]=="germline_diversity"][["AlleleNumber","Divergence_Average","Divergence_Sum"]].values.tolist()[0]
    mean_locus=df[df["ID"].str.contains("HLA")]["Divergence_Average"].mean()
    return total_n, total_avg,total_sum, mean_locus



@click.command()
@click.option('--input_dir',
              type=click.Path(exists=True),
              help="Input dir with summary data",
              required=True)
@click.option('--timing',
              type=click.Path(),
              help="Timing file",
              required=True) # This is only for internal use, not part of the manuscript
@click.option('--sample',
              help="Sample id",
              required=True)
@click.option('--path_fusions',
              type=click.Path(),
              help="Path with linx fusions",
              required=False)
@click.option('--path_lilac',
              type=click.Path(),
              help="Path base to lilac",
              required=False)
@click.option('--path_diversity',
              type=click.Path(),
              help="Path germline diversity calculation",
              required=False)
@click.option('--output_file',
              type=click.Path(),
              help="Output file)",
              required=True)

def run(input_dir,timing,sample,path_fusions,path_lilac,path_diversity,output_file):

    # load data
    df_somatic_muts, df_dels, df_amps, df_loh, df_p, df_muts_hla, df_indels_hla, df_lilac, df_lilac_qc  =load_data(input_dir,timing,sample,path_lilac)
    # start the processing
    list_results = [sample]
    list_names = ["sample"]
    # HLA-I locus alterations
    ## high level deletions
    dels_mhc_I = check_cnvs(df_dels, mhc_I)
    list_results += dels_mhc_I
    list_names += ["del_" + m for m in mhc_I]
    dels_mhc_I_other = check_cnvs_condensed(df_dels, mhc_I_other)
    list_results += [dels_mhc_I_other]
    list_names += ["del_MHC_I_other"]
    dels_mhc_II = check_cnvs_condensed(df_dels, mhc_II)
    list_results += [dels_mhc_II]
    list_names += ["del_MHC-II"]
    # LOH of hla
    loh_hla = check_loh(df_loh)
    list_results += loh_hla
    list_names += ["loh_" + m for m in mhc_I]
    loh_hla = check_loh_imbalance(df_loh)
    list_results += loh_hla
    list_names += ["loh_imbalance_" + m for m in mhc_I]
    # check ploidys
    loh_minor_p = check_ploidy(df_loh,column_cn="integer_minor")
    list_results += loh_minor_p
    list_names += ["ploidy_minor_" + m for m in mhc_I]
    loh_major_p = check_ploidy(df_loh, column_cn="integer_major")
    list_results += loh_major_p
    list_names += ["ploidy_major_" + m for m in mhc_I]

    # check events loh and imbalance
    loh_purple = is_loh_purple(df_loh)
    list_results += [loh_purple]
    list_names += ["loh_purple"]
    imbalance_purple = is_imbalance_purple(df_loh)
    list_results += [imbalance_purple]
    list_names += ["imbalance_purple"]
    loh_lohhla = is_loh_lohhla(df_loh)
    list_results += [loh_lohhla]
    list_names += ["loh_lohhla"]
    # somatic mutations of hla locus

    list_names += ["mut_" + m for m in mhc_I]
    list_results += check_mutations_hla(df_muts_hla,df_indels_hla)

    # Somatic muts
    muts_transport_mhc = check_mutations(df_somatic_muts, transport_mhc)
    list_results += [muts_transport_mhc]
    list_names += ["mut_transport_mhc"]
    muts_scaffold_mhc = check_mutations(df_somatic_muts, scaffold_mhc)
    list_results += [muts_scaffold_mhc]
    list_names += ["mut_scaffold_mhc"]
    muts_endopeptidases = check_mutations(df_somatic_muts, endopeptidases)
    list_results += [muts_endopeptidases]
    list_names += ["mut_endopeptidases"]
    muts_proteasome = check_mutations(df_somatic_muts, proteasome)
    list_results += [muts_proteasome]
    list_names += ["mut_proteasome"]
    muts_immunoproteasome = check_mutations(df_somatic_muts, immunoproteasome)
    list_results += [muts_immunoproteasome]
    list_names += ["mut_immunoproteasome"]
    muts_chaperones = check_mutations(df_somatic_muts, chaperones)
    list_results += [muts_chaperones]
    list_names += ["mut_chaperones"]
    muts_nmd = check_mutations(df_somatic_muts, nmd)
    list_results += [muts_nmd]
    list_names += ["mut_nmd"]
    muts_interferon = check_mutations(df_somatic_muts, interferon)
    list_results += [muts_interferon]
    list_names += ["mut_interferon"]
    muts_tfs = check_mutations(df_somatic_muts, tfs)
    list_results += [muts_tfs]
    list_names += ["mut_tfs"]
    muts_elongation = check_mutations(df_somatic_muts, elongation)
    list_results += [muts_elongation]
    list_names += ["mut_elongation"]
    muts_stop_readthrough = check_mutations(df_somatic_muts, stop_readthrough)
    list_results += [muts_stop_readthrough]
    list_names += ["mut_stop_readthrough"]
    muts_cd58 = check_mutations(df_somatic_muts, cd58_nk)
    list_results += [muts_cd58]
    list_names += ["mut_cd58"]
    # deletions
    dels_transport_mhc = check_cnvs_condensed(df_dels, transport_mhc)
    list_results += [dels_transport_mhc]
    list_names += ["del_transport_mhc"]
    dels_scaffold_mhc = check_cnvs_condensed(df_dels, scaffold_mhc)
    list_results += [dels_scaffold_mhc]
    list_names += ["del_scaffold_mhc"]
    dels_proteasome = check_cnvs_condensed(df_dels, proteasome)
    list_results += [dels_proteasome]
    list_names += ["del_proteasome"]
    dels_immunoproteasome = check_cnvs_condensed(df_dels, immunoproteasome)
    list_results += [dels_immunoproteasome]
    list_names += ["del_immunoproteasome"]
    dels_nmd = check_cnvs_condensed(df_dels, nmd)
    list_results += [dels_nmd]
    list_names += ["del_nmd"]
    dels_interferon = check_cnvs_condensed(df_dels, interferon)
    list_results += [dels_interferon]
    list_names += ["del_interferon"]
    dels_tfs = check_cnvs_condensed(df_dels, tfs)
    list_results += [dels_tfs]
    list_names += ["del_tfs"]
    dels_stop_readthrough = check_cnvs_condensed(df_dels, stop_readthrough)
    list_results += [dels_stop_readthrough]
    list_names += ["del_stop_readthrough"]
    dels_cd58 = check_cnvs_condensed(df_dels, cd58_nk)
    list_results += [dels_cd58]
    list_names += ["del_cd58"]


    # Amplifications
    amps_endopeptidases=check_cnvs_condensed(df_amps,endopeptidases)
    list_results+=[amps_endopeptidases]
    list_names+=["amps_endopeptidases"]

    amps_inhibitors = check_cnvs_condensed(df_amps, inhibitors)
    list_results += [amps_inhibitors]
    list_names += ["amps_inhibitors"]

    epi_inhibitors = check_cnvs_condensed(df_amps, epigenetic_regulators)
    list_results += [epi_inhibitors]
    list_names += ["epigenetic_regulators"]


    # Calculate diversity/divergence
    total_alleles,avg_diversity_germline,sum_diversity_germline,avg_diversity_locus = get_diversity(path_diversity)
    list_results += [total_alleles,avg_diversity_germline,sum_diversity_germline,avg_diversity_locus]
    list_names += ["n_germline_alleles","avg_divergence_germline","sum_diversity_germline","avg_diversity_locus"]



    # add purity

    list_results+=list(df_p.values[0])
    list_names+=list(df_p.columns.values)

    # add lilac information


    a = df_lilac.pivot(values=["Allele"], columns=["Allele_number"], index=["sample_id"])
    a.columns = [g[1] for g in a.columns.values]
    b = df_lilac.pivot(values=["TumorCopyNumber"], columns=["Allele_number_loh"], index=["sample_id"])
    b.columns = [g[1] for g in b.columns.values ]
    c = df_lilac.pivot(values=["nonsynonymous_mut"], columns=["Allele_number_mut"], index=["sample_id"])
    c.columns = [g[1] for g in c.columns.values]
    loh_lilac= get_loh_lilac(b)
    mut_lilac = get_mut_lilac(c)
    imbalance_lilac= get_imbalance_lilac(b)

    values = list(a.values[0]) + list(b.values[0]) + list(c.values[0]) + list([loh_lilac,mut_lilac,imbalance_lilac]) + list(df_lilac_qc[["Status", "HlaYAllele"]].values[0])
    cols = list(a.columns.values) + list(b.columns.values) + list(c.columns.values) + ["loh_lilac","mut_hla_lilac","imbalance_lilac","QC_LILAC","HLA_Y"]
    list_results+=values
    list_names+=cols

    loh_focal, loh_non_focal, loh_hfocal = "", "", ""
    if loh_lilac:
        path_cnv = f"/hpc/cuppen/shared_resources/HMF_data/DR-104-update4/somatics/{sample}/purple/{sample}.purple.cnv.gene.tsv" # this needs to adjusted to your location
        if sample.startswith("DO"): # PCAWG
            path_cnv = f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{sample}-from-jar/purplesoft3.3/{sample}T.purple.cnv.gene.tsv" # this needs to adjusted to your location
        cnv = pd.read_csv(path_cnv,sep="\t")
        cnv = cnv[cnv["gene"].isin(["HLA-A","HLA-B","HLA-C"])][["gene","minRegionStart","minRegionEnd"]]
        cnv["len"] = cnv.apply(lambda row: row["minRegionEnd"] - row["minRegionStart"],axis=1)
        # if hfocal
        cnv["is_hfocal"] = cnv.apply(lambda row: row["len"] < 3*10**6,axis=1)
        cnv["is_focal"] = cnv.apply(lambda row: row["len"] < (61*10**6)*0.75,axis=1) # chr6p lenght  61,000,000, 75% of chr6p lenght
        cnv["is_nonfocal"] = cnv.apply(lambda row: row["len"] >= (61*10**6)*0.75,axis=1) # chr6p lenght  61,000,000

        loh_hfocal=",".join(list(cnv[cnv["is_hfocal"]]["gene"].values))
        loh_focal=",".join(list(cnv[cnv["is_focal"]]["gene"].values))
        loh_non_focal=",".join(list(cnv[cnv["is_nonfocal"]]["gene"].values))

    imbalance_hfocal, imbalance_focal, imbalance_non_focal = "", "", ""
    if imbalance_lilac:
        path_cnv = f"/hpc/cuppen/shared_resources/HMF_data/DR-104-update4/somatics/{sample}/purple/{sample}.purple.cnv.gene.tsv" # this needs to adjusted to your location
        if sample.startswith("DO"): # PCAWG
            path_cnv = f"/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/{sample}-from-jar/purplesoft3.3/{sample}T.purple.cnv.gene.tsv" # this needs to adjusted to your location
        cnv = pd.read_csv(path_cnv,sep="\t")
        cnv = cnv[cnv["gene"].isin(["HLA-A","HLA-B","HLA-C"])][["gene","minRegionStart","minRegionEnd"]]
        cnv["len"] = cnv.apply(lambda row: row["minRegionEnd"] - row["minRegionStart"],axis=1)
        # if hfocal
        cnv["is_hfocal"] = cnv.apply(lambda row: row["len"] < 3*10**6,axis=1)
        cnv["is_focal"] = cnv.apply(lambda row: row["len"] < (61*10**6)*0.75,axis=1) # chr6p lenght  61,000,000, 75% of chr6p lenght
        cnv["is_nonfocal"] = cnv.apply(lambda row: row["len"] >= (61*10**6)*0.75,axis=1) # chr6p lenght  61,000,000

        imbalance_hfocal=",".join(list(cnv[cnv["is_hfocal"]]["gene"].values))
        imbalance_focal=",".join(list(cnv[cnv["is_focal"]]["gene"].values))
        imbalance_non_focal=",".join(list(cnv[cnv["is_nonfocal"]]["gene"].values))

    values = [loh_hfocal,loh_focal,loh_non_focal,imbalance_hfocal,imbalance_focal,imbalance_non_focal]
    cols = ["loh_hfocal","loh_focal","loh_nonfocal","imbalance_hfocal","imbalance_focal","imbalance_nonfocal"]

    list_results+=values
    list_names+=cols




    df_total = pd.DataFrame([list_results],columns=list_names)
    print (path_fusions)
    # Finally include number of fusions
    column_reported = "Reported"
    if os.path.exists(path_fusions):
        df_fus_total = pd.read_csv(path_fusions,sep="\t")
        if not (column_reported) in df_fus_total.columns.values:
            column_reported="reported"
        total_fusions = df_fus_total.shape[0]
        total_reported_fusions = df_fus_total[df_fus_total[column_reported]].shape[0]
        df_total["total_fusions"] = total_fusions
        df_total["total_reported_fusions"] = total_reported_fusions


    df_total = get_count_timing(df_total,timing)




    df_total.to_csv(output_file,sep="\t",index=False)

if __name__ == '__main__':
    run()