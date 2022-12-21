import os, sys
import pandas as pd
import allel
import numpy as np
import click


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


total = mhc_I + mhc_I_other + mhc_II + transport_mhc + scaffold_mhc + endopeptidases + endopeptidases + proteasome + immunoproteasome + chaperones + nmd + interferon + tfs + elongation + stop_readthrough + cd58_nk
somatic_mut_genes = list(set(total) - set(mhc_I) - set(mhc_II))


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
    mnvs=count_mnvs(df_muts["REF"],df_muts["ALT_1"])
    frameshfit=df_muts[df_muts["ANN"].str.contains("frameshift", na=False)].shape[0]
    missense=df_muts[df_muts["ANN"].str.contains("missense", na=False)].shape[0]
    tmb = df_muts.shape[0]
    df_muts = df_muts[~pd.isnull(df_muts["ANN"])]
    if df_muts.shape[0] >0:
        df_muts[["csq", "gene", "impact"]] = df_muts.apply(
            lambda row: pd.Series([row["ANN"].split("|")[1], row["ANN"].split("|")[3], row["ANN"].split("|")[2]]), axis=1)
        df_muts = df_muts[
            (df_muts["gene"].isin(total)) & ((df_muts["impact"].isin(["MODERATE", "HIGH"])) & (valid_csqn(df_muts["csq"])))]
    return df_muts,tmb,missense,frameshfit,mnvs

def load_germline(path_file,region):
    muts_germline = allel.read_vcf(path_file, fields=['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT',
                                                       'variants/ANN', 'variants/FILTER_PASS', 'calldata/GT',
                                                       'calldata/AD', "variants/ID"], region=region)
    if not(muts_germline):
        return pd.DataFrame([])
    df_germline = pd.DataFrame([])
    germline_status, genotype = [], []
    # process gt
    for row in muts_germline["calldata/GT"]:
        ref = row[0]
        genotype.append("/".join([str(x) for x in ref]))
        if ref[0] == ref[1] and ref[0] == 0:
            germline_status.append("reference")
        elif ref[0] == ref[1] and ref[0] == 1:
            germline_status.append("homozygous")
        elif ref[0] != ref[1] and ref[1] == 1:
            germline_status.append("heterozygous")
        elif ref[0] != ref[1] and ref[0] == 1:
            germline_status.append("heterozygous")
        else:
            germline_status.append("other")
    df_germline["germline_status"] = germline_status
    df_germline["genotype_germline"] = genotype
    # process ad
    ad = []
    for row in muts_germline["calldata/AD"]:
        ref = row[0]
        ad.append("/".join([str(x) for x in ref]))
    df_germline["AD"] = ad
    # prepare calldata
    for key in muts_germline:
        if "calldata" in key:
            continue
        if "ALT" in key:
            df_germline[key.split("/")[1]] = [alts[0] for alts in muts_germline[key]]
        else:
            df_germline[key.split("/")[1]] = muts_germline[key]
    return df_germline

def read_germline_muts(regions,path_germline):
    germline_muts = []
    for gene, chr_, s, e in regions[~regions["gene"].str.contains("HLA")][["gene", "chr", "pos_start", "pos_end"]].values: # iterate over the regions of itnerest
        # ger the region
        r = str(chr_) + ":" + str(s) + "-" + str(e)
        # Read it, saves a lot of memory
        df_muts_germline=load_germline(path_germline,r)
        if (not (df_muts_germline is None)) and (df_muts_germline.shape[0] > 0):
            x = df_muts_germline[(df_muts_germline["FILTER_PASS"]) & (~pd.isnull(df_muts_germline["ANN"]))&(df_muts_germline["germline_status"]!="reference")]
            if x.shape[0] >0:
                x=x[(x["ANN"].str.contains(gene)) & ((x["ANN"].str.contains("HIGH")) | (x["ANN"].str.contains("MODERATE")))]
            if x.shape[0] > 0:
                germline_muts.append(x)
    return pd.concat(germline_muts)


def calculate_caf(row,purity):
    ref_counts = int(row["t_ref_count"])
    alt_coutns = int(row["t_alt_count"])
    lp = float(row["maxCopyNumber"])
    purity = float(purity)
    # CAF = VAF * (Lp * Pt + 2 * (1 - Pt)) / (Lp * Pt)
    vaf = ref_counts / (ref_counts + alt_coutns)
    caf = vaf * (lp*purity + 2 * (1 - purity)) / (lp * purity)
    return float("{:.4f}".format(caf))


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
@click.option('--germline_vcf',
              type=click.Path(exists=True),
              help="Input path of vcf with germline mutations",
              required=False)
@click.option('--somatic_hla_mutations',
              type=click.Path(),
              help="Input path of tabulated file with somatic mutations in HLA locus",
              required=False)
@click.option('--somatic_hla_indels',
              type=click.Path(),
              help="Input path of tabulated file with somatic indels in HLA locus",
              required=False)
@click.option('--loh_hla',
              type=click.Path(exists=True),
              help="Input path of tabulated file with LOH events in HLA locus",
              required=False)
@click.option('--outpath',
              type=click.Path(),
              help="Output path (directory)",
              required=True)


def run(somatic_vcf, somatic_cnv, purity_info, germline_vcf, somatic_hla_mutations, somatic_hla_indels, loh_hla, outpath):

    # create output
    if not(os.path.exists(outpath)):
        os.mkdir(outpath)

    # Read sample ploidy and info
    df_purity = pd.read_csv(purity_info,  sep="\t")
    df_purity=df_purity[["purity", "ploidy", "diploidProportion", "gender", "wholeGenomeDuplication", "msStatus", "msIndelsPerMb","tml", "tmlStatus","tmbPerMb", "tmbStatus", "svTumorMutationalBurden"]]

    sample_ploidy = float(df_purity["ploidy"])


    # Read somatic cnv
    df_cnv = pd.read_csv(somatic_cnv,sep="\t")
    df_cnv = df_cnv[df_cnv["gene"].isin(total+inhibitors+epigenetic_regulators)]
    ## Determine AMP and DEL events based on HMF driver calling rules
    threshold = float(3 * sample_ploidy)
    df_amp = df_cnv[(df_cnv["minCopyNumber"] > threshold)]
    df_amp.to_csv(os.path.join(outpath, "AMP_genes.tsv.gz"), sep="\t", compression="gzip", index=False)
    df_loss = df_cnv[(df_cnv["minCopyNumber"] < 0.5)]
    df_loss.to_csv(os.path.join(outpath, "DEL_genes.tsv.gz"), sep="\t",compression="gzip", index=False)


    # Read somatic mutations
    df_muts,TMB,missense,frameshift,mnvs = load_somatic_mutations(somatic_vcf)


    df_purity["TMB"] = TMB
    df_purity["missense"] = missense
    df_purity["frameshift"] = frameshift
    df_purity["mnvs"] = mnvs
    df_purity.to_csv(os.path.join(outpath, "purity.tsv.gz"), sep="\t", compression="gzip", index=False)
    ## Include CN info, for the LOH
    df_muts = df_cnv[["gene","minCopyNumber","maxCopyNumber","minMinorAlleleCopyNumber"]].merge(df_muts,how="right")
    df_muts.to_csv(os.path.join(outpath, "somatic_mutations_APP.tsv.gz"), sep="\t", compression="gzip", index=False)

    # Read somatic_hla_mutations
    if somatic_hla_mutations and os.path.exists(somatic_hla_mutations):
        df_muts_hla = pd.read_csv(somatic_hla_mutations,sep="\t")
        if df_muts_hla.shape[0] >0:
            df_muts_hla["gene"] = df_muts_hla.apply(lambda row: str(row["individual"])[0:5].upper().replace("_", "-"), axis=1)
            df_muts_hla = df_cnv[["gene", "minCopyNumber", "maxCopyNumber", "minMinorAlleleCopyNumber"]].merge(df_muts_hla, how="right")
            #CAF = VAF * (Lp * Pt + 2 * (1 - Pt)) / (Lp * Pt)
            df_muts_hla["CAF"] = df_muts_hla.apply(lambda row: calculate_caf(row,df_purity["purity"]),axis=1)
        df_muts_hla.to_csv(os.path.join(outpath, "somatic_mutations_HLA.tsv.gz"), sep="\t", compression="gzip", index=False)
    # Read somatic_hla_mutations
    if somatic_hla_indels and os.path.exists(somatic_hla_indels):
        df_muts_hla = pd.read_csv(somatic_hla_indels, sep="\t")
        if df_muts_hla.shape[0] >0:
            df_muts_hla["gene"] = df_muts_hla.apply(lambda row: str(row["CHROM"])[0:5].upper().replace("_", "-"),axis=1)
        df_muts_hla.to_csv(os.path.join(outpath, "somatic_indels_HLA.tsv.gz"), sep="\t", compression="gzip",
                                index=False)
    # Read LOHHLA HLA
    if loh_hla and os.path.exists(loh_hla):
        error = False
        df_c = df_cnv[df_cnv["gene"].isin(["HLA-A","HLA-B","HLA-C"])][["gene", "minCopyNumber", "maxCopyNumber", "minMinorAlleleCopyNumber"]]
        try:
            df_lohhla = pd.read_csv(loh_hla,sep="\t")
        except:
            error = True
        if (not(error)) and df_lohhla.shape[0] > 0:
            df_lohhla = df_lohhla[df_lohhla["PVal_unique"] < 0.05][["LossAllele", "KeptAllele", "PVal_unique"]]
            if df_lohhla.shape[0] >0:
                df_lohhla["gene"] = df_lohhla.apply(lambda row: str(row["LossAllele"])[0:5].upper().replace("_", "-"),axis=1)
                df_lohhla = df_c.merge(df_lohhla, how="left")
                df_lohhla.to_csv(os.path.join(outpath, "LOHA_HLA.tsv.gz"), sep="\t", compression="gzip",
                                 index=False)
            else:
                df_lohhla = df_c
                df_lohhla.to_csv(os.path.join(outpath, "LOHA_HLA.tsv.gz"), sep="\t", compression="gzip",
                                 index=False)
        else:
            df_lohhla=df_c
            df_lohhla.to_csv(os.path.join(outpath, "LOHA_HLA.tsv.gz"), sep="\t", compression="gzip",
                                 index=False)
    else:
        df_cnv[df_cnv["gene"].isin(["HLA-A", "HLA-B", "HLA-C"])][
            ["gene", "minCopyNumber", "maxCopyNumber", "minMinorAlleleCopyNumber"]].to_csv(os.path.join(outpath, "LOHA_HLA.tsv.gz"), sep="\t", compression="gzip",
                         index=False)


    # Finally, read germline mutations, neets to be annotated!
    if germline_vcf and False:
        # Load the regions file
        df_regions = pd.read_csv("/hpc/local/CentOS7/cog/software/polysolver/data/all_genes.37.tsv.gz",sep="\t",usecols=list(range(0, 5)),
                                 names=["chr", "pos_start", "pos_end", "ensembl_gene", "gene"])
        regions = df_regions[df_regions["gene"].isin(total)].groupby(["gene", "chr", "ensembl_gene"],                                                                  as_index=False).agg(
            {"pos_start": np.nanmin, "pos_end": np.nanmax}).sort_values(["chr", "pos_start", "pos_end"])

        # Read germline mutations
        df_germline_muts = read_germline_muts(regions,germline_vcf)
        df_germline_muts[["csq", "gene", "impact"]] = df_germline_muts.apply(
            lambda row: pd.Series([row["ANN"].split("|")[1], row["ANN"].split("|")[3], row["ANN"].split("|")[2]]),
            axis=1)
        df_germline_muts.to_csv(os.path.join(outpath, "germline_mutations_APP.tsv.gz"), sep="\t", compression="gzip",
                           index=False)





if __name__ == '__main__':
    run()