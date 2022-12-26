import os
import pandas as pd
import numpy as np
import random
from tqdm import tqdm
from scipy import stats
import statsmodels
import matplotlib.pyplot as plt
import seaborn as sns
import click
import pybedtools
import glob


#os.system("export QT_QPA_PLATFORM='offscreen'") # this is important in case you run it on a cluster
#os.system("export XDG_RUNTIME_DIR=/home/cog/fmartinez/tmp_plots")


def plot_ttype(res,annotations_info, i, ttype, list_chunks, output_dir_plots,name_run):

    # calculate coordinates...
    xs, ys, xticks, xticklabels, annotations, pvalue, odds_ratio,y_top,y_bottom = [], [], [], [], [], [], [], [], []
    chr_p = ""
    d = []
    for j, x in enumerate(list_chunks):
        o = res.loc[x]
        xs.append(o["n_observed"])
        ys.append(o["n_mean_simulated"])
        y_top.append(o["Q3_simulated"])
        y_bottom.append(o["Q1_simulated"])
        if x in annotations_info:
            annotations.append((j, xs[-1] + 1, i.loc[x]["gene"]))
            d.append([ttype, i.loc[x]["gene"],x, o["pvalue_emp_local"], o["pvalue_ana_global"],o["q_value_ana_global"], o["odds_ratio_global"], xs[-1], ys[-1]])
        if x.split("_")[0] != chr_p:
            xticks.append(j)
            xticklabels.append(x.split("_")[0])
            chr_p = x.split("_")[0]

    # plot ttype
    typer = name_run.split("_")[1]
    d_colors={"loh":"#a8ddb5","deepdel":"#43a2ca","amp":"#de2d26"}
    color=d_colors[typer]
    fig, ax = plt.subplots(figsize=(20, 5))
    plt.plot(xs, color=color)
    ax.fill_between(range(0, len(xs)), xs, color=color)
    plt.plot(ys, color="black")
    ax.fill_between(y_bottom, y_top, color="black",alpha=0.2)
    for (x, y, s) in annotations:
        if len(s)==2:
            s=s[0]
        ax.annotate(xy=(x, y), text=s)
        ax.axvline(x=x,ymin=0,ymax=y,ls="--",lw=0.5,color="black")
    # ax.fill_between(range(0,len(xs)),ys,color="#999999")
    ax.set_xticks(xticks)
    _ = ax.set_xticklabels(xticklabels, rotation=90)
    ax.set_title(ttype, fontsize=18)
    name ="Number of samples"
    if "loh_focal" in name_run:
        name = "Number of samples with focal allelic LOH"
    elif "loh_hfocal" in name_run:
        name = "Number of samples with highly focal allelic LOH"
    elif "loh_nonfocal" in name_run:
        name = "Number of samples with allelic LOH"

    ax.set_ylabel(name)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(y=np.nanmean(ys),xmin=0,xmax=(len(xs)),ls="--",color="black")

    plt.savefig(f'{output_dir_plots}/{ttype}_{name_run}.pdf', dpi=800, bbox_inches="tight")
    plt.close()
    dff=pd.DataFrame(d,columns=["ttype", "gene","region", "pvalue_emp_locla", "p_value_ana_global","q_value_ana_global", "odds_ratio_global", "n_obs", "n_sim"])
    dff["mean_global"] = np.nanmean(ys)
    return dff



def plot_ttype_zoom(res,annotations_info, i, ttype, list_chunks, output_dir_plots,name_run):

    # calculate coordinates...
    xs, ys, xticks, xticklabels, annotations, pvalue, odds_ratio,y_top,y_bottom = [], [], [], [], [], [], [], [], []
    d = []
    for j, x in enumerate(list_chunks):
        o = res.loc[x]
        xs.append(o["n_observed"])
        ys.append(o["n_mean_simulated"])
        y_top.append(o["Q3_simulated"])
        y_bottom.append(o["Q1_simulated"])
        if x in annotations_info:
            genes=i.loc[x]["gene"]

            if not(isinstance(genes,str)):
                genes=";".join(list(i.loc[x]["gene"]))


            annotations.append((j, xs[-1] + 1,genes ))
            #.append([ttype, i.loc[x]["gene"],x, o["pvalue_ana_global"], o["q_value_ana"], o["odds_ratio"], xs[-1], ys[-1]])

        xticks.append(j)
        if j % 10 == 0:
            xticklabels.append(x)
        else:
            xticklabels.append("")

    # plot ttype
    typer = name_run.split("_")[1]
    d_colors = {"loh": "#a8ddb5", "deepdel": "#43a2ca", "amp": "#de2d26"}
    color = d_colors[typer]
    fig, ax = plt.subplots(figsize=(20, 5))
    plt.plot(xs, color=color)
    ax.fill_between(range(0, len(xs)), xs, color=color)
    plt.plot(ys, color="black")
    ax.fill_between(y_bottom, y_top, color="black", alpha=0.2)
    for (x, y, s) in annotations:
        if len(s)==2:
            s=s[0]
        ax.annotate(xy=(x, y), s=s)
        ax.axvline(x=x,ymin=0,ymax=y,ls="--",lw=0.5,color="grey")
    # ax.fill_between(range(0,len(xs)),ys,color="#999999")
    ax.set_xticks(xticks)
    _ = ax.set_xticklabels(xticklabels, rotation=90)
    ax.set_title(ttype, fontsize=18)
    name = "Number of samples"
    if "loh_focal" in name_run:
        name = "Number of samples with focal allelic LOH "
    elif "loh_hfocal" in name_run:
        name = "Number of samples with highly focal allelic LOH"
    elif "loh_nonfocal" in name_run:
        name = "Number of samples with allelic LOH"
    ax.set_ylabel(name)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(y=np.nanmean(ys),xmin=0,xmax=(len(xs)),ls="--",color="black")

    plt.savefig(f'{output_dir_plots}/{ttype}_{name_run}_zoom.pdf', dpi=800, bbox_inches="tight")
    plt.close()





@click.command()
@click.option('--chunks_representation',
              type=click.Path(exists=True),
              help="Input file with the chunks",
              required=True)
@click.option('--input_directory',
              type=click.Path(),
              help="Input directory with files to be considered...",
              required=True)

def run(chunks_representation, input_directory):
    # path genes
    df_genes = pd.read_csv(
        "../../../external_data/selected_regions.tsv",
        sep="\t")

    # chunks
    df_chunks = pd.read_csv(chunks_representation, sep="\t")
    df_chunks["chromosome"] = df_chunks["chromosome"].astype(int)
    list_chunks = list(df_chunks.sort_values(["chromosome", "start"])["name"].values)
    list_chunks_specific = list(df_chunks[(df_chunks["chromosome"]==6)&(df_chunks["start"]>25000001)&(df_chunks["end"]<50000001)]["name"].values)
    # genes to be highlighted
    highlight_genes_general = ["SETDB1","CD274"]


    df_genes_s = df_genes[df_genes["gene"].isin(highlight_genes_general)]

    query = pybedtools.BedTool.from_dataframe(df_chunks[["chromosome", "start", "end", "name"]])
    regions = pybedtools.BedTool.from_dataframe(df_genes_s[["chromosome", "start", "end", "gene"]])
    i = query.intersect(regions, wao=True, loj=True).to_dataframe()
    i.columns = ["chr_i", "start_i_q", "end_i_q", "name", "chr_i_r", "start_i_r", "end_i_r", "gene", "length"]
    i = i[i["gene"] != "."]
    i.set_index("name", inplace=True)
    annotations_info_general = set(i[i["gene"].isin(highlight_genes_general)].index)

    # Read inputs
    list_subfolders_with_paths = [f.path for f in os.scandir(input_directory) if f.is_dir()]
    for p in list_subfolders_with_paths:
        if not("amp" in p):
            continue
        print (p)
        ttypes_info = []
        output_dir_plots = os.path.join(p,"plots")
        if not(os.path.exists(output_dir_plots)):
            os.mkdir(output_dir_plots)

        for f in tqdm(glob.glob(p + "/*.tsv.gz")):
            name_run = os.path.dirname(f).split("/")[-1]
            ttype = os.path.basename(f).split(".")[-3]
            try:
                res = pd.read_csv(f,sep="\t")
            except:
                continue
            res.set_index("region", inplace=True)
            a = plot_ttype(res,annotations_info_general, i[i["gene"].isin(highlight_genes_general)], ttype, list_chunks, output_dir_plots,name_run)
            ttypes_info.append(a)
        # create a summary dataframe
        df_summary = pd.concat(ttypes_info)
        df_summary["n_obs"] = df_summary["n_obs"].astype(int)
        df_summary.to_csv(os.path.join(output_dir_plots,f"{name_run}_summary_significance.tsv.gz"), sep="\t", index=False,compression="gzip")


if __name__ == '__main__':
    run()
