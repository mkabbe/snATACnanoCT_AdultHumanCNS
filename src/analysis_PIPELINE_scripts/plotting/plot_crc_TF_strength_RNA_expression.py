import os
import csv
import scipy
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import style
import warnings
style.use("default")
warnings.filterwarnings('ignore')
#%config InlineBackend.figure_format = "retina"
os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE/")


## Plot GEX and TF strength (ranked plot)

rnadata = ad.read("../../CZI/episcanpy_analysis/AGG_ATAC_210218/data/feature_matrices/HCA_RNA_all_annot.h5ad")
rnadata.var_names = rnadata.var.gene_names
rnadata = rnadata[rnadata.obs.Tissue!="CB"]

cluster_dict = {"AST":['Astrocyte_1','Astrocyte_2','Astrocyte_3','Astrocyte_4'],
               "MOL":["Oligo"],
               "OPC":["OPC"],
               "MIGL":["Microglia-Macrophages_1","Microglia-Macrophages_2"],
               "CXEX":["Neuron_Ex_1","Neuron_Ex_2","Neuron_Ex_3"],
               "CXINH":["Neuron_In_1","Neuron_In_2"]}


fig,ax = plt.subplots(1,6,figsize=(24,5))
for i,c1 in enumerate(["MOL","OPC_ATAC","MIGL","AST","CXEX","CXINH_ATAC"]):
    if c1=="OPC_ATAC":
        c = "OPC"
    elif c1=="CXINH_ATAC":
        c ="CXINH"
    else:
        c = c1

    df = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/{c1}_roseSE/crc_out_3/{c1}_CRC_DEGREE_TABLE.txt",sep="\t")
    df[f"TF_strength_{c1}"] = df[f"In_Degree"] - df[f"Out_Degree"]
    df = df.sort_values(by=f"TF_strength_{c1}",ascending=False)
    df["TF_rank"] = [x for x in range(1,len(df)+1)]
    df["isHox"] = [True if x.startswith("HOX") else False for x in df.Tf]

    nodes = pd.read_csv(f"/date/gcb/GCB_MK/CRC_output/{c}_CRC_NODELIST.txt",header=None)
    node_list = [x for x in nodes[0] if x in set(rnadata.var_names)]

    subset = rnadata[rnadata.obs.clusters_named.isin(cluster_dict[c]),rnadata.var_names.isin(node_list)]
    node_list = [x for x in subset.var.gene_names] # preserve order
    df2 = pd.DataFrame()
    df2["nodes"] = node_list
    df2[f"{c}_mean_gex"] = list(np.mean(subset.X.toarray(),axis=0))
    gex_dict = dict(zip(df2.nodes,df2[f"{c}_mean_gex"]))
    df["mean_gex"] = [gex_dict[x] if x in gex_dict else np.nan for x in df.Tf]

    sns.scatterplot(data= df.sort_values(by="mean_gex",ascending=True), x="TF_rank", y=f"TF_strength_{c1}",
                    palette="Reds", size="mean_gex", sizes=(50, 500),
                    ec="None",ax=ax[i], alpha=0.5, hue=f"mean_gex")
    
    
    df = df.sort_values(by="mean_gex",ascending=False)
    df["GEX_rank"] = [_ for _ in range(1,len(df)+1)]
    df = df.sort_values(by="TF_rank",ascending=True)
    df = df.reset_index()

    for j in range(0,5):
        ax[i].annotate(df["Tf"][j], (df["TF_rank"][j], df[f"TF_strength_{c1}"][j]),fontsize=8)


    sns.despine()
    ax[i].legend("",frameon=False);
    ax[i].axhline(0,lw=1, color="black",linestyle="--")
    ax[i].set_ylabel(f"TF_strength_{c}")

plt.tight_layout()

plt.savefig(f"plots/CRC_TF_ntx_ranking_by_GEX_scatter_ALL.png", dpi=400, bbox_inches="tight")



## Plot GEX by TF strength scatter


fig,ax = plt.subplots(2,2,figsize=(8,8))
ax=ax.flatten()
    
for i,c1 in enumerate(["MOL","OPC_ATAC"]):
    
    if c1=="OPC_ATAC":
        c = "OPC"
    elif c1=="CXINH_ATAC":
        c ="CXINH"
    else:
        c = c1

    df = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/{c1}_roseSE/crc_out_3/{c1}_CRC_DEGREE_TABLE.txt",sep="\t")
    df[f"TF_strength_{c1}"] = df[f"In_Degree"] - df[f"Out_Degree"]
    df = df.sort_values(by=f"TF_strength_{c1}",ascending=False)
    df["TF_rank"] = [x for x in range(1,len(df)+1)]
    df["isHox"] = [True if x.startswith("HOX") else False for x in df.Tf]

    nodes = pd.read_csv(f"/date/gcb/GCB_MK/CRC_output/{c}_CRC_NODELIST.txt",header=None)
    node_list = [x for x in nodes[0] if x in set(rnadata.var_names)]

    subset = rnadata[rnadata.obs.clusters_named.isin(cluster_dict[c]),rnadata.var_names.isin(node_list)]
    node_list = [x for x in subset.var.gene_names] # preserve order


    df2 = pd.DataFrame()
    df2["nodes"] = node_list
    df2[f"{c}_mean_gex"] = list(np.mean(subset.X.toarray(),axis=0))

    gex_dict = dict(zip(df2.nodes,df2[f"{c}_mean_gex"]))

    df["mean_gex"] = [gex_dict[x] if x in gex_dict else 0 for x in df.Tf]
    df["wTF_strength"] = df[f"TF_strength_{c1}"] * df["mean_gex"]

    df.sort_values(by="mean_gex",ascending=False,inplace=True)
    df["GEX_rank"] = [x for x in range(1,len(df)+1)]


    
    sns.scatterplot(x=df[f"TF_strength_{c1}"],y=df["mean_gex"], ec=None,alpha=0.5,hue=df["mean_gex"],palette="viridis",s=150,ax=ax[i])
    sns.scatterplot(data= df.sort_values(by="isHox"), x=f"TF_strength_{c1}", y=f"mean_gex",palette=["gainsboro","red"],s=150, ec="None", alpha=0.5, hue="isHox",ax=ax[i+2])

    ax[i].legend("",frameon=False);
    ax[i].axvline(1,lw=3, linestyle="--",color="black",alpha=0.5)
    ax[i].set_xlabel(f"TF_strength_{c}")
    ax[i].set_ylabel("Average GEX (logNorm)")
    ax[i+2].legend("",frameon=False);
    ax[i+2].axvline(1,lw=3, linestyle="--",color="black",alpha=0.5)
    ax[i+2].set_xlabel(f"TF_strength_{c}")
    ax[i+2].set_ylabel("Average GEX (logNorm)")

sns.despine()
fig.tight_layout()

plt.savefig(f"plots/CRC_TF_ntx_ranking_by_GEX_scatter_HOX_OPC_MOL.png", dpi=400, bbox_inches="tight")
