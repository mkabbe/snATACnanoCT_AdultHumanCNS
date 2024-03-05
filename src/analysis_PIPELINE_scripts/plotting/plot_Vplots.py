import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import style
import warnings
style.use("default")
warnings.filterwarnings('ignore')
os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE")



import anndata as ad
adata = ad.read("mtx/ATAC_object_230726.h5ad")


#region = "chr2:176,140,400-176,141,100" ## HOXD4-HOXD6 TAD boundary
#region = "chr7:26,864,349-26,865,244" ## SKAP2 promoter
#region = "chr22:37,830,229-38,086,761" ##SOX10 to enhancer broad region
#region = "chr22:37,870,513-38,005,073" ## SOX10 reg
region = "chr22:37,891,500-37,895,000" ## SOX10 enh # 4kb
#region = "chr22:37,893,000-37,893,300" ## SOX10 enh zoom in # 300bp
#region = "chr3:170,412,544-170,491,247" # CLDN11
#region = "chr3:170,483,000-170,488,000" # CLDN11 enh zoom #5kb
#region = "chr3:170,485,850-170,486,000" # CLDN11 enh zoom # 150bp

#region ="chr22:37,892,132-37,894,257"
#region = "chr22:37,886,660-37,899,248"
#region = "chr22:37,892,982-37,893,318"


chrom=region.split(":")[0]
start=int(region.split(":")[-1].split("-")[0].replace(",",""))
end=int(region.split(":")[-1].split("-")[1].replace(",",""))


fig, ax= plt.subplots(2,1, figsize=[8,4], sharex=True, gridspec_kw={'height_ratios':[1,2]})
plt.subplots_adjust(hspace=None)

#ATAC MOL
df = pd.read_csv("/date/gcb/GCB_MK/ATAC_ctype_fragments_231219/chr22_mol.bed", sep="\t", header=None)
#df = pd.read_csv("/date/gcb/GCB_MK/ATAC_ctype_fragments_231219/chr3_mol.bed", sep="\t", header=None)
df["length"] = df[2] - df[1]
df["mid"] = (df[2] + df[1])/2

plot_df = df.loc[(df[1]>=start)&(df[2]<=end)&(df[0]==chrom)]

sns.kdeplot(x = plot_df[1].to_list() + plot_df[2].to_list(), color="black", lw=0.6, ax=ax[0])
density_scatter(np.array(plot_df["mid"]), np.array(plot_df["length"]),bins = [100,100], s=0.1,alpha=1, cmap="rainbow", ax=ax[1])
ax[0].set_yticks([])
ax[1].set_xlim(plot_df[1].iloc[0], plot_df[1].iloc[-1])
ax[1].set_ylim(0,300)

plt.savefig("plots/Vplot_SOX10_enhancer.png", dpi=600, bbox_inches="tight")

## Plot motif center


motif_center = (37_893_118+37_893_136)/2
chrom = "chr22"
start = motif_center - 100
end = motif_center + 100


fig, ax= plt.subplots(2,1, figsize=[8,4], sharex=True, gridspec_kw={'height_ratios':[1,2]})
plt.subplots_adjust(hspace=None)

#ATAC OLG

df = pd.read_csv("/date/gcb/GCB_MK/ATAC_ctype_fragments_231219/chr22_OLG.bed", sep="\t", header=None)
#df = pd.read_csv("/data/proj/GCB_MK/scCT/nanoCT_EAE/fragments/P29054/ctype/H3K27me3/chr22_OLG_me3.bed", sep="\t", header=None)
#df = pd.read_csv("/data/proj/GCB_MK/scCT/nanoCT_EAE/fragments/P29054/ctype/H3K27ac/chr22_OLG_ac.bed", sep="\t", header=None)

df["length"] = df[2] - df[1]
df["mid"] = (df[2] + df[1])/2

plot_df = df.loc[(df[1]>=start)&(df[2]<=end)&(df[0]==chrom)]

sns.kdeplot(x = plot_df[1].to_list() + plot_df[2].to_list(), color="black", lw=0.6, ax=ax[0])

sns.kdeplot(x=np.array(plot_df["mid"]), y=np.array(plot_df["length"]), fill=True, ax=ax[1],cmap="binary",levels=500,thresh=0)
ax[0].set_yticks([]) 
ax[0].set_ylabel("")

ax[1].set_xticks([start, motif_center-50, motif_center, motif_center+50, end],labels=["-100", "-50","0","+50","+100"])
ax[1].set_xlim(start,end)
ax[1].set_ylim(19,130)

plt.savefig("plots/TFAP2A_motif_Vplot.png", dpi=600, bbox_inches="tight")




