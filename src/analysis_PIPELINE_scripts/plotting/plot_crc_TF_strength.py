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



## plot OPC vs MOL TF network differences

c1 = "MOL"
c2 = "OPC"

df_1 = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/{c1}_roseSE/crc_out_3/{c1}_CRC_DEGREE_TABLE.txt",sep="\t")
df_2 = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/{c2}_roseSE/crc_out_3/{c2}_CRC_DEGREE_TABLE.txt",sep="\t")


#df_1 = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/OLG_CSC_roseSE/crc_out_3/{c1}_CRC_DEGREE_TABLE.txt",sep="\t")
#df_2 = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/OLG_BA4_roseSE/crc_out_3/{c2}_CRC_DEGREE_TABLE.txt",sep="\t")



overlap_tf = set(df_1.Tf).intersection(set(df_2.Tf))

df_1 = df_1.loc[df_1.Tf.isin(overlap_tf)]
df_2 = df_2.loc[df_2.Tf.isin(overlap_tf)]

df = pd.DataFrame(data=list(overlap_tf),columns=["Tf"]).sort_values(by="Tf")
df[f"In_{c1}"] = list(df_1["In_Degree"])
df[f"Out_{c1}"] = list(df_1["Out_Degree"])

df[f"In_{c2}"] = list(df_2["In_Degree"])
df[f"Out_{c2}"] = list(df_2["Out_Degree"])

df[f"TF_strength_{c1}"] = df[f"In_{c1}"] - df[f"Out_{c1}"]
df[f"TF_strength_{c2}"] = df[f"In_{c2}"] - df[f"Out_{c2}"]


#df["diff_SC"] = list(df_csc["Out_Degree"] - df_csc["In_Degree"]) # +ve value means TF is a strong REGULATOR
#df["diff_Ctx"] = list(df_ba["Out_Degree"] - df_ba["In_Degree"])

df["delta_Out"] = df[f"Out_{c2}"] - df[f"Out_{c1}"]
df["delta_In"] = df[f"In_{c2}"] - df[f"In_{c1}"]


lim=800
threshold = 200 #strength threshold


## Mark TFs that are REGULATORS in one, but TARGETS in the other
df["sig_diff"] = df[f'TF_strength_{c2}'] * df[f'TF_strength_{c1}'] < 0 ## Regulators vs Targets 
df["abs_diff"] = np.absolute(df[f'TF_strength_{c2}'] - df[f'TF_strength_{c1}'])
df["abs_diff_sig"] = df["abs_diff"] > threshold ## Minimum difference in strength 
df["sig"] = df["abs_diff_sig"] & df["sig_diff"] 

df["sig_class"] = "."
df.loc[(df.sig==True)&(df.TF_strength_MOL>0),"sig_class"] = "up_MOL"
df.loc[(df.sig==True)&(df.TF_strength_MOL<0),"sig_class"] = "up_OPC"


fig,ax = plt.subplots(1,2,figsize=(15,7))
ax=ax.flatten()
sns.scatterplot(x=df[f'TF_strength_{c1}'], y=df[f'TF_strength_{c2}'],
                s=50, palette=["silver","orange","purple"],ax=ax[1],ec=None,alpha=0.4,
                hue=df["sig_class"],size=df.abs_diff,sizes=(5,400))
sns.despine(trim=False)



sns.regplot(x=df[f'TF_strength_{c1}'], y=df[f'TF_strength_{c2}'],
            scatter_kws = {"color": "black","alpha": 0.7,"ec":None},
            line_kws = {"color": "red"},
            ci = 99,ax=ax[0])

for x in ax:
    x.set_xlim(-lim,lim)
    x.set_ylim(-lim,lim)
    x.set_ylabel("TF_strength_OPC")
#    x.set_ylim(-150,150)
#    x.set_ylim(-100,100)

    
#    ax[1].axline((0, 0), slope=1, color="black", ls="-", lw=0.25)
    ax[1].axline((0, threshold), slope=1, color="black", ls="--", lw=0.5)
    ax[1].axline((threshold, 0), slope=1, color="black", ls="--", lw=0.5)
    x.axhline((0) , color="gray", ls="--", lw=1,)
    x.axvline((0), color="gray", ls="--", lw=1)
    x.legend("",frameon=False);
plt.tight_layout()

plt.savefig(f"plots/CRC_TF_ntx_strength_MOL_v_OPC.png", dpi=400, bbox_inches="tight")


## Show HOX genes in the MOL v OPC TF plots

df["isHox"] = [True if x.startswith("HOX") else False for x in df.Tf]
fig,ax = plt.subplots(1,1,figsize=(7.5,7))
sns.scatterplot(data=df.sort_values(by="isHox",ascending=True),
                x=f'TF_strength_{c1}', y=f'TF_strength_{c2}',ax=ax,
                s=50, palette=["gainsboro","red"],ec=None,alpha=0.5,
                hue="isHox",)
ax.axline((0, threshold), slope=1, color="black", ls="--", lw=0.5)
ax.axline((threshold, 0), slope=1, color="black", ls="--", lw=0.5)
ax.axhline((0) , color="gray", ls="--", lw=1,)
ax.axvline((0), color="gray", ls="--", lw=1)
ax.set_ylabel("TF_strength_OPC")
sns.despine(trim=False)
plt.savefig(f"plots/CRC_TF_ntx_strength_MOL_v_OPC_HOX_only.png", dpi=400, bbox_inches="tight")




## plot TF strength

fig,ax = plt.subplots(1,6,figsize=(24,5))
#ax=ax.flatten()

for i,c1 in enumerate(["MOL","OPC","MIGL","AST","CXEX","CXINH"]):

    df = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/{c1}_roseSE/crc_out_3/{c1}_CRC_DEGREE_TABLE.txt",sep="\t")
    df[f"TF_strength_{c1}"] = df[f"In_Degree"] - df[f"Out_Degree"]
    df = df.sort_values(by=f"TF_strength_{c1}",ascending=False)
    df["TF_rank"] = [x for x in range(1,len(df)+1)]



    sns.scatterplot(data= df.sort_values(by="TF_rank"), x="TF_rank", y=f"TF_strength_{c1}",palette="RdBu_r", s=150, ec="None",ax=ax[i], alpha=0.4, hue=f"TF_strength_{c1}")
    ax[i].legend("",frameon=False);
    ax[i].axhline(0,lw=1, color="black",linestyle="--")
    if c1=="OPC_ATAC":
        ax[i].set_ylabel("TF_strength_OPC")
    elif c1=="CXINH_ATAC":
        ax[i].set_ylabel("TF_strength_CXEX")

    sns.despine()
    
plt.tight_layout()
plt.savefig(f"plots/CRC_TF_ntx_ranking_scatter_ALL.png", dpi=400, bbox_inches="tight")


## plot HOX genes in TF networks

for c1 in ["OPC","MOL"]:
    
    df = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/{c1}_roseSE/crc_out_3/{c1}_CRC_DEGREE_TABLE.txt",sep="\t")
    df[f"TF_strength_{c1}"] = df[f"In_Degree"] - df[f"Out_Degree"]
    df = df.sort_values(by=f"TF_strength_{c1}",ascending=False)
    df["TF_rank"] = [x for x in range(1,len(df)+1)]
    df["isHox"] = [True if x.startswith("HOX") else False for x in df.Tf]
    fig,ax=plt.subplots(1,1,figsize=(5,5))
    sns.scatterplot(data= df.sort_values(by="isHox"), x="TF_rank", y=f"TF_strength_{c1}",palette=["gainsboro","red"], 
                    size="isHox", sizes=(40,200),size_order=[True,False],s=150, ec="None", alpha=0.4, hue="isHox",ax=ax)
    sns.despine()
    ax.legend("",frameon=False);
    ax.axhline(0,lw=1, color="black",linestyle="--")
    ax.set_label("TF_strength_OPC")
    plt.savefig(f"plots/CRC_TF_ntx_ranking_scatter_HOX_{c1}.png",dpi=400, bbox_inches="tight")
    