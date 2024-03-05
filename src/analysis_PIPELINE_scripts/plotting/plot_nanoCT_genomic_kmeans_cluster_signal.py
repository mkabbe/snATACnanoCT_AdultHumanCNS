import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


## plot signal distribution in identified k-means clusters

for ctype in ["MOL","OPC","MIGL","AST","CXEX","CXINH"]:
    df = pd.read_csv(f"fragments/P29054/ctype/hPTM_HCA_TSS_only_{ctype}_profile_data.tsv",sep="\t")
    df = df.T
    df.index = [x for x in range(-1,2001)]
    df = df.drop(index=[0])
    df = df.drop(columns=[0])
    df_ac = df[[x for x in range(1,11)]].drop(index=[-1])
    df_me = df[[x for x in range(11,21)]].drop(index=[-1])
    df_me.columns = [x for x in range(1,11)]

    nrows=2
    ncols=5
    fig, ax = plt.subplots(nrows=nrows,ncols=ncols, figsize=[10,3]) 
    ax = ax.flatten()
    for c in range(ncols*nrows):
        sns.lineplot(data = df_ac, x = df_ac.index, y = df_ac[c+1], linewidth=1, color="green", alpha=1, ax=ax[c])
        sns.lineplot(data = df_me, x = df_me.index, y = df_me[c+1], linewidth=1, color="red", alpha=1, ax=ax[c])
        ax[c].set_xticks([0,500,1000])
        ax[c].set_xticklabels(["-5kb","TSS","+5kb"],fontsize=5)
#        if c==0:
#            ax[c].set_ylim([0, 30])
#        else:
#            ax[c].set_ylim([0, 17])
        if c==0 or c==ncols:
            ax[c].set_ylabel('Normalized signal', fontsize=8)
        else:
            ax[c].set_ylabel('')
        ax[c].set_title(f"cluster {c+1}",fontsize=8)
    sns.despine()
    fig.suptitle(f"{ctype} signal clusters", weight="bold")
    plt.tight_layout()
fig.savefig("plots/TSS_hPTM_cluster_celltype_ALL.png", dpi=400, bbox_inches="tight")