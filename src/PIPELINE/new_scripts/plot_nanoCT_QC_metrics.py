#!/bin/python

import os
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

## Generates QC plots for the 2 nanoCT modalities
## violinplots showing log10(fragments) and FRiP values 



proj = "CZI"

if proj == "CZI":
    acdata = ad.read("mtx/230123_P27505_H3K27ac_50kb_mtx.h5ad")
    medata = ad.read("mtx/230123_P27505_H3K27me3_50kb_mtx.h5ad")
elif proj == "EAE":
    acdata = ad.read("mtx/230127_P28305_H3K27ac_50kb_mtx.h5ad")
    medata = ad.read("mtx/230127_P28305_H3K27me3_50kb_mtx.h5ad")
else:
    print("Invalid project specified.")
    

for i,adata in enumerate([acdata, medata]):
    sc.pp.filter_cells(adata, min_counts=500)
    adata.obs["log_nb_features"] = [np.log10(x) for x in adata.obs["n_counts"]]
    adata.obs["GEM_ID"] = [x.split("-")[-1] for x in adata.obs.index]
    metadict = pd.read_csv(f"ref/metadata_{proj}_nanoCT.csv").to_dict(orient="index")
    if proj == "EAE":
        for k in ["NGI_ID","NGI_BATCH","MK_ID","Condition","Expt_Batch","Replicate"]:
            adata.obs[k] = [metadict[int(v)-1][k] for v in adata.obs["GEM_ID"]]
    elif proj == "CZI":
        for k in ["NGI_ID","NGI_BATCH","MK_ID","Tissue","10X_BATCH","caseNO","Sex"]:
            adata.obs[k] = [metadict[int(v)-1][k] for v in adata.obs["GEM_ID"]]
    if i == 0:
        ac_obs = adata.obs.copy()
    else:
        me3_obs = adata.obs.copy()


bcdct = [metadict[k]["NGI_ID"] for k in metadict.keys()]

for modality in ["H3K27ac","H3K27me3"]:

    for i,sample in enumerate(bcdct):
        tmp = pd.read_csv(f"singlecell_CT_filtered/{sample}_{modality}_singlecell_CT_filtered_relaxed.csv",)[["barcode","peak_ratio_CT"]]
        tmp["bcd"] = [x.split("-")[0]+f"-{i+1}" for x in tmp["barcode"]]
        if i==0:
            tmp1 = tmp.copy()
        else:
            tmp1 = tmp1.append(tmp)
    frip_dict = dict(zip(tmp1.bcd,tmp1.peak_ratio_CT))
    if modality == "H3K27ac":
#        ac_obs = adata.obs.copy()
        ac_obs["FRiP"] = [frip_dict[x] if x in frip_dict.keys() else 0 for x in ac_obs.index]
        ac_obs["modality"] = modality
    else:
#        me3_obs = adata.obs.copy()
        me3_obs["FRiP"] = [frip_dict[x] if x in frip_dict.keys() else 0 for x in me3_obs.index]
        me3_obs["modality"] = modality
        
df = ac_obs[["modality","log_nb_features","FRiP"]].append(me3_obs[["modality","log_nb_features","FRiP"]])
df.index = [x for x in range(len(df))]


fig, axs = plt.subplots(1,2)

sns.violinplot(data=df, y="log_nb_features",x="modality", bw=0.2, width=0.4, 
               palette=["green","turquoise"], gridsize=50,inner=None,ax=axs[0])

sns.violinplot(data=df, y="FRiP",x="modality", bw=0.2, width=0.4, 
               palette=["green","turquoise"], gridsize=50,inner=None, ax=axs[1])

axs[0].set(xlabel="modality", ylabel="log10(unique fragments)")
axs[1].set(xlabel="modality", ylabel="FRiP")
sns.despine()
plt.tight_layout()

plt.savefig(f"plots/nanoCT_{proj}_modality_QC_violinplots.pdf")
plt.savefig(f"plots/nanoCT_{proj}_modality_QC_violinplots.svg")