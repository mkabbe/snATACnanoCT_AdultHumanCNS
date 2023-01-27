#!/bin/python

import os
import anndata as ad
import pandas as pd
import seaborn as sns


os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE")
adata = ad.read("../../CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/220906_2kb_matrix_iterativeLSI_results_hiQ_annotated.h5ad")
annotations = pd.read_csv("/proj/tmp/mukund/MK_atac_RNAannotations.csv", index_col=0)
annotations = annotations.reindex(index=adata.obs.index)
annotation_dict = dict(zip(annotations.index, zip(annotations.MK_ATAC, annotations.LS_RNA, annotations.Nat19_RNA)))
adata.obs["LS_RNA"] = [annotation_dict[k][1] for k in adata.obs.index]
adata.obs["Nat19_RNA"] = [annotation_dict[k][2] for k in adata.obs.index]

## update a few cluster names 

adata.obs["annot_rough2"] = adata.obs["annot_rough2"].astype(str)

adata.obs.loc[adata.obs.annot_rough=="MIGL2", "annot_rough2"] = "ENDO"
adata.obs.loc[adata.obs.annot_rough=="OLIGO1", "annot_rough2"] = "MOL1"
adata.obs.loc[adata.obs.annot_rough=="OLIGO2", "annot_rough2"] = "MOL2"
adata.obs.loc[adata.obs.annot_rough=="OLIGO3", "annot_rough2"] = "MOL3"

adata.obs["annot_rough"] = adata.obs["annot_rough2"].astype("category")


## Heatmap showing annotation fidelity between metagene labels and CCA transferred labels  


df_total = (
    adata.obs.groupby("annot_rough").
        size().
        reset_index(name="n_total").
        set_index("annot_rough")
)

# Calculate number of cells for each pair of celltype annotation
df = (
    adata.obs.groupby(["annot_rough", "Nat19_RNA"]).
        size().
        reset_index(name="n").
        set_index("annot_rough").
        join(df_total).
        reset_index()
)

df_frac = df.assign(frac = lambda x: x.n / x.n_total)
df_wide = df_frac.set_index("annot_rough").pivot(columns="Nat19_RNA", values="frac")

ax = sns.heatmap(df_wide, cmap="inferno")
ax.set_title('Jäkel_Agirre_2019_annotation')
## Repeat with LS_RNA and change plot title to Seeker_2022_annotation
