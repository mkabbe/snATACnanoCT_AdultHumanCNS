import os
import pandas as pd
import anndata as ad

magic = ad.read("data/magic/210527_gene_metagene_activity_magic_matrix.h5ad")

### Annotating clusters ###
cluster_annot = {"MIGL_1":["3"], "OPC_1":["10"], 
"Oligo_1":["0","1","7","9"],
"Astro_1":["4","17"],
"CB_Ex_Neu_1":["2","8","13","18"],
"BA4_Ex_Neu_1":["5","11","12"],
"Inh_Neu_1":["14"],
"Unknown":["6","19"],
"doublet":["15","16"]}

magic.obs["annot_rough"] = ["NA" for x in magic.obs_names]

for annot in cluster_annot.keys():
    for cluster in cluster_annot[annot]:
        magic.obs.loc[magic.obs["louvain"] == cluster, "annot_rough"] = annot

### Binning ages ###
ol_adata.obs['age_bins'] = pd.cut(x=ol_adata.obs['Age'], 
                right=False, bins=[30, 40, 50, 60, 70, 80], 
                labels=["<40","<50","<60","<70","<80"])


ol_adata.obs['age_bins'] = pd.cut(x=ol_adata.obs['Age'],
                right=False, bins=[25, 35, 45, 55, 65, 75],
                labels=["25-34","35-44","45-54","55-64","65-74"])
