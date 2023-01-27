#!/bin/python


### Script for plotting the chromVar deviations


import os
import csv
import scipy
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import style
from matplotlib.patches import Patch
import pickle
import random

#import warnings
#style.use("default")
#warnings.filterwarnings('ignore')
#%config InlineBackend.figure_format = "retina"
os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE")



adata = ad.read("../../CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/220906_2kb_matrix_iterativeLSI_results_hiQ_annotated.h5ad")
adata.obs["predicted_doublets"] = adata.obs["predicted_doublets"].astype(str)
annotations = pd.read_csv("/proj/tmp/mukund/MK_atac_RNAannotations.csv", index_col=0)
annotations = annotations.reindex(index=adata.obs.index)
annotation_dict = dict(zip(annotations.index, zip(annotations.MK_ATAC, annotations.LS_RNA, annotations.Nat19_RNA)))
adata.obs["LS_RNA"] = [annotation_dict[k][1] for k in adata.obs.index]
adata.obs["Nat19_RNA"] = [annotation_dict[k][2] for k in adata.obs.index]

deviations_sparse = scipy.io.mmread("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/chromvar_analysis_FBP/deviation.mtx")
row_barcode_names = pd.read_csv("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/chromvar_analysis_FBP/barcode.txt",header=None)
col_motif_names = pd.read_csv("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/chromvar_analysis_FBP/motif_name.txt", header=None)

with open("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/ref/_sampleID.pickle","rb") as f:
    ngi_to_gem_dict = pickle.load(f)
gem_to_ngi_dict = {v:k for k,v in ngi_to_gem_dict.items()}

row_barcode_names["ngi_id"] = [x.split("#")[0] for x in row_barcode_names[0]]
row_barcode_names["sample"] = [gem_to_ngi_dict[k] for k in row_barcode_names["ngi_id"]]
row_barcode_names["barcode"] = [x.split("#")[1].split("-")[0] for x in row_barcode_names[0]]
row_barcode_names["barcode"] = row_barcode_names["barcode"] + "-" + row_barcode_names["sample"]

deviations_dense = sparsedeviations.toarray()

chromvar_df = pd.DataFrame(deviations_dense.T, columns=col_motif_names[0], index=row_barcode_names["barcode"])

var = "annot_rough" 
var_data = adata.obs.loc[:,var].sort_values()
cell_topic = chromvar_df.loc[var_data.index.to_list()]

df = pd.concat([cell_topic, var_data], axis=1, sort=False)

topic_order = df.groupby(var).mean().idxmax().sort_values().index.to_list()
cell_topic = cell_topic.loc[:, topic_order]
#df = df.drop("annot_rough",axis=1)

df_all = df.drop("annot_rough", axis=1)


hoxa_set = [f"HOXA{x}" for x in range(1,14)]
hoxb_set = [f"HOXB{x}" for x in range(1,14)]
hoxc_set = [f"HOXC{x}" for x in range(1,14)]
hoxd_set = [f"HOXD{x}" for x in range(1,14)]

hox_list = []
for x in [hoxa_set,hoxb_set,hoxc_set,hoxd_set]:
    for y in x:
        hox_list.append(y)
hox_set = set(hox_list)

hox_order = []
for x in hox_list:
    if x in cell_topic.columns:
        hox_order.append(x)
        

## This block sets the family name for the different motifs

df_tf_names = [x.split("_")[0] for x in df_all.columns]

tf_csv = pd.read_csv("../../CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/chromvar_analysis_FBP/TF_name_family.txt",sep="\t")
tf_merged_set = set(set(tf_csv.TF_Name).intersection(set(df_tf_names)))
tf_csv = tf_csv.loc[tf_csv.TF_Name.isin(tf_merged_set)].drop_duplicates()
motif_family_dict = dict(zip(tf_csv.TF_Name, tf_csv.Family_Name))

all_tf_table = pd.DataFrame({"jaspar_name":df_all.columns, "motif_name":df_tf_names})
all_tf_table = all_tf_table.loc[all_tf_table.motif_name.isin(tf_merged_set)]
all_tf_table["Family_Name"] = [motif_family_dict[k] for k in all_tf_table.motif_name]

all_tf_table.loc[all_tf_table.motif_name.isin(hox_set),"Family_Name"] = "HOX_family"
all_tf_table = all_tf_table.sort_values(by=["Family_Name","motif_name"])


###

df_all_common_tf = df_all[all_tf_table.jaspar_name]
lut = dict(zip(df.annot_rough.unique(), adata.uns["annot_rough_colors"]))
col_colors = var_data.map(lut)


tf_colors = list(map(lambda i: "#" + "%06x" % random.randint(0, 0xFFFFFF), range(len(set(all_tf_table.Family_Name)))))
tf_lut = dict(zip(all_tf_table.Family_Name.unique(), tf_colors))
row_colors = all_tf_table.Family_Name.map(tf_lut)


all_tf_table["family_color"] = [tf_lut[k] for k in all_tf_table.Family_Name]

new_lut = dict(zip(all_tf_table.jaspar_name, all_tf_table.family_color))
new_lut_2 = {k:new_lut[k] for k in df_all_common_tf.columns}
new_lut_3 = dict(zip(all_tf_table.Family_Name, all_tf_table.family_color))



vmax = np.percentile(df_all_common_tf, 99.9999)
vmin = np.percentile(df_all_common_tf, 0.0001)
g = sns.clustermap(df_all_common_tf.T, col_colors = col_colors.values, row_colors=list(new_lut_2.values()), vmax=vmax, vmin=vmin,
               figsize=(12, 16),yticklabels=False, xticklabels=False,
               row_cluster=True, col_cluster=False, cmap="viridis");
g.cax.set_position([1, 0.55, 0.05, 0.2])
handles = [Patch(facecolor=new_lut_3[name]) for name in new_lut_3]
plt.legend(handles, new_lut_3, title='TF Family',
           bbox_to_anchor=(0.1, 0.8), bbox_transform=plt.gcf().transFigure);
           

g.savefig("sns-heatmap.svg")


oligo_df = df[df.annot_rough.isin(["OPC","OLIGO1","OLIGO2","OLIGO3"])]
oligo_df = oligo_df.drop("annot_rough",axis=1)
oligo_hox_df = oligo_df[[x for x in df.columns if x.startswith("OLIG") or x.startswith("HOX")  or x.startswith("SOX")]] ### Also adding OLIG and SOX



motifdata = ad.AnnData(df_all)
sparse_X = scipy.sparse.csr_matrix(motifdata.X)
motifdata.X = sparse_X.copy()

sc.pp.neighbors(motifdata, use_rep="X") ## Not computing PCA. Directly using the Z score deviations as the "latent" matrix
sc.tl.umap(motifdata)
motifdata.obs["annot_rough"] = [annotation_dict[k][0] for k in motifdata.obs.index]
motifdata.obs["LS_RNA"] = [annotation_dict[k][1] for k in motifdata.obs.index]
motifdata.obs["Nat19_RNA"] = [annotation_dict[k][2] for k in motifdata.obs.index]
motifdata.write("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/chromvar_analysis_FBP/motif_deviations_chromVar.h5ad")




