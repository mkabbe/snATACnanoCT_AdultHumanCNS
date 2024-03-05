from scipy import sparse, io
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import random



adata = ad.read("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/chromvar_analysis_FBP/motif_deviations_chromVar_kNN_UMAP.h5ad")
## Modify the cell metadata
adata.obs["new_annot"] = adata.obs.annot_rough.tolist()
adata.obs.loc[adata.obs["annot_rough"].isin(["AST1","AST2"]), "new_annot"] = "AST"
adata.obs.loc[adata.obs["annot_rough"].isin(["CBEX","CBGRA"]), "new_annot"] = "CBEX"
adata.obs.loc[adata.obs["annot_rough"].isin(["OLIGO1","OLIGO2","OLIGO3"]), "new_annot"] = "MOL"
adata.obs.loc[adata.obs["annot_rough"]=="MIGL2", "new_annot"] = "ENDO"
adata.obs.loc[adata.obs["annot_rough"]=="MIGL1", "new_annot"] = "MIGL"
adata = adata[(adata.obs.annot_rough!="UNKNOWN")]



sparsedeviations = io.mmread("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/chromvar_analysis_FBP/deviation.mtx")

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

var = "new_annot" 
var_data = adata.obs.loc[:,var].sort_values()
cell_topic = chromvar_df.loc[var_data.index.to_list()]

df = pd.concat([cell_topic, var_data], axis=1, sort=False)

topic_order = df.groupby(var).mean().idxmax().sort_values().index.to_list()
cell_topic = cell_topic.loc[:, topic_order]
#df = df.drop("annot_rough",axis=1)

df_all = df.drop("new_annot", axis=1)

## This block sets the family name for the different motifs
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


# Now re-arrange df_all.columns in the order of all_tf_table.jaspar_name
# Then create color map for the TF families
# Plot clustermap

df_all_common_tf = df_all[all_tf_table.jaspar_name]

lut = dict(zip(df.new_annot.unique(), adata.uns["new_annot_colors"]))
col_colors = var_data.map(lut)

import random
tf_colors = list(map(lambda i: "#" + "%06x" % random.randint(0, 0xFFFFFF), range(len(set(all_tf_table.Family_Name)))))
tf_lut = dict(zip(all_tf_table.Family_Name.unique(), tf_colors))
row_colors = all_tf_table.Family_Name.map(tf_lut)

all_tf_table["family_color"] = [tf_lut[k] for k in all_tf_table.Family_Name]


vmax = np.percentile(df_all, 99.9999)
vmin = np.percentile(df_all, 0.0001)
g = sns.clustermap(df_all.T, col_colors = col_colors.values, vmax=1, vmin=-1,
               figsize=(8, 10),yticklabels=False, xticklabels=False,
               row_cluster=True, col_cluster=False, cmap="coolwarm");
g.fig.subplots_adjust(right=0.7)
g.ax_cbar.set_position((0.75, .3, .03, .2))

g.savefig("plots/ATAC_chromvar_deviation_matrix.png", dpi=400, bbox_inches="tight")
