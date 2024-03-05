import os
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import pyBigWig
import glob

os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE")

# read in AnnData
adata = ad.read("../../CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/220906_2kb_matrix_iterativeLSI_results_hiQ_annotated.h5ad")

##
## Might want to also read in the original 2kb matrix with NO bins filtered out
##
fulldata = ad.read("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/matrices/211215_2kb_matrix.h5ad")


## Modify the cell metadata
adata.obs["new_annot"] = adata.obs.annot_rough.tolist()
adata.obs.loc[adata.obs["annot_rough"].isin(["AST1","AST2"]), "new_annot"] = "AST"
adata.obs.loc[adata.obs["annot_rough"].isin(["OLIGO1","OLIGO2","OLIGO3"]), "new_annot"] = "MOL"
adata.obs.loc[adata.obs["annot_rough"]=="MIGL2", "new_annot"] = "ENDO"
adata.obs.loc[adata.obs["annot_rough"]=="MIGL1", "new_annot"] = "MIGL"
adata = adata[(adata.obs.annot_rough!="UNKNOWN")&(adata.obs.predicted_doublets==False)]

annotations = pd.read_csv("/proj/tmp/mukund/MK_atac_RNAannotations.csv", index_col=0)
annotations = annotations.reindex(index=adata.obs.index)
#annotations = annotations.reset_index()
annotation_dict = dict(zip(annotations.index, zip(annotations.MK_ATAC, annotations.LS_RNA, annotations.Nat19_RNA)))
adata.obs["LS_RNA"] = [annotation_dict[k][1] for k in adata.obs.index]
adata.obs["Nat19_RNA"] = [annotation_dict[k][2] for k in adata.obs.index]



## Add genomic interval info to feature metadata
adata.var["chrom"] = [x.split("_")[0] for x in adata.var.index]
adata.var["start"] = [int(x.split("_")[1]) for x in adata.var.index]
adata.var["end"] = [int(x.split("_")[2]) for x in adata.var.index]
#adata.var


## Subset the object to the celltype of interest
smalldata = adata[adata.obs.new_annot=="OPC"]


# Set the feature window for plotting tracks
window = "chr2-10_000_000-25_000_000"
c,s,e = window.split("-")

####
## Plot Single-cell tracks

# Subset object to features falling within the window
xdata = smalldata[:,(smalldata.var.chrom==c) & (smalldata.var.start>=int(s)) & (smalldata.var.end<=int(e))]

# Convert to dataFrame and sort cells by N_counts
df = pd.DataFrame.sparse.from_spmatrix(xdata.X)
df["sum"] = df.sum(axis=1)
df = df.sort_values(by="sum", ascending=False).iloc[:, :-1]


## Plot heatmap
plt.figure(figsize=(8,4))
sns.heatmap(df, cmap="binary",
            yticklabels=False, xticklabels=False,vmax=2,cbar=False,
            )# row_cluster=True)
        
## Plot heatmap as 2d array
def heatmap2d(arr: np.ndarray,
             cmap="viridis",
             vmax=2):
    plt.rcParams["figure.figsize"] = (13 ,3)
    plt.imshow(arr, cmap=cmap, vmax=vmax, aspect="auto")
    plt.colorbar()
    plt.show()

test_array = df.values#.reshape(50,50)
heatmap2d(test_array, cmap="binary", vmax=1.1)



## Plot Aggregate signal track

#BW file for celltype of interest
bw_file = "../../CZI/episcanpy_analysis/AGG_ATAC_210218/data/cluster_fragments/OLD/ALL_CELL_TYPES/BW_FILES/all_OPC_1.bed.rpkm.bw"


fig, ax = plt.subplots(1, 1, figsize=(20, 4))
#for c, ran in enumerate(ranges):
data = []
for i in range(1):
    with pyBigWig.open(bw_file) as bw:
        data.append(np.nan_to_num(bw.values(c, int(s), int(e), numpy=True),0))
#    data.append(np.nan_to_num(bw_file.values(ran[1], ran[2], ran[3], numpy=True),0))

x = np.fromiter(range(len(data[0])), dtype=int)
height = max([max(x) for x in data])
for i, v in enumerate(data):
    ax.plot(v, color="black", alpha=1)
    ax.fill_between(x,v, alpha=1, color = "black")
    ax.set_ylim(0,height)
#    ax.axis('off')
#ax.set_title(ran[0])
plt.tight_layout()




