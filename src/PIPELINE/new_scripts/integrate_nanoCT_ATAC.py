#!/bin/python 


## Annotate nanoCT data
## Integrate H3K27ac data with ATAC


import os 
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np

os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE")

adata = ad.read("mtx/raw_matrices_unprocessed/230612_P29054_H3K27ac_2kb_mtx.h5ad")
adata_ref = ad.read("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/220906_2kb_matrix_iterativeLSI_results_hiQ_annotated.h5ad")


sc.pp.filter_cells(adata, min_counts=1)
sc.pp.filter_genes(adata, min_cells=1)

## filter out lowest 10% of cells
min_features = np.quantile(adata.obs["n_counts"], 0.1)
sc.pp.filter_cells(adata, min_counts=min_features)
## select top 20,000 features
min_cells = np.sort(adata.var["n_cells"])[-100000]
adata = adata[:, adata.var["n_cells"] >= min_cells]


var_names = adata_ref.var_names.intersection(adata.var_names) ##get overlapping bins 
## subset both to contain only overlapping bins
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]


#sc.pp.filter_cells(adata, min_counts=500)
#sc.pp.filter_genes(adata, )
#sc.pp.filter_cells(adata, max_counts=10000)

metadata_dict = {
    "1": ["P27505_1001","SD041/19","CSC"], "2": ["P27505_1002","SD041/19","CSC"],
    "3": ["P27505_1003","SD042/14","CSC"], "4": ["P27505_1004","SD038/17","BA4"],
    "5": ["P29054_1001","SD030/18","BA4"], "6": ["P29054_1002","SD030/18","BA4"]}
adata.obs["agg_sample"] = [x.split("-")[-1] for x in adata.obs.index]
adata.obs["NGI_ID"] = [metadata_dict[x][0] for x in adata.obs.agg_sample]
adata.obs["caseNO"] = [metadata_dict[x][1] for x in adata.obs.agg_sample]
adata.obs["Tissue"] = [metadata_dict[x][2] for x in adata.obs.agg_sample]
adata.obs["modality"] = "H3K27ac"
#del adata.obs["agg_sample"] <- syntax may be specific for jupyter



## Only include the ATAC cells from same donors as in nanoCT object
adata_ref = adata_ref[adata_ref.obs.caseNO.isin(set(adata.obs.caseNO))]
adata_ref.obs["modality"] = "ATAC"



## store a copy of the whole ATAC dataset
whole_adata_ref = adata_ref.copy()

## Use only the ATAC data from the same donors and Tissue
ct_donors = set(adata.obs.caseNO) 
adata_ref = adata_ref[(adata_ref.obs["caseNO"].isin(ct_donors))&(adata_ref.obs["Tissue"]!="CB")]


print(f"Query dataset has shape {adata.shape}")
print(f"Reference dataset has shape {adata_ref.shape}")

sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref,metric="cosine")
sc.tl.umap(adata_ref)
#sc.pl.umap(adata_ref, color=["annot_rough"],frameon=False) ## jupyter

adata_ref.uns["pca"]["params"]["use_highly_variable"] = False
sc.tl.ingest(adata, adata_ref, obs='annot_rough')
adata.uns['annot_rough'] = adata_ref.uns['annot_rough_colors']  # fix colors


adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
adata_concat.obs.annot_rough = adata_concat.obs.annot_rough.astype('category')
adata_concat.obs.annot_rough.cat.reorder_categories(adata_ref.obs.annot_rough.cat.categories, inplace=True)  # fix category ordering
adata_concat.uns['annot_rough_colors'] = adata_ref.uns['annot_rough_colors']  # fix category colors


acdata = ad.read("mtx/230612_P29054_H3K27ac_10kb_mtx_RESULTS.h5ad")
annot_df = adata.obs[["annot_rough"]]
annot_df["barcode"] = annot_df.index.copy()
annot_bcd = set(annot_df.barcode)
acdata.obs["transferred_annot"] = [annot_df.loc[annot_df["barcode"]==x,"annot_rough"].iloc[0] if x in annot_bcd else "not_labelled" for x in acdata.obs.index]

acdata.write("mtx/230612_P29054_H3K27ac_10kb_mtx_RESULTS.h5ad")

## Add labels to H3K27me3 matrix

medata = ad.read("mtx/230612_P29054_H3K27me3_10kb_mtx_RESULTS.h5ad")

mdata_cells = medata.obs.index.copy()
adata_cells = acdata.obs.index.copy()
## Save shared cells as csv
shared_cells = set(mdata_cells.intersection(adata_cells))
pd.DataFrame(list(shared_cells), columns=["shared_barcode"]).to_csv("ref/P29054_nanoCT_shared_modality_barcodes.csv")

acdata.obs["shared_cell"] = ["True" if x in shared_cells else "False" for x in acdata.obs.index]
medata.obs["shared_cell"] = ["True" if x in shared_cells else "False" for x in medata.obs.index]

acdata = acdata[acdata.obs.shared_cell=="True"]
medata = medata[medata.obs.shared_cell=="True"]

acdata.obs["barcode"] = acdata.obs.index.copy()
medata.obs["transferred_annot"] = [acdata.obs.loc[acdata.obs.barcode==x, "transferred_annot"].iloc[0] for x in medata.obs.index]

acdata.write("mtx/230612_P29054_H3K27ac_10kb_mtx_RESULTS_SHARED_CELLS.h5ad")
medata.write("mtx/230612_P29054_H3K27me3_10kb_mtx_RESULTS_SHARED_CELLS.h5ad")

## Store OLG region barcodes
adata = ad.read("mtx/230612_P29054_H3K27ac_10kb_mtx_RESULTS_SHARED_CELLS.h5ad")
adata = adata[adata.obs.transferred_annot.isin(["OLIGO1", "OLIGO2", "OLIGO3","OPC"])]
adata.obs["OLG_tissue"] = "OLG_" + adata.obs["Tissue"].astype(str)
adata.obs[["barcode","OLG_tissue"]].to_csv("ref/OLG_region_H3K27ac_H3K27me3_P29054_barcodes.csv", header=False,index=None)

