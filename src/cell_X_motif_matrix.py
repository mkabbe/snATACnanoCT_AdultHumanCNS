import os
import anndata as ad
import episcanpy.api as epi
import numpy as np
import scipy

os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/data/motif")

## PEAK X MOTIF
pm_adata = ad.read("../feature_matrices/210604_peak_by_motif_matrix.h5ad")

if np.max(pm_adata.X > 1):
    pm_adata.layers["raw_motif_counts"] = pm_adata.X.copy()
    epi.pp.binarize(pm_adata)
    pm_adata.layers["binary_motif_counts"] = pm_adata.X.copy()


## CELL X PEAK 
cp_adata = ad.read("../feature_matrices/210510_LSI-2_RESULTS_merged_peak_matrix.h5ad")

cp_adata.layers["normalized"] = cp_adata.X.copy() ## store the normalized counts in a layer
cp_adata.X = cp_adata.layers["raw_peak_counts"] ## load the raw counts layer

## CELL X MOTIF 
cm_mtx = cp_adata.layers["raw_peak_counts"] * pm_adata.layers["binary_motif_counts"]    ## clear and explicit
#cm_mtx = cp_adata.X * pm_adata.X

cm = ad.AnnData(cm_mtx, 
                obs=cp_adata.obs.copy(), 
                var=pm_adata.var.copy(),
                obsm = cp_adata.obsm.copy(),
                uns=cp_adata.uns.copy())

##cm.write("") 