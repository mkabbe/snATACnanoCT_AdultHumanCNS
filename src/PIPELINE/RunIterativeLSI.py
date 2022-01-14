import os
import logging as logg
import anndata as ad
import scanpy as sc
import episcanpy.api as epi
import numpy as np
import pandas as pd
import scipy
from sklearn.decomposition import TruncatedSVD
from sklearn.utils import check_random_state
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
from time import time
from datetime import datetime
date = datetime.today().strftime('%Y-%m-%d')
#%config InlineBackend.figure_format = "retina"
#plt.rcParams["figure.figsize"] = (9.6,7.2)
os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE")

logg.basicConfig(filename="logs/220111_IterativeLSI.log", format='%(asctime)s: %(message)s', level=logg.ERROR)

## simple timer
def timer(func):
    def f(*args, **kwargs):
        before = time()
        rv = func(*args, **kwargs)
        after = time()
        logg.error(f'Total time elapsed: {np.round((after - before)/60)} minutes.')
        return rv
    return f

# Run TFIDF Normalization
def runTfIdf(
    adata_atac, 
    scale_factor=1e4, 
    copy=False):
    """
    TF-IDF normalization according to MultiVelo
    """
    npeaks = adata_atac.X.sum(1)
    npeaks_inv = scipy.sparse.csr_matrix(1.0/npeaks)
    tf = adata_atac.X.multiply(npeaks_inv)
    idf = scipy.sparse.diags(np.ravel(adata_atac.X.shape[0] / adata_atac.X.sum(0))).log1p()
    if copy:
        adata_atac_copy = adata_atac.copy()
        adata_atac_copy.X = tf.dot(idf) * scale_factor
        return adata_atac_copy
    else:
        adata_atac.X = tf.dot(idf) * scale_factor

## Run SVD
def runSVD(
    adata_atac, 
    algorithm="arpack", 
    n_components=30, 
    pca_key="X_pca"):
    random_state = np.random.RandomState(0)
    random_state = check_random_state(random_state)
    pca_ = TruncatedSVD(n_components=n_components, random_state=random_state, algorithm=algorithm)
    X_pca = pca_.fit_transform(adata_atac.X)
    
    adata_atac.obsm[pca_key] = X_pca
    adata_atac.uns["pca"] = {}
    adata_atac.uns["pca"]["params"] = {"zero_center": False}
    adata_atac.uns["pca"]["variance"] = pca_.explained_variance_
    adata_atac.uns["pca"]["variance_ratio"] = pca_.explained_variance_ratio_
    adata_atac.varm["PCs"] = pca_.components_.T

## Flatten a list of lists
def flatten(t):
    return [i for s in t for i in s]

## Select top variable features across each cluster 
def selectVarFeatures(adata_atac, use_key=None):
    var_features = []
    for cluster in set(adata_atac.obs[use_key]):
        minidata = adata_atac[adata_atac.obs[use_key]==cluster].copy()
        if minidata.shape[0] < 500:
            continue
        else:
            minidata.var["n_cells"] = np.sum(minidata.layers["raw_counts"], axis=0).T
            minidata.var["z_score_acc"] = scipy.stats.zscore(minidata.var["n_cells"])
            minidata.var["abs_z"] = [abs(x) for x in minidata.var["z_score_acc"]]
            var_features.append(minidata[:,minidata.var["abs_z"]>1.65].var_names.tolist())

    var_features = list(set(flatten(var_features)))
    return adata_atac[:,var_features].copy()

@timer
def iterativeLSI(
    adata_atac,
    n_pc_comps=30,
    n_iters=3,
    cluster_res=0.5,
    distance_metric="cosine",
    clustering_method="leiden",
    save_iter_umaps=False,
    file_format = "png",
    use_binary_counts=False):
    """
    Run iterative LSI on the AnnData object. 
    TFIDF -> SVD -> KNN --> Cluster --> variable_features --> TFIDF ... 
    """
    
    #check layers in adata
    if "raw_counts" not in adata_atac.layers:
        logg.error("Raw Counts not saved. Saving now.")
        if np.max(adata_atac.X) > 1:
            adata_atac.layers["raw_counts"] = adata_atac.X.copy()
        else:
            raise Exception("Raw counts missing")
        if use_binary_counts:
            logg.error("Using binary counts. Binarizing matrix.")
            epi.pp.binarize(adata_atac)
            adata_atac.layer["binary_counts"] = adata_atac.X.copy()
        else:
            logg.error("Not using binary counts for LSI. Proceeding with raw counts.")
    

    adata_i = adata_atac.copy()

    logg.error(f"Total number of iterations: {n_iters}")
    for i in range(1, n_iters+1):
        logg.error(f"**** Running iteration {i}")
        logg.error(f"** Using {str(adata_i.shape[-1])} features")

        ## run LSI
        logg.error(" Normalizing with TFIDF")
        runTfIdf(adata_i)
        logg.error(" Running SVD")
        runSVD(adata_i, n_components=n_pc_comps)
        
        ## build kNN graph
        logg.error(" Building graph")
        sc.pp.neighbors(adata_i,  n_pcs=n_pc_comps, metric=distance_metric)
        
        logg.error(" Clustering")
        ## Clustering
        if clustering_method == "louvain":
            cluster_key = "louvain_iter_"
            cluster_key_added = cluster_key + str(i)
            try:
                sc.tl.louvain(adata_i, resolution=cluster_res, key_added=cluster_key_added)
            except:           ## Not sure what exception to catch here --> check this!
                raise Exception("Can't run -  Fix the louvain dependency issue first")
        elif clustering_method == "leiden":
            logg.error(f" Using Leiden algorithm")
            cluster_key = "leiden_iter_"
            cluster_key_added = cluster_key + str(i)
            sc.tl.leiden(adata_i, resolution=cluster_res, key_added=cluster_key_added)
        
        logg.error(f" Saving cluster labels from this round")
        ## add clustering result to adata_atac
        adata_atac.obs[cluster_key_added] = adata_i.obs[cluster_key_added].copy()


        if i != n_iters:
            ## save intermediate UMAP if needed
            if save_iter_umaps:
                logg.error(f" Saving iteration {i} umap")
                sc.tl.umap(adata_i)
                sc.pl.umap(adata_i, color=cluster_key_added, size=3, save=f"_{date}_IterativeLSI_iteration_{i}.{file_format}", show=False)
            
            ## select variable features from original adata_atac dataset
            logg.error("** Selecting variable features")
            adata_i = selectVarFeatures(adata_atac, use_key=cluster_key_added)
            adata_i.X = adata_i.layers["raw_counts"].copy()

        else: ## Final iteration
            logg.error(f" Saving UMAP for final iteration")
            sc.tl.umap(adata_i)
            sc.pl.umap(adata_i, color=cluster_key_added, size=3, save=f"_{date}_IterativeLSI_iteration_{i}.{file_format}", show=False)
            adata_atac.obsp = adata_i.obsp.copy()
            adata_atac.uns = adata_i.uns.copy()
            adata_atac.obsm = adata_i.obsm.copy()
            logg.error("****** Finished iterative clustering ******")


adata = ad.read("data/matrices/211215_2kb_matrix.h5ad")

logg.error("Pre-processing AnnData object")
## Remove cells with <3,000 and >30,000 features
epi.pp.filter_cells(adata, min_features=1000)
epi.pp.filter_cells(adata, max_features=30000)
## Select only top 100,000 accessible features
epi.pp.filter_features(adata, min_cells=1)

## Remove bad QC samples identified by Jilin and the SD036/16 samples
bad_samples = ["P20057_1002","P20602_1005","P20603_1007","P20908_1008","P22551_1002","P22551_1004","P22551_1005","P22551_1006","P22551_1008","P20908_1017","P20602_1002"]
adata = adata[~adata.obs["NGI_ID"].isin(bad_samples)].copy()

## Select only top 100,000 accessible features
epi.pp.filter_features(adata, min_cells=1)
min_cells = np.sort(adata.var["n_cells"])[-300000]
adata = adata[:, adata.var["n_cells"] >= min_cells].copy()


adata.obs["log_nb_features"] = [np.log10(x) for x in adata.obs["nb_features"]]

### Do the iterative LSI
iterativeLSI(adata, n_iters=3, save_iter_umaps=True)

adata.write("data/matrices/211215_2kb_matrix_iterativeLSI_results__testing3iters.h5ad")
