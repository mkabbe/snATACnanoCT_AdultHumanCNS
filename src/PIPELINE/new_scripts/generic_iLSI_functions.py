#!/bin/python


import os
import csv
import scipy
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
from sklearn.feature_extraction.text import TfidfTransformer
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import style
import warnings
import pickle
import harmonypy
style.use("default")
warnings.filterwarnings('ignore')
%config InlineBackend.figure_format = "retina"
from lsi_functions import *
os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218")
plt.rcParams["figure.figsize"] = (9.6,7.2)


########################################################################
########################################################################

## Flatten a list of lists
def flatten(t):
    return [i for s in t for i in s]
    
## plot correlation between principal components and sequencing depth
def plotDepthCorr(adata, key = "X_pca", n_comps = 10):
    if key not in adata.obsm.keys():
        raise KeyError(f"{key} not found in .obsm.keys()")
            
    PC_depth_corr = []
    for PC in range(adata.obsm[key].shape[1]):
        corr = scipy.stats.pearsonr(adata.obs["nb_features"], adata.obsm[key][:,PC])[0]
        PC_depth_corr.append(corr)
    plt.plot([x for x in range(1,n_comps+1)], PC_depth_corr[:n_comps], 'o', color='black');
    plt.xticks(list(range(0,n_comps+2,2)));
    plt.yticks([1, 0.5, 0, -0.5,  -1]);
    plt.xlabel("PC");
    plt.ylabel("Correlation Coefficient");
    plt.title("Correlation between Components and Sequencing depth");

## Run TF-IDF normalization
def runTFIDF(adata, run_norm=False):
    X = adata.layers["raw_counts"]
    if run_norm:
        tfidf = TfidfTransformer(use_idf=True)
    else:
        tfidf = TfidfTransformer(use_idf=True, norm=False)
    X_norm = tfidf.fit_transform(X)
    adata.X = X_norm


## TFIDF according to MultiVelo
def tfidf_norm(adata_atac, scale_factor=1e4, copy=False):
    npeaks = adata_atac.X.sum(1)
    npeaks_inv = scipy.sparse.csr_matrix(1.0/npeaks)
    tf = adata_atac.X.multiply(npeaks_inv)
    idf = scipy.sparse.diags(np.ravel(adata_atac.X.shape[0] / adata_atac.X.sum(0)))s.log1p()
    if copy:
        adata_atac_copy = adata_atac.copy()
        adata_atac_copy.X = tf.dot(idf) * scale_factor
        return adata_atac_copy
    else:
        adata_atac.X = tf.dot(idf) * scale_factor

## Run SVD, build kNN graph and do leiden clustering
def runLSI(adata, cluster_key_added = "leiden"):
    sc.tl.pca(adata, zero_center=False, n_comps=31)
    sc.pp.neighbors(adata, n_pcs=30, metric="cosine")
    sc.tl.umap(adata)
    cluster_res = [0.5, 0.7, 0.9]
    for res in cluster_res:
        sc.tl.leiden(adata, resolution = res, key_added = f"{cluster_key_added}_{res}")

## Select top variable features across each cluster 
def selectVarFeatures(adata, use_key="leiden_0.5", threshold=1.65):
    var_features = []
    for cluster in set(adata.obs[use_key]):
        minidata = adata[adata.obs[use_key]==cluster].copy()
        if minidata.shape[0] < 500:
            continue
        else:
            minidata.var["n_cells"] = np.sum(minidata.layers["raw_counts"], axis=0).T
            minidata.var["z_score_acc"] = scipy.stats.zscore(minidata.var["n_cells"])
            minidata.var["abs_z"] = [abs(x) for x in minidata.var["z_score_acc"]]
            var_features.append(minidata[:,minidata.var["abs_z"]>threshold].var_names.tolist())

    var_features = list(set(flatten(var_features)))
    return adata[:,var_features].copy()
    

## plot stacked percentage bar plots
def plotClusterComp(adata, columns =["10X_BATCH","NGS_BATCH","caseNO","Tissue"], groups="leiden"):
    for x in columns:
        if x=="Tissue":
            cmap="coolwarm"
        else:
            cmap="tab20"
        grp_df = adata.obs.groupby([groups,x])["nb_features"].count().unstack().T
        plt_df = grp_df.div(grp_df.sum())
        p1 = plt_df.T.plot(kind='bar', stacked=True, rot=1, 
        figsize=(8,8),ylabel="fraction", edgecolor="black", 
        width = 1, colormap=cmap, alpha=0.5, title=x);
        p1.legend(loc='center left', bbox_to_anchor=(1, 0.5));


#def iterativeClustering(adata, n_iters=3, distance_metric="cosine", var_score_cutoff=1.65):
#    if "X_umap_original" not in adata.obsm.keys():
#        adata.obsm["X_umap_original"] = adata.obsm["X_umap"].copy()
#    else:
#        pass
    
    ## check if initial iteration has been done
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        