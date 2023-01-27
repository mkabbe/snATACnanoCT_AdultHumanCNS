import os
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import episcanpy.api as epi
from sklearn.feature_extraction.text import TfidfTransformer


def runTFIDF(adata, run_norm=False):
    X = adata.X
    if run_norm:
        tfidf = TfidfTransformer(use_idf=True)
    else:
        tfidf = TfidfTransformer(use_idf=True, norm=False)
    X_norm = tfidf.fit_transform(X)
    adata.X = X_norm


def runLSI(adata, tf_idf=True, svd=True, removeFPC=False, 
           distance="cosine", umap=True, run_norm=False,
           n_pc=50, n_n=15, umap_comp=2, louvain=False):
           """
           Runs LSI on the anndata object.
           Binarizes if necessary --> TF-IDF --> SVD
           Calls neighbors --> UMAP --> Clustering
           """

    ## check binary
    if np.max(adata.X) > 1:
        adata.layers["raw_counts"] = adata.X.copy()
        epi.pp.binarize(adata)
        adata.layers["binary_counts"] = adata.X.copy()
    
    ## TF-IDF normalization    
    if tf_idf:
        runTFIDF(adata, run_norm)
    
    ## run SVD
    if svd:
        sc.tl.pca(adata, zero_center=False, n_comps=51)
    
    #remove first PC ((NOT RUN BY DEFAULT))
    if removeFPC:
        removefirstPC(adata)
    
    sc.pp.neighbors(adata, n_pcs=n_pc, n_neighbors=n_n, metric=distance)
    sc.tl.umap(adata, n_components=umap_comp)
    
    if louvain:
        epi.tl.louvain(adata)
    else:
        epi.tl.leiden(adata)
