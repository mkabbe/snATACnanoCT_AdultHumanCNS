import os 
import pysam
import gzip
import csv
import numpy as np
import anndata as ad
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

def removefirstPC(adata):
    # store backup
    adata.obsm["X_pca_all"] = adata.obsm["X_pca"]
    adata.uns["pca"]["variance_ratio_all"] = adata.uns["pca"]["variance_ratio"]
    adata.uns["pca"]["variance_all"] = adata.uns["pca"]["variance"]
    # replace
    adata.obsm["X_pca"] = adata.obsm["X_pca"][:,1:]
    adata.uns["pca"]["variance_ratio"] = adata.uns["pca"]["variance_ratio"][1:]
    adata.uns["pca"]["variance"] = adata.uns["pca"]["variance"][1:]


def runLSI(adata, tf_idf=True, svd=True, removeFPC=False, 
           distance="cosine", umap=True, run_norm=False,
           n_pc=50, n_n=15, umap_comp=2, louvain=True):

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
    sc.tl.umap(adata, n_components=n_comps)
    
    if louvain:
        epi.tl.louvain(adata)
    

def parse_barcodes(barcode_file):
    bcd = {}
    with open(barcode_file) as bc_csv:
        barcodes = csv.reader(bc_csv)
        for row in barcodes:
            bcd[row[0]] = row[1]
    return bcd
    
def filter_fragments_by_cluster(fragment_file, bcd):
    cluster_frag = {}
    with gzip.open(fragment_file) as f:
        tbx = pysam.tabix_iterator(f, pysam.asBed())
        for line in tbx:
            if line.name in bcd:
                try:
                    cluster_frag[bcd[line.name]].append(str(line))
                except KeyError:
                    cluster_frag[bcd[line.name]] = [str(line)]
    return cluster_frag