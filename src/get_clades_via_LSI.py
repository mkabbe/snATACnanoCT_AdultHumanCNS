import os
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import scipy
import episcanpy.api as epi

featMatrix = "data/epi/5kb.feature_matrix.h5ad" ## get path to file

# read in matrix and binarize
adata = ad.read(featMatrix)
epi.pp.binarize(adata)

## custom TFIDF

def TFIDF(x):
    tf_x = np.divide(x, x.sum(0))
    idf_x = np.diag(np.log(1 + (np.shape(x)[1]) / x.sum(1)))
    tfidf_x = np.dot(tf_x.T, idf_x).T
    return tfidf_x

epi.pp.filter_cells(adata, min_features=100)
epi.pp.filter_features(adata, min_cells=10)

X = adata.X.T #transpose anndata matrix to get cells as columns
tfidf_X = TFIDF(X.toarray()) # convert to nparray and run TF-IDF

adata.X = scipy.sparse.csr_matrix(tfidf_X).T #convert back to scipy.sparse and transpose to get cells as rows

## TODO: test TF-IDF from sklearn (sklearn.feature_extraction.text.TfidfTransformer)
## TODO: run SVD (sklearn.decomposition.TruncatedSVD)
## TODO: figure out how to get a heatmap showing broad clusters 

