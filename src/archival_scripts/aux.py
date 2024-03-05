#!/bin/python


import os 
import pysam
import gzip
import csv
from sklearn.feature_extraction.text import TfidfTransformer
import anndata as ad
import scanpy as sc 
import episcanpy.api as epi


def runTFIDF(adata):
    X = adata.X
    tfidf = TfidfTransformer(use_idf=True)
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
        tbx = pysam.tabix_iterator(f,pysam.asBed())
        for line in tbx:
            if line.name in bcd:
                try:
                    cluster_frag[bcd[line.name]].append(str(line))
                except KeyError:
                    cluster_frag[bcd[line.name]] = [str(line)]
    return cluster_frag





