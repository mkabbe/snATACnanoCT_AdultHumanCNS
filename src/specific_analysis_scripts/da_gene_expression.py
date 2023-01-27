#!/bin/python

import os 
import anndata as ad
import scanpy as sc
import pandas as pd
%config InlineBackend.figure_format = "retina"
os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218")

rna_ol = ad.read("data/feature_matrices/HCA_oligo_and_opc_LS.h5ad")
rna_ol = rna_ol[rna_ol.obs.Tissue.isin(["CSC","BA4"])]
rna_ol.obs["label_tissue"] = rna_ol.obs["ol_clusters_named"].astype(str) +"_"+ rna_ol.obs["Tissue"].astype(str)
rna_ol.var_names = rna_ol.var["gene_names"].tolist()

genes = pd.read_csv("data/variable_peaks/annot_sig_diff_peaks_OPC_CSCvsBA4.bed", sep="\t")["Gene Name"].to_list()
da_genes = [x for x in genes if x in rna_ol.var_names]

sc.pl.stacked_violin(rna_ol, var_names=da_genes, groupby="label_tissue")