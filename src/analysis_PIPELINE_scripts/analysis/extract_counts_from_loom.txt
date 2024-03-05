#!/bin/python

import os
import loompy
import anndata as ad
import pandas as pd
import numpy as np
import scipy

os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE")

ks_annot = pd.read_csv("data/adult_human_RNA_siletti2022/annotation_silletti_2022.txt", sep="\t")

ls = []
for x in ks_annot["Cluster"]:
    try:
        ls.append(int(x))
    except ValueError:
        ls.append(x)
ks_annot["cluster_number"] = ls

selected_ctypes = ["OPC"] ## modify this to get multiple 
olg_clusters = [int(x) for x in ks_annot[ks_annot["Class auto-annotation"].isin(selected_ctypes)].cluster_number.tolist()]

# Connect to loom
with loompy.connect("data/adult_human_RNA_siletti2022/adult_human_20221007.loom",'r') as ds:

    ## Get column (cell) attributes for selected cells
    b00l = np.isin(ds.ca.Clusters,olg_clusters)
    print(np.sum(b00l))
    
    ca_dict = dict()
    for k,v in ds.ca.items():
        try:
            print(k)
            ca_dict[k] = v[b00l]
            
        except TypeError:
            continue
    
    ## Get count matrix
    ol = ds[:,b00l]
            
    ## Get all row (gene) attributes
    ra_dict = dict()
    for k,v in ds.ra.items():
        print(k)
        ra_dict[k] = v

## Generate anndata object
ol_X = scipy.sparse.csr_matrix(ol.T)    ## Transpose to get cells as rows (AnnData format)
ca_df = pd.DataFrame.from_dict(ca_dict)
ga_df = pd.DataFrame.from_dict(ra_dict)

adata_KS = ad.AnnData(ol_X, obs=ca_df, var=ga_df)

adata_KS.write("data/adult_human_RNA_siletti2022/adult_human_OPC_data.h5ad")

