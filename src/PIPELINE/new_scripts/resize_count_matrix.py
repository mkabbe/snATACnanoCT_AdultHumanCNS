#!/bin/python

### Date: 221210
### Author: Mukund Kabbe

## This script generates a new count matrix
## by resizing the genomic bins 
## Currently, converts a 2kb count matrix to a 5kb count matrix

## NOTES:
# This can also be done by building a new count matrix 
# from scratch with 5kb bins as the feature set, 
# but this implementation is faster
 


import os
import anndata as ad
import numpy as np
import pandas as pd
import scipy
from numba import jit

os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE")


def extractFeatureCoordinates(adata):
    """
    get coordinates of features from AnnData object.
    """

    raw_adata = adata.copy()
    raw_adata_features = {}
    feature_index = 0
    for line in raw_adata.var_names.tolist():
        line = line.split('_')
        if line[0] not in raw_adata_features.keys():
            raw_adata_features[line[0]] = [[int(line[1]),int(line[2]), feature_index]]
        else:
            raw_adata_features[line[0]].append([int(line[1]),int(line[2]), feature_index])
        feature_index += 1
    return raw_adata_features

def extractFeatureCoordinates_2(bin_csv, valid_chroms):
    """
    New function -> uses a DF. 
    get coordinates of features from DataFrame.
    """
    df = pd.read_csv(bin_csv, sep="\t", header=None)
    df = df[df[0].isin(valid_chroms)]
    df["name"] = df[0]+"_"+df[1].astype(str)+"_"+df[2].astype(str)

    raw_adata_features = {}
    feature_index = 0
    for line in df.name.tolist():
        line = line.split('_')
        if line[0] not in raw_adata_features.keys():
            raw_adata_features[line[0]] = [[int(line[1]),int(line[2]), feature_index]]
        else:
            raw_adata_features[line[0]].append([int(line[1]),int(line[2]), feature_index])
        feature_index += 1
    return raw_adata_features

@jit(parallel=True)
def resizeCountMatrix(old_features, new_features, chrom_indices, count_array):
    ## iterate through the chrom_indices
    ## outer loop is for NEW features

    old_features = np.array(old_features)
    new_features = np.array(new_features)
    chrom_indices = np.array(chrom_indices)
    count_array = np.array(count_array)
    
    gene_index = []
    activity_X = []
    
    for chrom in chrom_indices:
        chrom_index = 0
        previous_features_index = 0
        
        for gene in new_features[chrom]:
            gene_values = []
            gene_start = gene[0]
            gene_end = gene[1]
            if gene[-1]%1000 == 0:
                print(f"{gene[-1]} features complete")
            
            for feature in old_features[chrom]:
                feature_index = 0
                if (feature[1]<= gene_start): # the window is before the feature. test the next window.
                    continue
                elif (gene_end <= feature[0]): # the window is after the feature. test the next feature.
                    break
                else: # the window overlaps the feature.
                    gene_values.append(count_array[chrom][:,feature[2]]) 
            
            if gene_values != []:
                activity_X.append(np.array(gene_values).sum(0))
                #gene_index.append(gene[-1])
            else: ## add row of zeros for the new feature 
                activity_X.append(np.zeros((1,adata.X.shape[0])).flatten())
    
    #activity_X = np.concatenate(tuple(activity_X), axis=-1)
    activity_X = np.array(activity_X)
    return activity_X ## NOTE:  shape is (n_bins x n_cells) . Need to transpose for AnnData creation!


## Read in annotated AnnData object
oldata = ad.read("../../CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/220906_2kb_matrix_iterativeLSI_results_hiQ_annotated.h5ad")
## Filter adata to keep only OLGs
ol_barcodes = set(oldata[oldata.obs.annot_rough.isin(["OPC","OLIGO1","OLIGO2","OLIGO3"])].obs.index)
## Read in AnnData containing ALL 2kb features
del oldata
adata = ad.read("../../CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/matrices/211215_2kb_matrix.h5ad")
adata = adata[adata.obs.index.isin(ol_barcodes)]
del adata.obs["nb_features"]
del adata.obs["log_nb_features"]
adata

## current bins as dict (2kb)
bins_2kb_dict = extractFeatureCoordinates(adata)
valid_chroms = bins_2kb_dict.keys()
bins_5kb_dict = extractFeatureCoordinates_2(bin_csv, valid_chroms)

## store everything as a list of list to allow 

## map chrom_name to int in list -> "chr1"=0..."chrY"=23
chrom_as_int = [x for x in range(len(valid_chroms))]
chrom_map = dict(zip(valid_chroms,chrom_as_int))

## store list of features per chrom in a new list --> [[[0,2000,0],[0,4000,1]..], [[0,2000,0],[0,4000,1]..]..]
##                                                              chr1            ,          chr2            ..
## each element in the list is a LIST of features (also as list) in that chromosome.
list_2kb_bins = []
for k in (valid_chroms):
    list_2kb_bins.append(bins_2kb_dict[k])

## repeat for the 5kb bins
list_5kb_bins = [] 
for k in (valid_chroms):
    list_5kb_bins.append(bins_5kb_dict[k])


## convert adata.X to dense numpy array - chromosome-wise!
## X_arr is a list containing .X as a numpy array 
## Each element in X_arr is the .X matrix from EACH chromosome.

adata.var["chrom"] = [x.split("_")[0] for x in adata.var.index]
adata.var["chrom_index"] = [chrom_map[k] for k in adata.var["chrom"]]

chrom_lengths = []
for chrom in chrom_as_int:
    length = len(list_2kb_bins[chrom])
    chrom_lengths.append(length)

X_arr = []
for chrom in chrom_as_int:
    start = 0
    if not chrom==0: ## change start if not the first chromosome
        start += np.sum(chrom_lengths[:chrom])
    end = start + chrom_lengths[chrom]
    print(f"Slicing -> adata[:,{start}:{end}]")
    
    ## slice AnnData and convert matrix to Array
    try:
        print(f"Converting chromosome {chrom}")
        chrom_X = adata[:,start:end].X.toarray()
        print(f"Appending array to list")        
        X_arr.append(chrom_X)
    except MemoryError:
        print("Too much memory required. Reduce the number of chromosomes.")
        print(f"Managed to complete {chrom} chromosomes")


## Actually build the new matrix
for WHICH_CHROM in range(len(chrom_as_int)):
    
    ## LOAD the old matrix
    X_arr = []
    for chrom in chrom_as_int[WHICH_CHROM:WHICH_CHROM+1]:
        start = 0
        if not chrom==0: ## change start if not the first chromosome
            start += np.sum(chrom_lengths[:chrom])
        end = start + chrom_lengths[chrom]
        print(f"Slicing -> adata[:,{start}:{end}]")

        try:
            print(f"Converting chromosome {chrom}")
            chrom_X = adata[:,start:end].X.toarray()
            print(f"Appending.")        
            X_arr.append(chrom_X)
        except MemoryError:
            print("Running out of memory.")

    ## BUILD the new matrix
    new_count_arr = resizeCountMatrix(old_features = [list_2kb_bins[WHICH_CHROM]], 
                                      new_features = [list_5kb_bins[WHICH_CHROM]], 
                                      chrom_indices = [chrom_as_int[WHICH_CHROM]],
                                      count_array = [X_arr[0]],
                                      s_pos =list_2kb_bins[WHICH_CHROM][0][-1])

    ## SAVE the new matrix
    if len(list_5kb_bins[WHICH_CHROM]) == new_count_arr.shape[0]:

        new_count_arr = scipy.sparse.csr_matrix(new_count_arr)   ## convert to AnnData
        minidata = ad.AnnData(X=new_count_arr.T, obs=adata.obs.copy())
        minidata.write(f"../../CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/matrices/ATAC_5kb_chrom_objects/chr_index_{str(WHICH_CHROM)}.h5ad")

    else:
        print("*************")
        print("*************\n")

        print("*** ERROR ***\n")

        print("Number of bins expected does not match number of bins in the matrix!")
        print(" Not converting to csr, not writing to file.\n")

        print("*************")
        print("*************")
        break




## Store feature names
names_5kb_bins = {}
for k in bins_5kb_dict.keys():
    for feature in bins_5kb_dict[k]:
        s = feature[0]
        e = feature[1]
        try:            
            names_5kb_bins[k].append(f"{k}_{s}_{e}")
        except KeyError:
            names_5kb_bins[k] = [f"{k}_{s}_{e}"]
            
## read in all AnnData objects, add var_names, store in dict for merging
ad_dict = {}
ad_dict = {}
for chrom in chrom_map.keys():
    idx = chrom_map[chrom] # this is needed to call the AnnData objects from file
    minidata = ad.read(f"../../CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/matrices/ATAC_5kb_chrom_objects/chr_index_{idx}.h5ad")
    minidata.var.index = names_5kb_bins[chrom]
    sc.pp.filter_genes(minidata, min_counts=10) # remove features with low counts
    ad_dict[idx] = minidata
    
## Merging AnnData objects 
ad_list = []
for key in ad_dict:
    ad_list.append(ad_dict[key])


adata = ad.concat(ad_list, axis=1)

adata.write("../../CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/matrices/221210_5kb_matrix_OLG_hiQ_annotated.h5ad")
