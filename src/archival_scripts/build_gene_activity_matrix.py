#!/bin/python
import os
import csv
import scipy
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi


os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218")


def extractGeneAnnot(gtf_file, upstream = 2000, feature_type = "gene", annotation = "HAVANA"):
    """
    get coordinates and annotations of genes from a GTF file. Returns dictionary indexed by chromosome
    """
    gtf = {}
    
    with open(gtf_file) as f:
        for line in f:
            if line[0:2] != "##" and "\t"+feature_type+"\t" in line and "\t"+annotation+"\t" in line:
                line = line.rstrip("\n").split("\t")
                if line[6] == "+":
                    if line[0] not in gtf.keys():
                        gtf[line[0]] = [[int(line[3])-upstream, int(line[4]), line[-1].split(";")[:-1]]]
                    else:
                        gtf[line[0]].append([int(line[3])-upstream, int(line[4]), line[-1].split(";")[:-1]])
                else:
                    if line[0] not in gtf.keys():
                        gtf[line[0]] = [[int(line[3]), int(line[4])+upstream, line[-1].split(";")[:-1]]]
                    else:
                        gtf[line[0]].append([int(line[3]), int(line[4])+upstream, line[-1].split(";")[:-1]])
    return gtf

def extractPromoterCoord(gtf_file, window = 2000, feature_type = "gene", annotation = "HAVANA"):
    gtf = {}
    extnd = window//2
    with open(gtf_file) as f:
            for line in f:
                if line[0:2] != "##" and "\t"+feature_type+"\t" in line and "\t"+annotation+"\t" in line:
                    line = line.rstrip("\n").split("\t")
                    if line[6] == "+":
                        if line[0] not in gtf.keys():
                            gtf[line[0]] = [[int(line[3])-extnd, int(line[3])+extnd, line[-1].split(";")[:-1]]]
                        else:
                            gtf[line[0]].append([int(line[3])-extnd, int(line[3])+extnd, line[-1].split(";")[:-1]])
                    else:
                        if line[0] not in gtf.keys():
                            gtf[line[0]] = [[int(line[4])-extnd, int(line[4])+extnd, line[-1].split(";")[:-1]]]
                        else:
                            gtf[line[0]].append([int(line[4])-extnd, int(line[4])+extnd, line[-1].split(";")[:-1]])
    return gtf


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

def buildGeneActivityMatrix(adata, raw_adata_features, gtf, feature_type):
    """
    find the features overlaping the genes. and build the count matrix
    """
    gene_index = []
    gene_activity_X = []
    raw_adata = adata.copy()
    for chrom in gtf.keys():
        if chrom in raw_adata_features.keys():
            #print(chrom)
            chrom_index = 0
            previous_features_index = 0
            for gene in gtf[chrom]:
                gene_values = []
                gene_start = gene[0]
                gene_end = gene[1]
                for feature in raw_adata_features[chrom]:
                    feature_index = 0
                    if (feature[1]<= gene_start): # the window is before the gene. we need to test the next window.
                        continue
                    elif (gene_end <= feature[0]): # the window is after the gene. we need totest the next gene.
                        break
                    else: # the window is overlapping the gene. 
                        gene_values.append(raw_adata.X[:,feature[2]].todense())
                if gene_values != []:
                    gene_activity_X.append(np.sum(gene_values, axis=0))
                    gene_index.append(gene[-1])
                

    gene_activity_X = np.concatenate(tuple(gene_activity_X), axis=-1)

    # get the variable metadata
    if feature_type=='transcript':
        gene_name = [x[7].lstrip(' transcript_name "').rstrip('"') for x in gene_index]
    else:
        gene_name = [x[2].lstrip(' gene_name "').rstrip('"') for x in gene_index]

    metadata_genes = {'gene_id' : [],
                      'transcript_id' : [],
                      'gene_type' : [],
                      'gene_name' : [],
                      'transcript_type' : [],
                      'transcript_name' : [],
                      'protein_id' : []}

    for line in gene_index:
        dico_line = {}
        for element in line:
            if ' "' in element:
                dico_line[element.rstrip('"').lstrip(" ").split(' "')[0]] = element.rstrip('"').lstrip(" ").split(' "')[1]
    
        for key in metadata_genes.keys():
            if key in dico_line.keys():
                metadata_genes[key].append(dico_line[key])
            else:
                metadata_genes[key].append('NA')  
            
    dataframe_genes = pd.DataFrame.from_dict(metadata_genes)
    dataframe_genes.index = gene_name

    #adata.layers[layer_name] = ad.AnnData(gene_activity_X, var=dataframe_genes, obs=raw_adata.obs)
    gene_adata = ad.AnnData(gene_activity_X, var=dataframe_genes, obs=raw_adata.obs)
    gene_adata.uns = adata.uns.copy()
    gene_adata.obsm = adata.obsm.copy()
    gene_adata.obsp = adata.obsp.copy()
    return(gene_adata)


adata = ad.read("data/210421_merged_peak_matrix_RESULTS.h5ad")
#gtf = extractGeneAnnot(gtf_file="ref/gencode.v35.annotation.gtf")
gtf = extractPromoterCoord(gtf_file="ref/gencode.v35.annotation.gtf")
adata_features = extractFeatureCoordinates(adata)

gene_adata = buildGeneActivityMatrix(adata=adata, gtf=gtf, raw_adata_features=adata_features, feature_type="gene")

#gene_adata.write("210407_LSI-2_gene_activity_matrix.h5ad")
gene_adata.write("data/210421_promoter_activity_matrix.h5ad")

