#!/bin/python
import anndata as ad
import pandas as pd
import numpy as np
import sys

## refactored the _geneactivity() function from episcanpy.api
##TODO -> add the gene activity matrix as a new layer to the existing anndata object

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


def regionsToFeatures(coords, bedout):
    '''
    write the extracted regions/gene coordinates to a .bed file
    '''
    if len(coords.keys()) == 0:
        print("fragments file not in the specified path\n")
        sys.exit(1)

    with open(bedout,"w") as _out_:
        for chrom in gtf.keys():
            for gene in gtf[chrom]:
                meta = gene[-1]
                gene_name = meta[2].lstrip(' gene_name "').rstrip('"')
                _out_.write(chrom + "\t" + str(gene[0]) + "\t" + str(gene[1]) + "\t" + gene_name + "\n")



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
    return adata_features
    

def buildGeneActivityMatrix(adata, raw_adata_features, gtf, feature_type="gene"):
    """
    find the features overlaping the genes and build the count matrix
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


def addMetageneScores(adata, gene_modules):
    '''
    assigns a metagene celltype score to each cell. score is the sum of the reads in the marker genes for each celltype.
    gene_modules should be a dict: {"OPC": ["SOX10", "PDGFRA", "PTPRZ1"], "ASTRO": ["AQP4", "GFAP"]} 
    '''
    
    celltype_marker_index = {}
    
    if "index" not in adata.var.keys():
        adata.var["index"] = [i for i in range(0,len(adata.var))]
    
    
    for module in gene_modules.keys():
        celltype_marker_index[module] = []
        for gene in gene_modules[module]:
            celltype_marker_index[module].append(adata.var.loc[adata.var["gene_name"] == gene, "index"].iloc[0])

        adata.obs[module+"_score"] = np.sum(adata.X[:,celltype_marker_index[module]], axis=1)
