import os
import anndata as ad
import scanpy
import episcanpy.api as epi
import argparse
import numpy as np
import pandas as pd
from metagene_functions import *


parser = argparse.ArgumentParser(description="Build a Gene Activity matrix")
parser.add_argument("-a", "--anndata", help="path to 2kb feature matrix")
parser.add_argument("-r", "--rna", help="path to RNA dataset")
parser.add_argument("-G", "--gadata", help="GA anndata file name (automatically writes to GA directory)")
parser.add_argument("-g", "--genebody", action="store_true", help="Use gene body instead of promoter")
parser.add_argument("-M", "--magic", action="store_true", help="Run MAGIC")
args = parser.parse_args()


adata = ad.read(args.anndata)

## Check if the GA or PA matrix already exists
if not args.magic:
    if "X_promoter_activity" in adata.obsm.keys() and not args.genebody:
        raise SystemExit(f"PromoterActivity already exists in {args.anndata}. Exiting.")
    elif "X_gene_activity" in adata.obsm.keys() and args.genebody:
        raise SystemExit(f"GeneActivity already exists in {args.anndata}. Exiting.")

else:
    if "X_magic_promoter_activity" in adata.obsm.keys() or "X_magic_gene_activity" in adata.obsm.keys():
        raise SystemExit(f"MAGIC imputation already exists in {args.anndata}. Exiting.")
    else:
        pass


if np.max(adata.X) != 1:
    adata.layers["normalized_counts"] = adata.X.copy()
    adata.X = adata.layers["binarized_counts"].copy()
else:
    pass

if args.genebody:
    gtf = extractGeneAnnot(gtf_file="ref/gencode.v35.annotation.gtf")
else:    
    gtf = extractPromoterCoord(gtf_file="ref/gencode.v35.annotation.gtf")

adata_features = extractFeatureCoordinates(adata)

# build gene activity matrix
gene_adata = buildGeneActivityMatrix(adata=adata,
                                     gtf=gtf, 
                                     raw_adata_features=adata_features)

if args.genebody:
    adata.obsm["X_gene_activity"] = gene_adata.X.copy()
    adata.uns["var_gene_activity"] = gene_adata.var.copy()
else:
    adata.obsm["X_promoter_activity"] = gene_adata.X.copy()
    adata.uns["var_promoter_activity"] = gene_adata.var.copy()

rna = ad.read(args.rna)
gene_modules  = extractMarkerGenes(rna_adata = rna)
addMetageneScores(adata=gene_adata, gene_modules=gene_modules)

adata.obs = gene_adata.obs.copy()                        

## run MAGIC:
if args.magic:
    rna.var_names = rna.var["gene_names"]
    var_list = rna.var_names.intersection(gene_adata.var_names)
    gene_adata_magic = imputeGA(gadata = gene_adata, gene_list=var_list)
    
    if args.genebody:
        adata.obsm["X_magic_gene_activity"] = gene_adata_magic.X.copy()
        adata.uns["var_magic_gene_activity"] = gene_adata_magic.var.copy()
    else:
        adata.obsm["X_magic_promoter_activity"] = gene_adata_magic.X.copy()
        adata.uns["var_magic_promoter_activity"] = gene_adata_magic.var.copy()

#write to files
adata.write(args.anndata)
gene_adata.write(os.path.join("data/GA_matrices",args.gadata))
