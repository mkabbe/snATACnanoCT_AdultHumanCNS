library(Seurat)
library(scater)
library(dplyr)
library(reticulate)
#library(loomR)


## Read in RNA dataset
luise_rna <- readRDS("luise_Seurat/HCA_rough_annotated_all.RDS")


## load anndata
anndata <- import("anndata", convert = FALSE)

## save gene names
gene_names <- rownames(luise_rna[["RNA"]]@counts)

## build the anndata object
adata <- anndata$AnnData(
  X = t(GetAssayData(object = luise_rna)),
  obs = data.frame(luise_rna@meta.data),
  var = data.frame(gene_names),
  obsm  = list(
    "pca" = Embeddings(luise_rna[["pca"]]),
    "umap" = Embeddings(luise_rna[["umap"]]),
    "tsne" = Embeddings(luise_rna[["tsne"]])
  )
)

## save as .h5ad
anndata$AnnData$write(adata, 'HCA_rough_annotated_all.h5ad')