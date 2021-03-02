library(Signac)
library(Seurat)
library(dplyr)
library(reticulate)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(rtracklayer)

setwd("//Volumes/shared/Molekylär Neurobiologi/Castelo-Branco/Mukund/CZI_ADULT/CZI_ATAC/")
getwd()

remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

DATA_PATH <- "../../../NGSDATA/scATAC_human_adult_CZI_1/Processed_data/AGG_ATAC_210218/outs/"
FRAGMENTS_PATH <- paste0(DATA_PATH,"fragments.tsv.gz")
CELLS <- read.csv("data/epi/LSI-2_barcodes.txt", header=F)
FEATURES_PATH <- "data/epi/merged_louvain_cluster_peaks.bed"
FEATURES <- rtracklayer::import(FEATURES_PATH)
FRAGMENTS <- CreateFragmentObject(FRAGMENTS_PATH)

count.matrix.peaks <- FeatureMatrix(fragments = FRAGMENTS,
                                    features = FEATURES,
                                    cells = CELLS)

peak_matrix <- FeatureMatrix(fragments = FRAGMENTS,
                             features = FEATURES,
                             cells = CELLS)
packageVersion("Signac")

