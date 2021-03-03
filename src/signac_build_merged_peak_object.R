devtools::install_github("timoast/signac", ref = "develop")
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


DATA_PATH <- "../../../NGSDATA/scATAC_human_adult_CZI_1/Processed_data/AGG_ATAC_210218/outs/"
FRAGMENTS_PATH <- paste0(DATA_PATH,"fragments.tsv.gz")
CELLS <- read.csv("data/epi/LSI-2_barcodes.txt", header=F)
FEATURES_PATH <- "data/epi/merged_louvain_cluster_peaks.bed"
FEATURES <- rtracklayer::import(FEATURES_PATH)
FRAGMENTS <- CreateFragmentObject(FRAGMENTS_PATH)

cell_bc = c()
for (i in CELLS) {
  cell_bc <- append(cell_bc, i)
}

count.matrix.peaks <- FeatureMatrix(fragments = FRAGMENTS,
                                    features = FEATURES,
                                    cells = cell_bc)

chrom_assay <- CreateChromatinAssay(
  counts = count.matrix.peaks,
  sep = c(":","-"),
  genome = "hg38",
  fragments = FRAGMENTS_PATH,
  min.cells = 1,
  min.features = 10
)

adata.obs <- read.csv(file = "data/epi/LSI-2_RESULTS_metadata.csv")
rownames(adata.obs) <- adata.obs$X
adata.obs <- subset(adata.obs, select = -X)

czi_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = adata.obs
)

saveRDS(czi_atac, file = "data/signac/seurat_objects/210303_merged_peak_matrix.rds")




