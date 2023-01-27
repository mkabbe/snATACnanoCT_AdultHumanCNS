library(Signac)
library(Seurat)
library(dplyr)
library(reticulate)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(rtracklayer)

setwd("/Volumes/projects/Castelo_Branco/Mukund/CZI_ADULT/CZI_ATAC/data/AGG_ATAC_210507")

FRAGMENTS_path <- "outs/fragments.tsv.gz"
CELLS <- read.csv("outs/barcodes_for_LSI-2.csv", header=F)
FEATURES_bed <- "outs/210510_merged_louvain_cluster_peaks.bed"
FEATURES <- rtracklayer::import(FEATURES_bed)

cell_bc = c()
for (i in CELLS) {
  cell_bc <- append(cell_bc, i)
}

FRAGMENTS <- CreateFragmentObject(FRAGMENTS_path, cells= cell_bc)

count.matrix.peaks <- FeatureMatrix(fragments = FRAGMENTS,
                                    features = FEATURES,
                                    cells = cell_bc)
#library(Matrix)
#writeMM(count.matrix.peaks, file = "count_peak_matrix.mtx")

chrom_assay <- CreateChromatinAssay(
  counts = count.matrix.peaks,
  sep = c(":","-"),
  genome = "hg19",
  fragments = FRAGMENTS,
  min.cells = 1,
  min.features = 10
)

adata.obs <- read.csv(file = "outs/metadata_for_LSI2.csv")
obs_df <- data.frame(row.names = "X", adata.obs)

czi_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = obs_df
)

saveRDS(czi_atac, file = "signac_out/210609_lsi2_merged_peak_matrix.rds")