library(Signac)
library(Seurat)
library(dplyr)
library(reticulate)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(rtracklayer)
library(Matrix)

setwd("/Volumes/projects/Castelo_Branco/Mukund/CZI_ADULT/CZI_ATAC/data/AGG_ATAC_210507")

FRAGMENTS_path <- "outs/fragments.tsv.gz"
CELLS <- read.csv("outs/OLG_barcodes.csv", header=F)
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

#writeMM(count.matrix.peaks, file = "count_peak_matrix.mtx")
#count.matrix.peaks <- readMM(file = "count_peak_matrix.mtx")

chrom_assay <- CreateChromatinAssay(
  counts = count.matrix.peaks,
  sep = c(":","-"),
  genome = "hg38",
  fragments = FRAGMENTS,
  min.cells = 1,
  min.features = 10
)

adata.obs <- read.csv(file = "outs/OLG_metadata_for_LSI2.csv")
obs_df <- data.frame(row.names = "X", adata.obs)

ol_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = obs_df
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(ol_atac) <- annotations

## SAVE
saveRDS(ol_atac, file = "signac_out/210624_olg_lsi2_merged_peak_matrix.rds")

### Normalization and DimRed
ol_atac <- RunTFIDF(ol_atac)
ol_atac <- FindTopFeatures(ol_atac, min.cutoff = 'q0')
ol_atac <- RunSVD(object = ol_atac)
ol_atac <- RunUMAP(object = ol_atac,reduction = 'lsi',dims = 2:30)
ol_atac <- FindNeighbors(object = ol_atac,reduction = 'lsi',dims = 2:30)
ol_atac <- FindClusters(object = ol_atac,algorithm = 3,
                        resolution = 1.2,verbose = FALSE)
DimPlot(object = ol_atac, label = TRUE, group.by="caseNO")

### GENE ACTIVITY MATRIX 
gene.activities <- GeneActivity(ol_atac)
# add the gene activity matrix to the Seurat object as a new assay
ol_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
ol_atac <- NormalizeData(object = ol_atac, assay = 'RNA',
                         normalization.method = 'LogNormalize',
                         scale.factor = median(ol_atac$nCount_RNA))

DefaultAssay(ol_atac) <- 'RNA'
FeaturePlot(object = ol_atac,
            features = c('PDGFRA','SOX10',"CSPG4","PLP1","CNP","GAPDH"),
            pt.size = 0.1, max.cutoff = 'q95', ncol = 3)

### Label transfer
luise_rna <- readRDS("/Users/mukkab/Downloads/srt_oligos_and_opcs_LS.RDS")

transfer.anchors <- FindTransferAnchors(
  reference = luise_rna,
  query = ol_atac,
  reduction = 'cca',
  dims = 1:40
)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = luise_rna$subclass,
  weight.reduction = ol_atac[['lsi']],
  dims = 2:30
)
ol_adata <- AddMetaData(object = ol_adata, metadata = predicted.labels)
