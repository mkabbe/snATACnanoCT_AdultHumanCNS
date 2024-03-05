library(Signac)
library(Seurat)
#BiocManager::install("JASPAR2020")
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)

setwd("/Volumes/shared/Molekylär Neurobiologi/Castelo-Branco/Mukund/CZI_ADULT/CZI_ATAC/")

czi_atac <- readRDS("data/signac/seurat_objects/210303_RESULTS_merged_peak_matrix.rds")
czi_atac
DefaultAssay(czi_atac) <- 'peaks'

p1 <- DimPlot(czi_atac, label = TRUE, pt.size = 0.1) + NoLegend()
p1

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
czi_atac <- AddMotifs(
  object = czi_atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

DefaultAssay(czi_atac) <- 'peaks'
# differential peaks between OPC vs rest
da_peaks <- FindMarkers(
  object = czi_atac,
  group.by = "predicted.id",
  ident.1 = 'Neuron_Ex',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

# test enrichment
enriched.motifs <- FindMotifs(
  object = czi_atac,
  features = top.da.peak
)
enriched_motifs_Neuron_Ex <- enriched.motifs

MotifPlot(
  object = czi_atac,
  motifs = head(rownames(enriched.motifs))
)

czi_atac <- RunChromVAR(
  object = czi_atac,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(czi_atac) <- 'chromvar'

saveRDS(czi_atac, file = "data/signac/seurat_objects/210323_RESULTS_merged_peak_matrix.rds")


# look at the activity of different motifs



p2 <- FeaturePlot(
  object = czi_atac,
  features = "MA1513.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p2
plot2

chromVar_mtx <- czi_atac[["chromvar"]]@data

write.csv(x = enriched_motifs_Oligos_bad, file = "data/signac/enriched_motifs_Oligos-bad.csv")
