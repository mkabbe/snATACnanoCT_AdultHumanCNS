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

luise_RNA <- "../../titration_ATAC/cellranger_output/S10-0050_stock/ArchR/luise_Seurat/HCA_rough_annotated_all.RDS"
czi_atac <- readRDS("data/signac/seurat_objects/210303_merged_peak_matrix.rds")

dim(czi_atac)
czi_atac[["peaks"]]


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(czi_atac) <- annotations


# compute nucleosome signal score per cell
czi_atac <- NucleosomeSignal(object = czi_atac)
czi_atac$nucleosome_group <- ifelse(czi_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = czi_atac, group.by = 'nucleosome_group')

# compute TSS enrichment score per cell
czi_atac <- TSSEnrichment(object = czi_atac, fast = FALSE)
czi_atac$high.tss <- ifelse(czi_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(czi_atac, group.by = 'high.tss') + NoLegend()

czi_atac <- RunTFIDF(czi_atac)
czi_atac <- FindTopFeatures(czi_atac, min.cutoff = 'q0')
czi_atac <- RunSVD(czi_atac)

DepthCor(czi_atac)

czi_atac <- RunUMAP(object = czi_atac, reduction = 'lsi', dims = 2:30)
czi_atac <- FindNeighbors(object = czi_atac, reduction = 'lsi', dims = 2:30)
czi_atac <- FindClusters(object = czi_atac, verbose = FALSE, algorithm = 3)
DimPlot(object = czi_atac, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(czi_atac)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
czi_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
czi_atac <- NormalizeData(
  object = czi_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(czi_atac$nCount_RNA)
)


DefaultAssay(czi_atac) <- 'RNA'

FeaturePlot(
  object = czi_atac,
  features = c('SOX10', 'OLIG2', 'MBP', 'TBP', 'AIF1', 'GFAP'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

luise_rna <- readRDS(luise_RNA)

transfer.anchors <- FindTransferAnchors(
  reference = luise_rna,
  query = czi_atac,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = luise_rna$rough_annot,
  weight.reduction = czi_atac[['lsi']],
  dims = 2:30
)

czi_atac <- AddMetaData(object = czi_atac, metadata = predicted.labels)

plot1 <- DimPlot(object = luise_rna,
                 group.by = 'rough_annot',
                 label = TRUE,
                 repel = TRUE) + NoLegend() + ggtitle('CZI_snRNA_luise_rough_annot')

plot2 <- DimPlot(object = czi_atac, 
                 group.by = 'predicted.id',
                 label = TRUE,
                 repel = TRUE) + NoLegend() + ggtitle('CZI_snATAC_210303')

plot1 + plot2
plot2

p <- DimPlot(object = czi_atac,
             group.by = 'prediction.score.Oligo',
             label = TRUE,
             repel = TRUE)+ NoLegend() + ggtitle('LSI-2_louvain_clusters')
p
czi_atac
DimPlot(object = czi_atac, label = TRUE) + NoLegend()


colnames(czi_atac@meta.data)

saveRDS(czi_atac, file = "data/signac/seurat_objects/210303_RESULTS_merged_peak_matrix.rds")















