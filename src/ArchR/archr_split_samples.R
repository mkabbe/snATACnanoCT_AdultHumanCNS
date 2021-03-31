#Run always with new R sessions
#devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.1", repos = BiocManager::repositories())
library(ArchR)
#setwd('~/Documents/PhD/snATAC-seq/trials/titration_data/snATAC_CR_out/P18861_1003')
addArchRThreads(threads = 1)

setwd("/Volumes/shared/Molekylär Neurobiologi/Castelo-Branco/Mukund/CZI_ADULT/CZI_ATAC/data/")
#install.packages("~/Downloads/BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar", repos = NULL, type = "source")
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRGenome("hg38")

#load an ArchR project
projATAC <- loadArchRProject(path = "Save-ProjATAC//", 
                 force = FALSE, 
                 showLogo = TRUE) 
# save the ArchR object
saveArchRProject(ArchRProj = projATAC, 
                 outputDirectory = "Save-ProjATAC", 
                 load = FALSE)



# Create Arrow file
REL_PATH = "../../../../NGSDATA/scATAC_human_adult_CZI_1/Processed_data/"
inputFragments <- c(paste0(REL_PATH, "P20056_1001/outs/fragments.tsv.gz"),
                    paste0(REL_PATH, "P20056_1002/outs/fragments.tsv.gz"),
                    paste0(REL_PATH, "P20056_1003/outs/fragments.tsv.gz")  ,
                    paste0(REL_PATH, "P20056_1004/outs/fragments.tsv.gz"),
                    paste0(REL_PATH, "P20057_1001/outs/fragments.tsv.gz"),
                    paste0(REL_PATH, "P20057_1002/outs/fragments.tsv.gz"),
                    paste0(REL_PATH, "P20057_1003/outs/fragments.tsv.gz"))

sampleName <- c("P20056_1001", "P20056_1002", "P20056_1003", "P20056_1004",
                "P20057_1001", "P20057_1002", "P20057_1003")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFragments,
  sampleNames = sampleName,
  filterTSS = 4,
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles <- c("P20056_1001.arrow", "P20056_1002.arrow",
                "P20056_1003.arrow", "P20056_1004.arrow",
                "P20057_1001.arrow", "P20057_1002.arrow",
                "P20057_1003.arrow")

projATAC <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "~/Downloads/PhD/CZI_analysis/ArchR_210218",
  copyArrows = FALSE
)

projATAC
getAvailableMatrices(projATAC)
head(projATAC$Sample)
quantile(projATAC$TSSEnrichment)
meta <- projATAC@cellColData
projATAC@sampleMetadata



df <- getCellColData(projATAC, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
p
p1 <- plotGroups(
  ArchRProj = projATAC, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1

p2 <- plotGroups(
  ArchRProj = projATAC, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges",
  alpha = 0.6,
  addBoxPlot = TRUE
)
p2

p3 <- plotGroups(
  ArchRProj = projATAC, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p3

p4 <- plotGroups(
  ArchRProj = projATAC, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges",
  alpha = 0.6,
  addBoxPlot = TRUE
)
p4  

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", 
        ArchRProj = projATAC, addDOC = FALSE, 
        width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = projATAC)
p1 
p2 <- plotTSSEnrichment(ArchRProj = projATAC)
p2

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", 
        ArchRProj = projATAC, addDOC = FALSE, 
        width = 5, height = 5)


addArchRThreads(threads = 1)

projATAC <- addIterativeLSI(
  ArchRProj = projATAC,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

projATAC <- addHarmony(
  ArchRProj = projATAC,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)


projATAC

projATAC <- addClusters(
  input = projATAC,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

head(projATAC$Clusters)
table(projATAC$Clusters)

cM <- confusionMatrix(paste0(projATAC$Clusters), paste0(projATAC$Sample))
cM

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black",
  cluster_cols = TRUE
)
p

projATAC <- addUMAP(
  ArchRProj = projATAC, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = projATAC, colorBy = "cellColData", 
                    name = "Sample", embedding = "UMAP")
p1
p2 <- plotEmbedding(ArchRProj = projATAC, colorBy = "cellColData", 
                    name = "Clusters", embedding = "UMAP")
p2
ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "UMAP_clusters.pdf", 
        ArchRProj = projATAC, addDOC = FALSE, 
        width = 5, height = 5)

projATAC <- addUMAP(
  ArchRProj = projATAC, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p3 <- plotEmbedding(ArchRProj = projATAC, colorBy = "cellColData", 
                    name = "Sample", embedding = "UMAPHarmony")
p3
p4 <- plotEmbedding(ArchRProj = projATAC, colorBy = "cellColData", 
                    name = "Clusters", embedding = "UMAPHarmony")
p4
ggAlignPlots(p3, p4, type = "h")

plotPDF(p3,p4, name = "UMAPHarmony_clusters.pdf", 
        ArchRProj = projATAC, addDOC = FALSE, 
        width = 5, height = 5)


markersGS <- getMarkerFeatures(
  ArchRProj = projATAC, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, 
                         cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
head(markerList$C18)

markerGenes <- c(
  "ALDH1L1", "AQP4", "S100B",  #ASTRO
  "OLIG2", "PDGFRA", "SOX10", "PTPRZ1", #OPC 
  "MAG", "MYRF", "PLP1", #OLIGO 
  "RBFOX3", "DLG4", "MAP2", #NEU
  "TMEM119", "AIF1", "CD11B", #MIGL 
  "SLC17A7", "SLC6A1" # testing
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot",
                     cluster_rows =TRUE)


projATAC <- addImputeWeights(projATAC, k=7,
                             sampleCells = 500,
                             td = 10)


p5 <- plotEmbedding(
  ArchRProj = projATAC,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(projATAC)
)

p5$AQP4


setwd("~/Downloads/PhD/CZI_analysis/ArchR_210218/")
doubScores <- addDoubletScores(
  input = projATAC,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  outDir = "~/Downloads/PhD/CZI_analysis/ArchR_210218/"
)
doubScores@cellColData

deDoubled_proj <- filterDoublets(doubScores)

de_doublet_bcd <- deDoubled_proj$cellNames
write.csv(de_doublet_bcd, file = "de_doublet_bcd.csv")
write.csv(projATAC$cellNames, file = "all_bcd.csv")

deDoubled_proj



czi_rna <- readRDS("../../../titration_ATAC/cellranger_output/S10-0050_stock/ArchR/luise_Seurat/HCA_small.rds")

czi_rna

DimPlot(czi_rna, reduction = "umap")

colnames(czi_rna@meta.data)

czi_rna@meta.data$rough_annot

projATAC <- addGeneIntegrationMatrix(
  ArchRProj = projATAC, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = czi_rna,
  addToArrow = FALSE,
  groupRNA = "rough_annot",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

getAvailableMatrices(projATAC)

p1 <- plotEmbedding(
  projATAC, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un",
  embedding = "UMAPHarmony"
  )
p1

pathToMacs2 <- "/opt/anaconda3/bin/MACS2"


projATAC <- addGroupCoverages(ArchRProj = projATAC, 
                               groupBy = "predictedGroup_Un"
							   )


							   projATAC <- addReproduciblePeakSet(
							     ArchRProj = projATAC, 
							     groupBy = "predictedGroup_Un", 
							     pathToMacs2 = pathToMacs2
							   )


cell_meta <- projATAC@cellColData
write.csv(cell_meta, file = "epi/archr_metadata.csv", sep=",")


## Took forever to run, and then failed. Issue with directory.
projATAC <- addReproduciblePeakSet(
  ArchRProj = projATAC, 
  groupBy = "predictedGroup_Un", 
  pathToMacs2 = pathToMacs2
)

					   							   
saveArchRProject(ArchRProj = projATAC, 
                 outputDirectory = "Save-ProjATAC", 
                 load = FALSE
				 )
				 
			