#devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.1", repos = BiocManager::repositories())
#Run always with new R sessions
library(ArchR)
addArchRThreads(threads = 20)

setwd("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/Archr_files")
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRGenome("hg38")

# Create Arrow files
base_path = "/data/proj/GCB_MK/10XATAC/data/scatac/cellranger_COUNT_output/"

inputFragments <- c(paste0(base_path,"P20056_1001/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20056_1002/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20056_1003/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20056_1004/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20057_1001/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20057_1002/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20057_1003/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20602_1001/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20602_1002/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20602_1003/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20602_1004/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20602_1005/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20602_1006/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20602_1007/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20603_1001/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20603_1002/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20603_1003/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20603_1004/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20603_1005/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20603_1006/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20603_1007/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1001/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1002/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1003/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1004/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1005/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1006/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1007/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1008/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1009/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1010/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1011/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1012/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1013/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1014/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1015/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1016/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1017/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1018/outs/fragments.tsv.gz"),
                    paste0(base_path,"P20908_1019/outs/fragments.tsv.gz"))

sampleName <- c("P20056_1001", "P20056_1002", "P20056_1003", "P20056_1004",
                "P20057_1001", "P20057_1002", "P20057_1003",
                "P20602_1001", "P20602_1002", "P20602_1003", "P20602_1004",
                "P20602_1005", "P20602_1006", "P20602_1007", 
                "P20603_1001", "P20603_1002", "P20603_1003", "P20603_1004", 
                "P20603_1005", "P20603_1006", "P20603_1007",
                "P20908_1001", "P20908_1002", "P20908_1003", "P20908_1004",
                "P20908_1005", "P20908_1006", "P20908_1007", "P20908_1008",
                "P20908_1009", "P20908_1010", "P20908_1011", "P20908_1012",
                "P20908_1013", "P20908_1014", "P20908_1015", "P20908_1016",
                "P20908_1017", "P20908_1018", "P20908_1019")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFragments,
  sampleNames = sampleName,
  filterTSS = 3,
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, 
    knnMethod = "UMAP",
    LSIMethod = 1
)

projATAC <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Projects",
  copyArrows = TRUE
)
