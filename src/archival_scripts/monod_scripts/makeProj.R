library(ArchR)
addArchRThreads(threads = 20)

setwd("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/Archr_files")
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRGenome("hg38")

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

ArrowFiles <- c()
for(i in 1:length(sampleName)){
  arrow <- paste0(sampleName[i],".arrow")
  ArrowFiles <- c(ArrowFiles, arrow)
}


projATAC <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Projects",
  copyArrows = FALSE
)

projATAC
saveArchRProject(ArchRProj = projATAC, outputDirectory = "../data/archr/Save-ProjATAC1", load = FALSE)
