#this is specific to one of your K27ac objects. Maybe you have to change some layers.

library(dplyr)
library(reticulate)
library(anndata)
library(Signac)
library(EnsDb.Hsapiens.v86)

setwd("/data/proj/GCB_MK/scCT/nanoCT_EAE/mtx/")

ad_file = "230612_P29054_H3K27ac_10kb_mtx_RESULTS_SHARED_CELLS.h5ad"
rds_file = "230612_P29054_H3K27ac_10kb_mtx_RESULTS_SHARED_CELLS.rds"

#K27ac_5kb <- read_h5ad("/data/proj/GCB_MK/scCT/nanoCT_EAE/mtx/221028_P27505_H3K27ac_5kb_mtx_RESULTS.h5ad")
K27ac_5kb <- read_h5ad(ad_file)
 
K27ac_5kb.mat <- t(as.matrix(K27ac_5kb$layers[["raw_counts"]]))
 
K27ac_5kb.annot <- cbind.data.frame(K27ac_5kb$obs[["barcode"]] , K27ac_5kb$obs[["transferred_annot"]] ,  K27ac_5kb$obs[["caseNO"]] )
 
#this depends on the features format ( these were peaks)
new.rownames <- rownames(K27ac_5kb.mat)
new.rownames <- gsub("(.*)[_](.*)[_](.*)", "\\1:\\2-\\3", (new.rownames) )
 
rownames(K27ac_5kb.mat) <- new.rownames
 
###############################
 
##no need this part
#grange.counts <- StringToGRanges(rownames(new_mat), sep = c(":", "-"))
#grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
#K27ac_5kb.counts <- new_mat[as.vector(grange.use), ]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "hg38"
seqlevelsStyle(annotations) <- 'UCSC'
 
################################
frag.file <- "fragments.tsv.gz"
 
 new_mat <- K27ac_5kb.mat %>% as.data.frame()  %>% replace(is.na(.), 0) #%>% as.numeric()
 
chrom_assay <- CreateChromatinAssay(
   counts = new_mat ,
   sep = c(":", "-"),
  # genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
 
HS_CT_MK.seur <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "K27ac_10Kb_bins",
  meta.data = K27ac_5kb.annot
)

saveRDS(HS_CT_MK.seur, rds_file)
