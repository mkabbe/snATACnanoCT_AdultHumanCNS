
# Convert to obs/exp

#hicTransform -m P27454_hOPC_REdeep/P27454_REdeep_hOPC.mapped.pairs.contact-map.5000.cool --method obs_exp --chromosomes chr17 -o obs_exp_hOPC_chr17_5000.cool
#hicTransform -m P27454_hOPC_REdeep/P27454_REdeep_hOPC.mapped.pairs.contact-map.5000.cool --method obs_exp --chromosomes chr2 -o obs_exp_hOPC_chr2_5000.cool
#hicTransform -m P27454_hOPC_REdeep/P27454_REdeep_hOPC.mapped.pairs.contact-map.5000.cool --method obs_exp --chromosomes chr7 -o obs_exp_hOPC_chr7_5000.cool
#hicTransform -m cools/hOPC_100kb.cool --method obs_exp -o obs_exp_hOPC_100kb.cool
#hicTransform -m cools/bCell3_100kb.cool --method obs_exp -o cools/obs_exp_bCell3_100kb.cool
#hicTransform -m cools/bCell1_100kb.cool --method obs_exp -o cools/obs_exp_bCell1_100kb.cool
#hicTransform -m cools/bCell2_100kb.cool --method obs_exp -o cools/obs_exp_bCell2_100kb.cool
#hicTransform -m cools/bCell_merged_100kb.cool --method obs_exp -o cools/obs_exp_bCell_merged_100kb.cool
#hicTransform -m cools/bCell_merged_5kb.cool --method obs_exp -o cools/obs_exp_bCell_merged_5kb.cool

## Get eigen vectors 

#hicPCA -m cools/hOPC_100kb.cool -o hOPC_100kb_pca1.bigwig hOPC_100kb_pca2.bigwig -f bigwig
#hicPCA -m cools/hOPC_100kb.cool -o hOPC_100kb_pca1.bedgraph hOPC_100kb_pca2.bedgraph -f bedgraph
#hicPCA -m cools/hOPC_10kb.cool -o hOPC_10kb_pca1.bedgraph hOPC_10kb_pca2.bedgraph -f bedgraph

#hicPCA -m cools/bCell1_100kb.cool -o bCell_1_100kb_pca1.bigwig bCell_1_100kb_pca2.bigwig -f bigwig
#hicPCA -m cools/bCell2_100kb.cool -o bCell_2_100kb_pca1.bigwig bCell_2_100kb_pca2.bigwig -f bigwig
#hicPCA -m cools/bCell3_100kb.cool -o bCell_3_100kb_pca1.bigwig bCell_3_100kb_pca2.bigwig -f bigwig


## Merge matrices

#hicSumMatrices -m P27454_1002/P27454_1002.PT.pairs.contact-map.100000.unbal.cool P27454_1001/P27454_1001.PT.pairs.contact-map.100000.unbal.cool P27454_1003/P27454_1003.PT.pairs.contact-map.100000.unbal.cool -o cools/bCell_merged_100kb.cool
#cooler balance cools/bCell_merged_100kb.cool

#hicSumMatrices -m P27454_1002/P27454_1002.PT.pairs.contact-map.5000.unbal.cool P27454_1001/P27454_1001.PT.pairs.contact-map.5000.unbal.cool P27454_1003/P27454_1003.PT.pairs.contact-map.5000.unbal.cool -o cools/bCell_merged_5kb.cool
#cooler balance cools/bCell_merged_5kb.cool
#hicTransform -m cools/bCell_merged_5kb.cool --method obs_exp -o cools/obs_exp_bCell_merged_5kb.cool

hicSumMatrices -m P27454_1002/P27454_1002.PT.pairs.contact-map.1000.unbal.cool P27454_1001/P27454_1001.PT.pairs.contact-map.1000.unbal.cool P27454_1003/P27454_1003.PT.pairs.contact-map.1000.unbal.cool -o cools/bCell_merged_1kb.cool
cooler balance cools/bCell_merged_1kb.cool
#hicTransform -m cools/bCell_merged_1kb.cool --method obs_exp -o cools/obs_exp_bCell_merged_1kb.cool

## Differential loops

#/date/gcb/GCB_MK/microC/mustache_hic/mustache/mustache/diff_mustache.py \
#-f1 P27454_1001/P27454_1001.PT.pairs.contact-map.5000.cool \
#-f2 P27454_hOPC_REdeep/P27454_REdeep_hOPC.mapped.pairs.contact-map.5000.cool \
#-o DiffLoops/diffLoops_1001_v_hOPC.tsv \
#-r 5000 \
#-cz /data/proj/GCB_MK/micro-C/hg38.valid.chrom.sizes \
#-pt 0.01 \
#-pt2 0.01 \
#-p 8
#--normalization KR 



#/date/gcb/GCB_MK/microC/mustache_hic/mustache/mustache/diff_mustache.py \
#-f1 P27454_1003/P27454_1003.PT.pairs.contact-map.5000.cool \
#-f2 P27454_1001/P27454_1001.PT.pairs.contact-map.5000.cool \
#-o DiffLoops/diffLoops_1003_v_1001.tsv \
#-r 5000 \
#-cz /data/proj/GCB_MK/micro-C/hg38.valid.chrom.sizes \
#-pt 0.01 \
#-pt2 0.01 \
##-p 8
#--normalization KR
