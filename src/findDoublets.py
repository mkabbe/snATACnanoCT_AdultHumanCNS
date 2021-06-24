import os
import time
import logging as logg
import anndata as ad
import episcanpy.api as epi
import concurrent.futures
import scrublet as scr


logg.basicConfig(filename='logs/Scrublet.log', format='%(asctime)s: %(message)s', level=logging.DEBUG)

datapath = "/data/proj/GCB_MK/10XATAC/data/scatac/cellranger_COUNT_output"
#annot = epi.ct.load_features("ref/2kb_binned_hg38.bed")
annot = epi.ct.load_features("data/peaks/210510_merged_louvain_cluster_peaks.bed")


samples = []
for fn in os.listdir(datapath):
        samples.append(f"{datapath}/{fn}/outs")


def extract(sample, anndata = False):
    ## returns fragments, singlecell, and outfile name
    frag = sample+"/fragments.tsv.gz"
    csv = sample+"/singlecell_EPI.csv"
    save = sample+"/peak_mtx.h5ad"
    if anndata:
        return(save)
    return(frag,csv,save)

def build_mtx(sample_id):
    logg.info(f"***Starting {sample_id.split('/')[-2]}")
    f,c,s = extract(sample_id)
    epi.ct.bld_mtx_fly(tsv_file=f, 
                         csv_file=c, 
                         annotation=annot, 
                         save=s)
    logg.info(f"***Finished {sample_id.split('/')[-2]}")

def runScrublet(sample_id):
    ann = sample_id+"/peak_mtx.h5ad"
    adata = ad.read(ann)
    logg.info(f"***Running scrublet for {sample_id.split('/')[-2]}")
    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(mean_center=False)
    logg.info(f"***Finished running scrublet for {sample_id.split('/')[-2]}")
    adata.obs["doublet_scores"] = doublet_scores.tolist()
    adata.obs["predicted_doublets"] = predicted_doublets.tolist()
    adata.obs.to_csv(sample_id+"/doublets.csv")
    logg.info(f"***Saved doublet scores for {sample_id.split('/')[-2]}")


start = time.perf_counter()
with concurrent.futures.ThreadPoolExecutor() as executor:
    executor.map(build_mtx, samples)
print(f"BUILD_MATRIX TOTAL TIME: {round((time.perf_counter()-start),2)} seconds")

start = time.perf_counter()
for s in samples:
    runScrublet(s)
print(f"SCRUBLET TOTAL TIME: {round((time.perf_counter()-start),2)} seconds")