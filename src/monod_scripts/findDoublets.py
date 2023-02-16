import os
import time
import logging as logg
import concurrent.futures
import anndata as ad
import pandas as pd
import scrublet as scr


logg.basicConfig(filename='logs/Scrublet.log', format='%(asctime)s: %(message)s', level=logg.DEBUG)

datapath = "/data/proj/GCB_MK/10XATAC/data/scatac/cellranger_COUNT_output"

samples = []
for fn in os.listdir(datapath):
    if not fn.startswith("P2005"):
        samples.append(f"{datapath}/{fn}/outs")

def runScrublet(sample_id):
    #ann = extract(sample_id, anndata=True)
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
#with concurrent.futures.ProcessPoolExecutor() as executor:
#    executor.map(runScrublet, samples)

for s in samples:
    runScrublet(s)

print(f"{round((time.perf_counter()-start),2)} seconds")
