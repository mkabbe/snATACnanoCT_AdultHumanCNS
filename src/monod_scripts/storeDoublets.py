import os 
import anndata as ad
import pickle
import csv

ngi_to_id = {v:k for k,v in pickle.load(open("ref/_sampleID.pickle", "rb")).items()}

datapath = "/data/proj/GCB_MK/10XATAC/data/scatac/cellranger_COUNT_output"

scores_csv = []
for fn in os.listdir(datapath):
    scores_csv.append(f"{datapath}/{fn}/outs/doublets.csv")

doublet_dict = {}

for csvfile in scores_csv:
    with open(csvfile, "r") as handle:
        ngi_id = csvfile.split("/")[-3]
        reader = csv.reader(handle)
        for row in reader:
            if row[0] == "": ## skip header
                continue
            bcd,score,is_doub = row
            bcd = bcd.split("-")[0]+"-"+ngi_to_id[ngi_id]
            doublet_dict[bcd] = [score, is_doub]

with open("ref/doublet_scores_lsi2_cells.pickle", "wb") as f:
    pickle.dump(doublet_dict, f)
