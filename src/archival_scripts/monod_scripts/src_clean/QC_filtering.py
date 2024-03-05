import os
import sys
import numpy as np
import pandas as pd
import pickle 

dct = {}
min_TSS_score = 3
min_features = 10000
min_filtered_cells = 10

qc_dict = {"NGI_ID":[], "total_cells":[], "QC_cells":[], "pass_percentage":[], "QC_sample":[]}

DATA_PATH = "/data/proj/GCB_MK/10XATAC/data/scatac/cellranger_COUNT_output/"


## mark "bad cells" with low reads and low TSS scores
for file_ in os.listdir("data/SAMPLE_DATA/"):
    if file_.endswith("_TSS_obs.csv"):
        sample_ID = file_.split("_TSS")[0]
        df = pd.read_csv("data/SAMPLE_DATA/" + file_, header=0)
        bad_cells = df.loc[(df["TSS_score"] < min_TSS_score) | (df["nb_features"] < min_features), "barcode"].tolist()
        dct[sample_ID] = bad_cells
        pass_ratio = (len(df)-len(bad_cells))/len(df)
        qc_dict["NGI_ID"].append(sample_ID)
        qc_dict["total_cells"].append(len(df))
        qc_dict["QC_cells"].append(len(df)-len(bad_cells))
        qc_dict["pass_percentage"].append(round(pass_ratio*100,2))
        qc_dict["QC_sample"].append(pass_ratio > 0.5)

## update singlecell.csv to remove the "bad cells"
for sample in dct.keys():
    df = pd.read_csv(DATA_PATH + sample + "/outs/singlecell.csv", header=0)
    for cell in dct[sample]:
        df.loc[df["barcode"] == cell, "is__cell_barcode"] = 0
    df.to_csv(DATA_PATH + sample + "/outs/singlecell_qc_10k.csv", header = True, index_label=False)        


## write QC stats to csv and pickle
qc_metrics_df = pd.DataFrame.from_dict(qc_dict)
qc_metrics_df.to_csv("data/SAMPLE_DATA/QC_stats_NEW_3x10000.csv", index=False, header=True)
pickle.dump(qc_dict, open("data/SAMPLE_DATA/QC_stats_NEW_3x10000.pickle", "wb"))
