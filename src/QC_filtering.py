import os
import pandas as pd
import pickle 


dct = {}
min_TSS_score = 3
min_features = 1000 
min_filtered_cells = 500

qc_dict = {"Columns": ["Total_Cells", "failed_QC", "fail_percent", "is_good_sample"]}

DATA_PATH = "/data/proj/GCB_MK/10XATAC/data/scatac/cellranger_COUNT_output/"


## mark "bad cells" with low reads and low TSS scores
for file_ in os.listdir("data/SAMPLE_DATA/"):
    if file_.endswith("_TSS_obs.csv"):
        sample_ID = file_.split("_TSS")[0]
        df = pd.read_csv("data/SAMPLE_DATA/" + file_, header=0)
        bad_cells = df.loc[(df["TSS_score"] < min_TSS_score) | (df["nb_features"] < min_features), "barcode"].tolist()
        dct[sample_ID] = bad_cells
        qc_dict[sample_ID] = [len(df), len(bad_cells), round((len(bad_cells)/len(df))*100,2)]
        if len(bad_cells)/len(df) > 0.5:
            qc_dict[sample_ID].append(0)
        else:
            qc_dict[sample_ID].append(1)

## update singlecell.csv to remove the "bad cells"
for sample in dct.keys():
    df = pd.read_csv(DATA_PATH + sample + "/outs/singlecell.csv", header=0)
    for cell in dct[sample]:
        df.loc[df["barcode"] == cell, "is__cell_barcode"] = 0
    df.to_csv(DATA_PATH + sample + "/outs/singlecell.csv", header = True, index_label=False)


pickle.dump(qc_dict, open("data/SAMPLE_DATA/QC_metrics.pickle", "wb"))
