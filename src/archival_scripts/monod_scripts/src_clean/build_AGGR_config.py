import os 
import pandas as pd
import pickle

## read in QC_stats file written from QC_filtering.py
df = pd.read_csv("data/SAMPLE_DATA/QC_stats_NEW_3x10000.csv")
#df = df[df.QC_sample==True]

ngi_id_list = df.NGI_ID.tolist()
cr_path = "/data/proj/GCB_MK/10XATAC/data/scatac/cellranger_COUNT_output/"
dct = {"fragments":"fragments.tsv.gz",
      "cells":"singlecell_qc_10k.csv",
      "peaks":"peaks.bed"}


config_df = pd.DataFrame(columns=["library_id","fragments","cells","peaks"])
config_df["library_id"] = ngi_id_list
for col in dct.keys():
    config_df[col] = [f"{cr_path}{x}/outs/{dct[col]}" for x in ngi_id_list]

config_df.to_csv("ref/cr_agg_config_10k.csv",sep=",",index=False, header=True)

sampleID_dict = {}
for i, sample in enumerate(ngi_id_list):
    sampleID_dict[str(i+1)] = sample

#pickle.dump(sampleID_dict, open("ref/_TEST_sampleID.pickle", "wb"))
