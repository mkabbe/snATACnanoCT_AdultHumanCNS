

import pandas as pd

project = "P27505"
modality = "H3K27ac"
#modality = "H3K27me3"
df_agg = pd.read_csv(f"/proj/tmp/mukund/ct_eae/{project}_{modality}_agg_relax/outs/singlecell_EPI.csv")
df_agg["agg_sample"] = [x.split("-")[-1] for x in df_agg.barcode]
df_agg["bcd_nosuff"] = [x.split("-")[0] for x in df_agg.barcode]

for sample in ["1001","1002","1003","1004"]:
    df_1 = pd.read_csv(f"singlecell_CT_filtered/{project}_{sample}_{modality}_singlecell_CT_filtered_relaxed.csv")
    df_1["bcd_nosuff"] = [x.split("-")[0] for x in df_1.barcode]
    unique_set = set(df_1.bcd_nosuff)
    ## Set is_cell_barcode for all custom selected barcodes to 1
    df_agg.loc[(df_agg.agg_sample=="1") & (df_agg.bcd_nosuff.isin(unique_set)),"is_cell_barcode"] = 1

# write to new singlecell file. Use this to build the matrix.
df_agg.loc[df_agg.is_cell_barcode==1].to_csv(f"/proj/tmp/mukund/ct_eae/{project}_{modality}_agg_relax/outs/singlecell_EPI_custom_selection.csv")
