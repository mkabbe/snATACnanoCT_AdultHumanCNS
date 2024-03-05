import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df_csc = pd.read_csv("/data/proj/GCB_MK/scCT/rose_SE/rose/OLG_CSC_roseSE/crc_out_3/CSC_OLG_CRC_DEGREE_TABLE.txt",sep="\t")
df_ba = pd.read_csv("/data/proj/GCB_MK/scCT/rose_SE/rose/OLG_BA4_roseSE/crc_out_3/BA4_OLG_CRC_DEGREE_TABLE.txt",sep="\t")
overlap_tf = set(df_ba.Tf).intersection(set(df_csc.Tf))

df_csc = df_csc.loc[df_csc.Tf.isin(overlap_tf)]
df_ba = df_ba.loc[df_ba.Tf.isin(overlap_tf)]

df = pd.DataFrame(data=list(overlap_tf),columns=["Tf"]).sort_values(by="Tf")
df["In_CSC"] = list(df_csc["In_Degree"])
df["Out_CSC"] = list(df_csc["Out_Degree"])
df["In_Ctx"] = list(df_ba["In_Degree"])
df["Out_Ctx"] = list(df_ba["Out_Degree"])

df["Regulator_strength_CSC"] = df["In_CSC"] - df["Out_CSC"]
df["Regulator_strength_Ctx"] = df["In_Ctx"] - df["Out_Ctx"]

#df["diff_SC"] = list(df_csc["Out_Degree"] - df_csc["In_Degree"]) # +ve value means TF is a strong REGULATOR
#df["diff_Ctx"] = list(df_ba["Out_Degree"] - df_ba["In_Degree"])

df["delta_Out"] = df["Out_CSC"] - df["Out_Ctx"]
df["delta_In"] = df["In_CSC"] - df["In_Ctx"]


lim=750
threshold = 100 #strength threshold

## Mark TFs that are strong REGULATORS in one region, but weak in the other
df["sig_diff"] = df['Regulator_strength_CSC'] * df['Regulator_strength_Ctx'] < 0 ## Regulators vs Targets 
df["abs_diff"] = np.absolute(df['Regulator_strength_CSC'] - df['Regulator_strength_Ctx'])
df["abs_diff_sig"] = df["abs_diff"] > threshold ## Minimum difference in strength 
df["sig"] = df["abs_diff_sig"] & df["sig_diff"] 


fig,ax = plt.subplots(1,1,figsize=(7,7))
sns.scatterplot(x=df["Regulator_strength_Ctx"], y=df["Regulator_strength_CSC"],
                s=20, palette=["gainsboro","red"],ax=ax,ec=None,alpha=0.9,hue=df["sig"])
sns.despine(trim=False)
ax.set_xlim(-lim,lim)
ax.set_ylim(-lim,lim)
ax.axline((0, 0), slope=1, color="black", ls="-", lw=1)
ax.axline((0, threshold), slope=1, color="black", ls="--", lw=1)
ax.axline((threshold, 0), slope=1, color="black", ls="--", lw=1)

ax.axhline((0) , color="gray", ls="--", lw=1,)
ax.axvline((0), color="gray", ls="--", lw=1)
ax.legend("",frameon=False);

plt.savefig("/data/proj/GCB_MK/scCT/nanoCT_EAE/plots/CRC_TFnetwork_CSC_vs_BA4_OLG_scatter.png",dpi=600, bbox_inches="tight")