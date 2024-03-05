import pandas as pd

df = pd.read_csv("sample_randomizer_metadata.csv")

# Sample processing date and sequencing batch info
df["10X_BATCH"] = [x.split("M")[0] if type(x)==str else "-" for x in df["MK_ID"]]
df["NGS_BATCH"] = [x.split("_")[0] if type(x)==str else "-" for x in df["NGI_ID"]]

# Remove samples that weren't processed
df = df[df["10X_BATCH"]!="-"]

# split ages into different bins
df['age_bins'] = pd.cut(x=df['Age'],
                right=False, bins=[25, 35, 45, 55, 65, 75],
                labels=["25-34","35-44","45-54","55-64","65-74"])

#store as categorical
for col in ["Age","PMI","Box"]:
    df[col] = df[col].astype(int).astype("category")

#save file
df.to_csv("sample_metadata_COMPLETE.csv", index=False)
