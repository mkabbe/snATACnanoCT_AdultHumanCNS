import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad 

adata = ad.read("data/epi/feature_matrices/adata.h5ad")

tissue = []
cells = []
caseNO = []

NGI_ID = adata.obs["NGI_ID"].unique().to_list()
cell_counts = adata.obs["NGI_ID"].value_counts().to_dict()

for sample in NGI_ID:
    tissue.append(adata.obs.loc[adata.obs["NGI_ID"] == sample, "Tissue"].unique()[0])
    caseNO.append(adata.obs.loc[adata.obs["NGI_ID"] == sample, "caseNO"].unique()[0])
    cells.append(cell_counts[sample])

data = [cells, caseNO, tissue]
df = pd.DataFrame(data).T
df.columns=["cells", "caseNO", "tissue"]

axes = sns.scatterplot(x = df.tissue, y = df.caseNO, size = df.cells, hue = df.caseNO)