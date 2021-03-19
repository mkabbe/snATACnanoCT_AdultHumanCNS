import os 
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# read in matrix
adata = ad.read("data/epi/feature_matrices/210302_LSI-2_RESULTS_merged_peak_matrix.h5ad")


## stacked bar plots
p1 = sns.displot(adata.obs, x="louvain", hue="Tissue", multiple="stack")
p2 = sns.displot(adata.obs, x="louvain", hue="NGI_ID", multiple="stack")
p3 = sns.displot(adata.obs, x="louvain", hue="caseNO", multiple="stack")


# pie chart
labels = ['CB', 'BA4', 'CSC']
sizes = [adata.obs['Tissue'].value_counts()[0],
         adata.obs['Tissue'].value_counts()[1],
         adata.obs['Tissue'].value_counts()[2]]
# print(sizes) # adds up to 1433, which is the total number of participants
fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=False)
ax1.axis('equal')
plt.show()

