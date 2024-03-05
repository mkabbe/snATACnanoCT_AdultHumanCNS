import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pickle
import plotly.graph_objects as go
import ast

os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE/")


def flatten(xss):
    return [x for xs in xss for x in xs]

rnadata = ad.read("../../CZI/episcanpy_analysis/AGG_ATAC_210218/data/feature_matrices/HCA_RNA_all_annot.h5ad")
rnadata.var_names = rnadata.var.gene_names
rnadata = rnadata[rnadata.obs.Tissue!="CB"]

cluster_dict = {"AST":['Astrocyte_1','Astrocyte_2','Astrocyte_3','Astrocyte_4'],
               "MOL":["Oligo"],
               "OPC":["OPC"],
               "MIGL":["Microglia-Macrophages_1","Microglia-Macrophages_2"],
               "CXEX":["Neuron_Ex_1","Neuron_Ex_2","Neuron_Ex_3"],
               "CXINH":["Neuron_In_1","Neuron_In_2"]}
                   
c1 = sys.argv[1]


df = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/{c1}_roseSE/crc_out_3/{c1}_CRC_DEGREE_TABLE.txt",sep="\t")

df[f"TF_strength_{c1}"] = df[f"In_Degree"] - df[f"Out_Degree"]
df = df.sort_values(by=f"TF_strength_{c1}",ascending=False)
df["TF_rank"] = [x for x in range(1,len(df)+1)]

c=c1

df = pd.read_csv(f"/data/proj/GCB_MK/scCT/rose_SE/rose/{c1}_roseSE/crc_out_3/{c1}_CRC_DEGREE_TABLE.txt",sep="\t")
df[f"TF_strength_{c1}"] = df[f"In_Degree"] - df[f"Out_Degree"]
df = df.sort_values(by=f"TF_strength_{c1}",ascending=False)
df["TF_rank"] = [x for x in range(1,len(df)+1)]
df["isHox"] = [True if x.startswith("HOX") else False for x in df.Tf]

nodes = pd.read_csv(f"/date/gcb/GCB_MK/CRC_output/{c}_CRC_NODELIST.txt",header=None)
node_list = [x for x in nodes[0] if x in set(rnadata.var_names)]

subset = rnadata[rnadata.obs.clusters_named.isin(cluster_dict[c]),rnadata.var_names.isin(node_list)]
node_list = [x for x in subset.var.gene_names] # preserve order


df2 = pd.DataFrame()
df2["nodes"] = node_list
df2[f"{c}_mean_gex"] = list(np.mean(subset.X.toarray(),axis=0))

gex_dict = dict(zip(df2.nodes,df2[f"{c}_mean_gex"]))

df["mean_gex"] = [gex_dict[x] if x in gex_dict else 0 for x in df.Tf]
df["wTF_strength"] = df[f"TF_strength_{c1}"] * df["mean_gex"]


df = df.sort_values(by="wTF_strength",ascending=False)


with open(f"/data/proj/GCB_MK/scCT/rose_SE/rose/{c1}_roseSE/crc_out_3/{c1}_CRC.ntx", "rb") as ntxfile:
        network_dict_of_lists = pickle.load(ntxfile)
G = nx.from_dict_of_lists(network_dict_of_lists)

n_nodes = 15
G = nx.subgraph(G, df.head(n_nodes).Tf)
#G = nx.subgraph(G, clique_nodes)
layout = nx.spring_layout(G)

## Edges ##
edge_x = []
edge_y = []

for reporter, partner in G.edges():
    x0, y0 = layout[reporter]
    x1, y1 = layout[partner]

    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None)
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

edge_trace = go.Scatter(
    x=edge_x, y=edge_y,
    line=dict(width=1, color='silver'),
    hoverinfo='none',
    mode='lines')

## Nodes ##
df_node_attributes = pd.DataFrame([{'TF': node, **G.nodes[node]} for node in G.nodes()])

node_x = []
node_y = []
node_text = []

for node in df_node_attributes['TF']:
    x, y = layout[node]
    text = node
    node_x.append(x)
    node_y.append(y)
    node_text.append(text)

node_trace = go.Scatter(
    x=node_x, 
    y=node_y,
    text=node_text,
    mode='text',
    hoverinfo=None,
    textfont=dict(
        size=40, color="black"  # Set the text size
    ),
    marker=dict(
        showscale=True,
        # colorscale options
        #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
        #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
        #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
        colorscale='Reds',
        reversescale=False,
        color=[],
        size=10,
        colorbar=dict(
            thickness=10,
            title='Node Connections',
            xanchor='left',
            titleside='right'
        ),
        line_width=5
        ))
        
# Create Network Graph

fig = go.Figure(layout=go.Layout(
                title='',
                showlegend=False, 
                hovermode=None,
                margin=dict(b=20,l=5,r=5,t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )

fig.add_trace(edge_trace)
fig.add_trace(node_trace)
             
fig.update_layout(
    autosize=False,
    width=1500,
    height=1000,
    plot_bgcolor='white'
)
       
fig.show()

