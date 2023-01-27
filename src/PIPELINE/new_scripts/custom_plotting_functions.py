#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
import seaborn as sns



def convert_to_rgb(cmap):
    rgb_cmap = {}
    for var in cmap:
        color = [mcolors.to_rgb(x) for x in cmap[var].values()]
        rgb_cmap[var] = dict(zip(cmap[var].keys(), color))
    return rgb_cmap

def custom_cTop_heatmap(cistopic_obj: 'CistopicObject',
                       variables: Optional[List[str]] = None,
                       remove_nan: Optional[bool] = True,
                       scale: Optional[bool] = False,
                       cluster_topics: Optional[bool] = False,
                       color_dict: Optional[Dict[str, Dict[str, str]]] = {},
                       seed: Optional[int] = 555,
                       legend_loc_x: Optional[float] = 1.2,
                       legend_loc_y: Optional[float] = -0.5,
                       legend_dist_y: Optional[float] = -1,
                       figsize: Optional[Tuple[float, float]] = (6.4, 4.8),
                       selected_topics: Optional[List[int]] = None,
                       selected_cells: Optional[List[str]] = None,
                       harmony: Optional[bool] = False,
                       save: Optional[str] = None,
                       cmap="viridis"):

    model = cistopic_obj.selected_model
    if harmony:
        cell_topic = model.cell_topic_harmony
    else:
        cell_topic = model.cell_topic
    cell_data = cistopic_obj.cell_data

    if selected_topics is not None:
        cell_topic = cell_topic.loc[['Topic' + str(x) for x in selected_topics], ]
    if selected_cells is not None:
        cell_topic = cell_topic.loc[:, selected_cells]
        cell_data = cell_data.loc[selected_cells]

    if scale:
        cell_topic = pd.DataFrame(sklearn.preprocessing.StandardScaler().fit_transform(
            cell_topic), index=cell_topic.index.to_list(), columns=cell_topic.columns)

    if (remove_nan) & (sum(cell_data[variables].isnull().sum()) > 0):
        cell_data = cell_data[variables].dropna()
        cell_topic = cell_topic.loc[:, cell_data.index.tolist()]

    cell_topic = cell_topic.transpose()

    var = variables[0]
    var_data = cell_data.loc[:, var].sort_values()
    cell_topic = cell_topic.loc[var_data.index.to_list()]
    df = pd.concat([cell_topic, var_data], axis=1, sort=False)
    topic_order = df.groupby(var).mean().idxmax().sort_values().index.to_list()
    cell_topic = cell_topic.loc[:, topic_order].T
    # Color dict
    col_colors = {}
    if variables is not None:
        for var in variables:
            var_data = cell_data.loc[:, var].sort_values()
            categories = set(var_data)

            try:
                color = [color_dict[var][x] for x in color_dict[var].keys()]
                color = [mcolors.to_rgb(x) for x in color]
                color_dict[var] = dict(zip(categories, color))
                #color_dict = color_dictionary[var]       ## THIS LINE WILL FORCE EXCEPT BLOCK EXECUTION
                
            #except BaseException:
            except:
                #random.seed(seed)
                color = list(map(
                    lambda i: "#" + "%06x" % random.randint(0, 0xFFFFFF), range(len(categories))
                ))
                color = [mcolors.to_rgb(x) for x in color]
                color_dict[var] = dict(zip(categories, color))

            col_colors[var] = var_data.map(color_dict[var])
            
        col_colors = pd.concat([col_colors[var]
                                for var in variables], axis=1, sort=False)

        g = sns.clustermap(cell_topic,
                           row_cluster=cluster_topics,
                           col_cluster=False,
                           col_colors=col_colors,
                           #cmap=cm.viridis,
                           cmap=cmap
                           xticklabels=False,
                           figsize=figsize)

        cbar = g.cax
        cbar.set_position([legend_loc_x, 0.55, 0.05, 0.2])
        g.ax_col_dendrogram.set_visible(True)
        g.ax_row_dendrogram.set_visible(True)

        pos = legend_loc_y
        for key in color_dict:
            patchList = []
            for subkey in color_dict[key]:
                data_key = mpatches.Patch(
                    color=color_dict[key][subkey], label=subkey)
                patchList.append(data_key)
            legend = plt.legend(
                handles=patchList,
                bbox_to_anchor=(
                    legend_loc_x,
                    pos),
                loc="center",
                title=key)
            ax = plt.gca().add_artist(legend)
            pos += legend_dist_y
    else:
        g = sns.clustermap(cell_topic,
                           row_cluster=cluster_topics,
                           col_cluster=True,
#                           cmap=cm.viridis,
                           cmap=cmap,
                           xticklabels=False,
                           figsize=figsize)

    if save is not None:
        g.savefig(save, bbox_inches='tight')
    plt.show()
    

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    
    Scatter plot of log_nb_fragments vs TSS_score (density colored)
    
    source: https://stackoverflow.com/a/53865762
    
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z, **kwargs )

    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel('Density')

    return ax