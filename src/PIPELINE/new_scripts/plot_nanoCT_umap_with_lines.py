



def plot_umaps_with_lines(
umap_coords_modality1, umap_coords_modality2, 
cluster_labels_modality, legend_dict, 
line_subset_ratio=1.0):
    # Convert the input lists to NumPy arrays and stack them to create 2D arrays
    umap_coords_modality1 = np.array(umap_coords_modality1)
    umap_coords_modality2 = np.array(umap_coords_modality2)

    # Check if the input arrays have the same number of points
    if umap_coords_modality1.shape != umap_coords_modality2.shape:
        raise ValueError("Both UMAP arrays should have the same number of points.")

    # Determine the number of points to draw lines for
    num_points = len(umap_coords_modality1)
    num_lines = int(num_points * line_subset_ratio)

    # Randomly select indices for drawing lines
    line_indices = np.random.choice(num_points, num_lines, replace=False)

    # Calculate the shift value for the second UMAP
    shift_value = 15

    # Create a new figure and set up the subplot
    plt.figure(figsize=(12, 7))

    # Remove the box around the UMAP plots
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().tick_params(axis='both', which='both', length=0)

    # Plot the first UMAP with colors based on cluster_labels_modality1
    sc_modality1 = plt.scatter(umap_coords_modality1[:, 0], umap_coords_modality1[:, 1], s=10, 
                               c=cluster_labels_modality, cmap='Paired',alpha=0.7, 
                               edgecolors='silver', linewidth=0.1, label="H3K27ac")

    # Plot the second UMAP with shifted x-coordinates and colors based on cluster_labels_modality2
    sc_modality2 = plt.scatter(umap_coords_modality2[:, 0] + shift_value, umap_coords_modality2[:, 1], s=10, 
                               c=cluster_labels_modality, cmap='Paired', alpha=0.7,
                               edgecolors='silver', linewidth=0.1, label="H3K27me3")

    # Draw lines connecting corresponding points from both plots for the selected subset
    for idx in line_indices:
        x1, y1 = umap_coords_modality1[idx]
        x2, y2 = umap_coords_modality2[idx]

        # Get the color of the dots in both scatter plots
        color_modality1 = cluster_labels_modality[idx]

        # Draw lines that connect the dots between the two modalities with the corresponding colors
        plt.plot([x1, x2 + shift_value], [y1, y2], linestyle='-', 
                 color=color_modality1, alpha=0.1)

    # Set labels 
    plt.title('H3K27ac                                                                      H3K27me3')
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    
    #Set legend
    patchList = []
    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)
    plt.legend(handles=patchList,frameon=True)
    

    # Remove the x and y tick labels
    plt.xticks([])
    plt.yticks([])

    # Show the plot
    plt.tight_layout()
    plt.show()
    
    plt.savefig("plots/230723_nanoCT_umaps_joint.png", dpi=400)
    plt.savefig("plots/230723_nanoCT_umaps_joint.pdf")
    
    
def main():
    
    import os
    import pickle
    import numpy as np
    import anndata as ad
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    #%config InlineBackend.figure_format = "retina"

    adata = ad.read("mtx/230722_P29054_H3K27ac_2kb_mtx_variable_bins_RESULTS_SHARED_CELLS.h5ad")
    medata = ad.read("mtx/230722_P29054_H3K27me3_10kb_mtx_variable_bins_RESULTS_SHARED_CELLS.h5ad")

    umap_1 = adata[adata.obs.shared_cell=="True"].obsm["X_umap"]
    umap_2 = medata[medata.obs.shared_cell=="True"].obsm["X_umap"]

    ann = list(set(adata[adata.obs.shared_cell=="True"].obs.transferred_annot))
    col = adata.uns["transferred_annot_colors"].tolist()
    col_dict = dict(zip(ann,col))
    ctype_label = [col_dict[k] for k in adata[adata.obs.shared_cell=="True"].obs.transferred_annot]


    plot_umaps_with_lines(umap_1, umap_2, ctype_label, col_dict, line_subset_ratio=0.008)

if __name__ == "__main__":
    
    main()
