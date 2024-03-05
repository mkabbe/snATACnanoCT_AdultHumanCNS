

def plot_signal_enrichment_profile(df, lw=3, a=1):
    fig, ax = plt.subplots(1,4, figsize=[11,2])

    legend_dict = {"H3K27ac":"green","H3K27me3":"red"}
    title_dict = {0:"H3K27me3", 1:"H3K27ac", 2:"ATAC", 3:"Nott et.al. OLG H3K27ac"}

    

    #H3K27me3
    sns.lineplot(data = df, x = "idx", y = "ac_at_me3", linewidth=lw, color="green", alpha=a,
                 ax=ax[0])
    sns.lineplot(data = df, x = "idx", y = "me3_at_me3", linewidth=lw, color="red", alpha=a,
                 ax=ax[0])

    #H3K27ac
    sns.lineplot(data = df, x = "idx", y = "ac_at_ac", linewidth=lw, color="green", alpha=a, 
                 ax=ax[1])
    sns.lineplot(data = df, x = "idx", y = "me3_at_ac", linewidth=lw, color="red", alpha=a,
                 ax=ax[1])

    #ATAC
    sns.lineplot(data = df, x = "idx", y = "ac_at_atac", linewidth=lw, color="green", alpha=a,
                 ax=ax[2])
    sns.lineplot(data = df, x = "idx", y = "me3_at_atac", linewidth=lw, color="red", alpha=a,
                 ax=ax[2])

    #Nott_OLG_H3K27ac
    sns.lineplot(data = df, x = "idx", y = "ac_at_nottolgac", linewidth=lw, color="green", alpha=a,
                 ax=ax[3])
    sns.lineplot(data = df, x = "idx", y = "me3_at_nottolgac", linewidth=lw, color="red", alpha=a,
                 ax=ax[3])


    for k in range(4):
        ax[k].set_title(title_dict[k],fontsize=8)
        ax[k].set_xlabel('Distance (in kb)', fontsize=8)
        if k==0:
            ax[k].set_ylabel('Normalized signal', fontsize=8)
        else:
            ax[k].set_ylabel('')
            

    #sns.despine(trim=True)
    sns.despine()

    # Set labels 
    #ax[0].axes.set_title('Signal at H3K27ac peaks')



    #Set legend
    patchList = []
    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)
    plt.legend(handles=patchList,frameon=False, bbox_to_anchor=(1,0.5,0.5,0.3),fontsize=7)
    plt.savefig("plots/230723_nanoCT_signal_enrichment.png", dpi=400, bbox_inches='tight')
    plt.savefig("plots/230723_nanoCT_signal_enrichment.pdf", bbox_inches='tight')
    plt.tight_layout()

# Remove the x and y tick labels
#plt.xticks([]);
#plt.yticks([])

def main():
    import os
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns
    os.chdir("/data/proj/GCB_MK/scCT/nanoCT_EAE")
    
    ls = ["ac_at_ac","me3_at_ac",
          "ac_at_me3","me3_at_me3", 
          "ac_at_atac", "me3_at_atac", 
          "ac_at_nottolgac", "me3_at_nottolgac"]

    dct = {}    
    for l in ls:
        dct[l] = pd.read_csv(f"fragments/P27505_OLG/{l}_profile.csv",sep="\t")[1:].T[2:][1]
    
    data = pd.DataFrame.from_dict(dct)
    data["idx"] = [x/100 for x in range(int(-len(data)/2),int(len(data)/2))]
    
    plot_signal_enrichment_profile(data, lw=2)
    

if __name__ =="__main__":
    main()