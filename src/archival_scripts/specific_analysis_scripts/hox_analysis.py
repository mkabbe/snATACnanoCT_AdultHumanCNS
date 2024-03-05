#####
#####
#   HOX ANALYSIS CODE   #



def extractPromoterCoord(gtf_file, window = 2000, feature_type = "gene", annotation = "HAVANA"):
    gtf = {}
    extnd = window//2
    with open(gtf_file) as f:
            for line in f:
                if line[0:2] != "##" and "\t"+feature_type+"\t" in line and "\t"+annotation+"\t" in line:
                    line = line.rstrip("\n").split("\t")
                    if line[6] == "+":
                        if line[0] not in gtf.keys():
                            gtf[line[0]] = [[int(line[3])-extnd, int(line[3])+extnd, line[-1].split(";")[:-1]]]
                        else:
                            gtf[line[0]].append([int(line[3])-extnd, int(line[3])+extnd, line[-1].split(";")[:-1]])
                    else:
                        if line[0] not in gtf.keys():
                            gtf[line[0]] = [[int(line[4])-extnd, int(line[4])+extnd, line[-1].split(";")[:-1]]]
                        else:
                            gtf[line[0]].append([int(line[4])-extnd, int(line[4])+extnd, line[-1].split(";")[:-1]])
    return gtf


os.system("cat ref/gencode.v35.annotation.gtf | grep HOX > ref/hox_gencode.gtf")
hox_dict = extractPromoterCoord("ref/hox_gencode.gtf")

hox_genes = []
gene_index = []
for chrom in hox_dict.keys():
    for gene in hox_dict[chrom]:
        gene_start = gene[0]
        gene_end = gene[1]
        gene_name = gene[-1][2].lstrip(' gene_name "').rstrip('"')
        hox_genes.append([chrom, gene_start, gene_end, gene_name])

hoxa_genes = []
hoxb_genes = []
hoxc_genes = []
hoxd_genes = []
for x in hox_info:
    if x[-1][:4]=="HOXA":
        hoxa_genes.append(x)
    elif x[-1][:4]=="HOXB":
        hoxb_genes.append(x)
    elif x[-1][:4]=="HOXC":
        hoxc_genes.append(x)
    elif x[-1][:4]=="HOXD":
        hoxd_genes.append(x)

hoxa_df = pd.DataFrame.from_records(hoxa_genes)
hoxb_df = pd.DataFrame.from_records(hoxb_genes)
hoxc_df = pd.DataFrame.from_records(hoxc_genes)
hoxd_df = pd.DataFrame.from_records(hoxd_genes)

hoxa_df.to_csv("data/cluster_fragments/OL_typeTissue/hox_regions/hoxa_promoters.bed", header=False, index=False, sep="\t")
hoxb_df.to_csv("data/cluster_fragments/OL_typeTissue/hox_regions/hoxb_promoters.bed", header=False, index=False, sep="\t")
hoxc_df.to_csv("data/cluster_fragments/OL_typeTissue/hox_regions/hoxc_promoters.bed", header=False, index=False, sep="\t")
hoxd_df.to_csv("data/cluster_fragments/OL_typeTissue/hox_regions/hoxd_promoters.bed", header=False, index=False, sep="\t")