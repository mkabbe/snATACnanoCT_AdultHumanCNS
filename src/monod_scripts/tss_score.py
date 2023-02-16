import pickle
barcode_tss_dict = {}
in_file = "tss_overlap.bed"

with open(in_file, "r") as _handle:
    for line in _handle:
        line = line.strip().split("\t")[7:12]
        barcode = line[3]
        if barcode not in barcode_tss_dict.keys():
            barcode_tss_dict[barcode] = [line]
        else:
            barcode_tss_dict[barcode].append(line)

nb_tss_fragments = {}

for barcode in barcode_tss_dict.keys():
    nb_tss_fragments[barcode] = len(barcode_tss_dict[barcode])

with open("tss_scores.pickle", "wb") as handle:
    pickle.dump(nb_tss_fragments, handle, protocol=pickle.HIGHEST_PROTOCOL)
