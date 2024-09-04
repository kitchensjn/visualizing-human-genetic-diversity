import numpy as np
import gzip


global_mafs = []

with gzip.open("assets/Dryad/new_1kg_nyc_hg38_filt_total.biallelic_snps.pops.freq.total.txt.gz",'r') as f:
    header = f.readline()
    i = 0
    for line in f:
        if i < 100:
            global_mafs.append(line.split()[5])
        else:
            break
        i+= 1

global_mafs = np.array(global_mafs).astype(np.float64)
np.histogram(global_mafs)