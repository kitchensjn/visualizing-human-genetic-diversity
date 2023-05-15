import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pkg_resources
from geovar import *
from tqdm import tqdm
import gzip


def geovar_codes_streaming_fixed(geovar_obj, freq_mat_file):
    """Version of GeoVar code generation algorithm that streams through file to avoid memory overflow.
    Args:
        freq_mat_file (:obj:`string`): filepath to
        frequency table file (see example notebook for formatting).
    """
    assert geovar_obj.bins is not None
    geovar_codes = []
    # Setting up the testing bins
    test_bins = np.array([x[1] for x in geovar_obj.bins])
    with gzip.open(freq_mat_file,'r') as f:
        header = f.readline()
        # Take the population labels currently
        pops = np.array(header.split()[6:])
        geovar_obj.pops = pops
        for line in tqdm(f):
            # Split after the 6th column ...
            maf_vector = np.array(line.split()[6:]).astype(np.float64)
            cur_geovar = np.digitize(maf_vector, test_bins, right=True)
            cur_geovar_code = "".join([str(i) for i in cur_geovar])
            geovar_codes.append(cur_geovar_code)
    # Setting the variables here
    geovar_obj.geovar_codes = np.array(geovar_codes)
    geovar_obj.n_variants = geovar_obj.geovar_codes.size
    geovar_obj.n_populations = geovar_obj.pops.size


geovar = GeoVar(bins=[(0,0), (0,0.01), (0.01,0.05), (0.05,0.1), (0.1,1.0)])
geovar_codes_streaming_fixed(geovar_obj=geovar, freq_mat_file="assets/Dryad/new_1kg_nyc_hg38_filt_total.biallelic_snps.pops.freq.total.txt.gz")
counts_output = geovar.count_geovar_codes()
geovar_code_counts = pd.DataFrame({"codes":counts_output[0],"counts":counts_output[1]})
geovar_code_counts.to_csv("assets/new_1kg_nyc_hg38_filt_total.biallelic_snps.pops.ncat_000_001_005_010_100.geodist.total.txt", sep=' ', header=False, index=False)
