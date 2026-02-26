import os, h5py, re
import numpy as np
import math
import warnings
import pandas as pd
import argparse, os, sys, multiprocessing, itertools
import pandas as pd
import numpy as np


def calculate_genomic_coordinates(locus, num_bins = 896):
    locus = locus.split('_')
    start = int(locus[1])
    end = int(locus[2])
    midn = math.ceil((start + end) / 2)
    # from position choose center bin
    center_ind = midn - 1
    center_bin = center_ind // 128
    half_bins = num_bins // 2
    start_bin = center_bin - half_bins
    end_bin = center_bin + half_bins
    if num_bins % 2 != 0: # if num_bins is odd
        end_bin += 1
    return([start_bin, end_bin])

def calculate_genomic_coordinates2(locus, bin_size = 128):
    locus = locus.split('_')
    start = int(locus[1])
    end = int(locus[2])
    start_bin = start // bin_size
    end_bin = end // bin_size
    return([start_bin, end_bin])

def extract_enformer_predictions(enfref_dir, chr_num, locs):
    with h5py.File(os.path.join(enfref_dir, f"chr{chr_num}_cat.h5"), "r") as f:
        epigen = f[f'chr{chr_num}'][locs[0]:locs[1], :] 
    return epigen


def slice_bins(locus, bin_size=128, nbins=896):
    locus = locus.split('_')
    start = int(locus[1])
    end = int(locus[2])
    midn = math.ceil((start + end) / 2)
    nstart = midn - ((bin_size*nbins) / 2) # (128*896) / 2
    sstart = (start - nstart)
    send = sstart + ((end - start) - 1)
    bins = list(range(0, nbins*bin_size, bin_size))
    out = np.digitize([sstart, send], bins=bins).tolist()
    if((end - start) <= 128):
        pass
    elif(((end - start) > 128) or (end - start) % bin_size > 0):
        out[1] = out[1] + 1
    # if [448], make it [448, 449]
    # if [448, 448] make it [448, 449]
    if len(out) == 1:
        out.append(out[0] + 1)
    elif len(out) == 2:
        if out[0] == out[1]:
            out[1] = out[1] + 1
    #print(f"INFO - Bins to take: {out}")
    return(out)
    

def aggregate_and_collect_epigenome(locus, reference_epigenome_dir, pad_bins = 1):
    output = dict()
    # get the chrom
    chrom = locus.split('_')[0]
    chrom = int(re.sub('chr', '', chrom))
    # get bins coords
    bins_coords = calculate_genomic_coordinates(locus)
    enf_pred = extract_enformer_predictions(enfref_dir=reference_epigenome_dir, chr_num=chrom, locs=bins_coords)
    if pad_bins == -1:
        # return all the bins
        output['locus'] = locus
        output['values'] = enf_pred
    else:
        bins_to_take = slice_bins(locus)
        bins_to_aggregate = [bins_to_take[0] - pad_bins, bins_to_take[1] + pad_bins] # add one more bin upstream and downstrea
        print(f"INFO - After padding bins with {pad_bins}: {bins_to_aggregate}")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            try:
                output['locus'] = locus
                output['values'] = enf_pred[bins_to_aggregate[0]:bins_to_aggregate[1], :].mean(axis = 0)
            except RuntimeWarning as rw: # this probably won't work since I am ignoring the warning anyway
                print(f'ERROR - {locus} has encountered an error')
                output['locus'] = None
                output['values'] = None
    return(output)



# def save_output(dict_object, output_file_basename):
#     import pickle 
#     with open(f'{output_file_basename}.pkl', 'wb') as ff:
#         pickle.dump(dict_object, ff)

# def read_output(fileX):
#     import pickle
#     with open(fileX, 'rb') as f:
#         output_dictionary = pickle.load(f)
#     return(output_dictionary)


# def extract_enformer_predictions(enfref_dir, chr_num, locs):
#     with h5py.File(os.path.join(enfref_dir, f"chr{chr_num}_cat.h5"), "r") as f:
#         epigen = f[f'chr{chr_num}'][locs[0]:locs[1], :] 
#     return epigen

# parameters = dict()
# parameters['input_locus'] = 'chr10_102949239_104380410'
# parameters['reference_epigenome_dir'] = '/beagle3/haky/data/enformer-reference-epigenome'
# bin_indices_to_extract = calculate_genomic_coordinates(parameters['input_locus'], num_bins = 128)
# extracted_epigenome = extract_enformer_predictions(enfref_dir=parameters['reference_epigenome_dir'], chr_num=10, locs=bin_indices_to_extract)
# extracted_epigenome = pd.DataFrame(extracted_epigenome)
# extracted_epigenome.columns = [f'f_{i}' for i in range(1, extracted_epigenome.shape[1]+1)]
# extracted_epigenome.shape
# #np.savetxt('/beagle3/haky/users/temi/misc/test_extracted_epigenome.csv', extracted_epigenome, delimiter='\t')
# extracted_epigenome.to_csv('/beagle3/haky/users/temi/misc/test_extracted_epigenome.csv.gz', index=False, sep='\t', compression = 'gzip')

# 815471 - 804290

# 809945 - 809817


# # I want to select 896 bins from the midpoint of the locus
# # but you can select more by increasing num_bins
# calculate_genomic_coordinates(parameters['input_locus'], num_bins = 896)

# # for a 2MB window, this would mean 2MB/128 = 15625 bins
# # so you would do
# calculate_genomic_coordinates(parameters['input_locus'], num_bins = 11181.0234375)

# #
# calculate_genomic_coordinates2(parameters['input_locus'], bin_size = 128)

# 1431171 / 128
# start = 102949239
# end = 104380410
# end - start

# midn = math.ceil((start + end) / 2)
# # from position choose center bin
# center_ind = midn - 1
# center_bin = center_ind // 128
# half_bins = num_bins // 2
# start_bin = center_bin - half_bins
# end_bin = center_bin + half_bins
# if num_bins % 2 != 0: # if num_bins is odd
#     end_bin += 1

# # read in an hdf5 file
# with h5py.File('/beagle3/haky/data/enformer-reference-epigenome/chr22_cat.h5', "r") as bigfile:
#     shh = bigfile[f'chr22'].shape
