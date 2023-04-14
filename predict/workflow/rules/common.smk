# Usage: Common snakefile
# Author: Temi
# Date: Wed Mar 29 2023

import pandas as pd
import os, glob, sys, re
from snakemake.io import glob_wildcards

sys.path.append('workflow/scripts')

# try:
#     import module
# except ModuleNotFoundError as merr:
#     raise Exception("ERROR - modules not found")

cfile = "../config/pipeline.yaml"
configfile: cfile

metadata_dt = pd.read_csv(config['metadata'])
rscript = config['rscript']

if not os.path.isdir(config['enformer']['prediction_directives']['metadata_dir']):
    os.makedirs(config['enformer']['prediction_directives']['metadata_dir'])

config_file_path = os.path.abspath("config/pipeline.yaml")

# print(config_file_path)
# print(metadata_dt)

# directories
BEDFILES_DIR = 'bed_files'
HOMERFILES_DIR = 'homer_files'

details = []
for row in metadata_dt.itertuples():
    #print(row)
    r = [row.transcription_factor, row.tissue]
    details.append(r)
print(details)


BED_FOLDERS = ['_'.join(a) for a in details] # e.g. 'FOXA1_Breast/'

def return_cistrome_ids(TF, tissue, cistrome_db):
    db = pd.read_table(cistrome_db)
    ids = db.loc[(db['Factor'] == TF) & (db['Tissue_type'] == tissue)].loc[:,'DCid'].tolist()
    return(ids)

def return_cistrome_bedfiles(TF, tissue, cistrome_db, cistrome_db_folder, suffix=None, prefix=None):

    bedfiles_ids = return_cistrome_ids(TF=TF, tissue=tissue, cistrome_db = cistrome_db)

    patterns = [f'{os.path.join(cistrome_db_folder, f"{id}_*.bed")}' for id in bedfiles_ids]
    #patterns = [f"{id}_*.bed" for id in bedfiles_ids]
    patterns = [glob.glob(p) for p in patterns]
    patterns = [os.path.basename(p) for p in patterns for p in p]
    return({f'{TF}_{tissue}': patterns})

def link_cistrome_bedfiles(TF, tissue):

    bedfiles_ids = return_cistrome_bedfiles(TF=TF, tissue=tissue, cistrome_db_folder = config['cistrome_data_dir'], cistrome_db = config['TF_table'])


    for key, values in bedfiles_ids.items():
        dst_folder = os.path.join('data', BEDFILES_DIR, key)
        if not os.path.isdir(dst_folder):
            os.makedirs(dst_folder)
        for value in values:
            dst_file = os.path.join(dst_folder, value)
            if not os.path.isfile(dst_file):
                os.symlink(os.path.join(config['cistrome_data_dir'], value), dst_file)
    return(0)

# homer motifs
def link_homer_motifs(TF, tissue):
    homer_data = os.path.join(config['homer']['dir'], 'data/knownTFs/motifs')
    #print(homer_data)
    dst_folder = os.path.join('data', HOMERFILES_DIR, TF)
    if not os.path.isdir(dst_folder):
        os.makedirs(dst_folder)
    #tf = TF.lower()
    if TF.lower() == 'ar':
        use_tf = 'ar-half'
    else:
        use_tf = TF.lower()
    search_motif = os.path.join(homer_data, f'{use_tf}*.motif')
    motif_files = glob.glob(search_motif)
    if len(motif_files) >= 1:
        #print(f"WARNING - more than one motif found for {TF}")
        for f in motif_files:
            dst_file = os.path.join(dst_folder, os.path.basename(f))
            if not os.path.isfile(dst_file):
                os.symlink(f, dst_file)
    return(0)    

# link all the bed files
files_linked = [link_cistrome_bedfiles(TF=d[0], tissue=d[1]) for d in details]
files_linked_2 = [link_homer_motifs(TF=d[0], tissue=d[1]) for d in details]
homer_wildcards = glob_wildcards(os.path.join('data', HOMERFILES_DIR, '{tfs}/{motif_files}.motif'))

def group_tf_motif_files(dts):
    grouping_dict = {}
    for dt in dts:
        if dt[0] not in grouping_dict.keys():
            grouping_dict[dt[0]] = [dt[1]]
        else:
            grouping_dict[dt[0]].append(dt[1])
    return(grouping_dict)

def group_tf_tissues(dts):
    grouping_dict = {}
    for dt in dts:
        if dt[0] not in grouping_dict.keys():
            grouping_dict[dt[0]] = [dt[1]]
        else:
            grouping_dict[dt[0]].append(dt[1])
    return(grouping_dict)


dts = zip(homer_wildcards.tfs, homer_wildcards.motif_files)
motif_grouping_dict = group_tf_motif_files(dts)
tissue_grouping_dict = group_tf_tissues(details)
valid_chromosomes = [f'chr{i}' for i in range(1,23)] + ['chrX']

TF_list = [d[0] for d in details]
tissue_list = [d[1] for d in details]

# print(motif_grouping_dict)
# print(tissue_grouping_dict)


rule all:
    input:
        # expand(os.path.join('../prepare/data', 'predictor_files', '{tf}_{tissue}_predictor_regions.txt'), zip, tf = TF_list, tissue = tissue_list),
        # expand(os.path.join('../prepare/data', 'predictor_files', '{tf}_{tissue}_predictors.txt'), zip, tf = TF_list, tissue = tissue_list),
        # expand(os.path.join('../prepare/data', 'predictor_files', '{tf}_{tissue}_ground_truth.txt'), zip, tf = TF_list, tissue = tissue_list),
        #expand(os.path.join('metadata', 'enformer_config', f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join('metadata', 'enformer_config', f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf = TF_list, tissue = tissue_list)
        #/aggregation_config_{prediction_data_name}_{prediction_id}.json'

rule predict_with_enformer:
    input:
        pfile = os.path.join('../prepare/data', 'predictor_files', '{tf}_{tissue}_predictors.txt')
    output:
        os.path.join('metadata', 'enformer_config', f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json')
    params:
        jobname = '{tf}_{tissue}',
        cf = os.path.join('../metadata', 'enformer_config', f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json')
    message: 
        "working on {params.jobname}"
    # resources:
    #     gpus = config['enformer']['prediction_directives']['parsl_parameters']['num_of_full_nodes'],
    #     partition="beagle3"
    shell:
        """
            echo `which python3`
            # source ~/.bashrc
            # conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
            # export LD_LIBRARY_PATH=/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib:$LD_LIBRARY_PATH
            # export PATH=$PATH:/project2/haky/temi/software/homer/bin

            python3 workflow/scripts/parallel_enformer/enformer_predict.py --param_config {params.cf}
        """
    
