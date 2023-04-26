# Description: Snakemake pipeline to predict run ENFORMER at inference
# Author: Temi
# Date: Wed Mar 29 2023

import pandas as pd
import os, glob, sys, re
from snakemake.io import glob_wildcards

sys.path.append('workflow/scripts')
sys.path.append('modules')
import module

print(os.getcwd())

cfile = "config/pipeline.yaml"
configfile: cfile

metadata_dt = pd.read_csv(config['metadata'])
valid_dt = pd.read_csv(os.path.join(os.path.dirname(config['metadata']), 'valid_TFs.csv'))

print(valid_dt)
rscript = config['rscript']

config_file_path = os.path.abspath(cfile)
# directories
DATA_DIR = 'data'
PREDICTION_PARAMS_DIR = os.path.join(DATA_DIR, 'prediction_parameters')

# every other thing is unnecessary; so, this globs the enformer config files
# and uses them to generate the aggregation config files
# agg_wildcards = glob_wildcards(os.path.join(PREDICTION_PARAMS_DIR, f"enformer_parameters_{config['dataset']}_{{tfs}}_{{tissues}}.json"))

# filter for only those available in the metadata
tfs = valid_dt['transcription_factor'].tolist()
tissues = [m.replace(' ', '-') for m in valid_dt['tissue'].tolist()]

# print(tfs)
# print(tissues)

rule all:
    input:
        expand(os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf = tfs, tissue = tissues)

rule predict_with_enformer:
    input:
        pfile = os.path.join(DATA_DIR, 'predictor_files', '{tf}_{tissue}_predictors.txt')
    output:
        os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json')
    params:
        jobname = '{tf}_{tissue}',
        cf = os.path.join(PREDICTION_PARAMS_DIR, f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json')
    message: 
        "working on {params.jobname}"
    # resources:
    #     gpus = config['enformer']['prediction_directives']['parsl_parameters']['num_of_full_nodes'],
    #     partition="beagle3"
    shell:
        """
            python3 predict/workflow/scripts/enformer_predict.py --param_config {params.cf}
        """
    
