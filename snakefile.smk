# Description: Given a TF and tissue (or context) , this pipeline trains logistic elastic net models of that TF binding activity in that tissue
# Author: Temi
# Date: Wed Mar 29 2023
# Usage: --> see README

import pandas as pd
import os, glob, sys, re, yaml, subprocess
from snakemake.io import glob_wildcards
import numpy as np
import itertools
from collections import Iterable

sys.path.append('workflow/src')
sys.path.append('workflow/modules')

import helpers

print_progress = False

runname = config['dataset']
rundate = config['date']
run = f'{runname}_{rundate}'

# directories

#DATA_DIR = os.path.join('data') 
DATA_DIR = 'data' # here I want to have some common files (like the motif files; this should not need to run everytime if already available)
HOMERFILES_DIR = os.path.join(DATA_DIR, 'homer_instances')

METADATA_DIR = 'metadata'

RUN_DIR = os.path.join(DATA_DIR, f"{runname}_{rundate}") 
BEDLINKS_DIR = os.path.join(RUN_DIR, 'bed_links')
SORTEDBEDS_DIR = os.path.join(RUN_DIR, 'sortedbeds')
PREDICTORS_DIR = os.path.join(RUN_DIR, 'predictor_files')
PREDICTION_PARAMS_DIR = os.path.join(RUN_DIR, 'prediction_parameters')
PREDICTIONS_DIR = os.path.join(config['scratch_dir'], 'predictions_folder') if os.path.exists(config['scratch_dir']) else os.path.join(RUN_DIR, 'predictions_folder')
AGGREGATION_DIR = os.path.join(RUN_DIR, 'aggregation_folder')
MODELS_DIR = os.path.join(RUN_DIR, 'models') #'output/models'
MODELS_EVAL_DIR = os.path.join(RUN_DIR, 'evaluation') #'output/models_eval'

# prepare input ======
metadata_dt = pd.read_table(config['models_metadata'], dtype={'assay': 'string', 'context': 'string'})
metadata_dt = metadata_dt.fillna('none')
##print(metadata_dt)

# reading YAML file
with open(config['models_config']) as stream:
    try:
        model_config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# verify details
tp = tuple(metadata_dt.itertuples(index=False, name = None))
vfy = helpers.verify_model_details(tp, model_config, print_info = True)

if not vfy:
    print("ERROR - [FATAL] Please check the metadata file and the model configuration file for consistency. The assay column should be present in the models configuration yaml file. Exiting...")
    sys.exit(1)

metadata_dt = metadata_dt[metadata_dt.assay.isin(vfy)]
model_config = {k: model_config[k] for k in vfy}

details = []
for row in metadata_dt.itertuples():
    #print(row)
    r = [row.assay, row.context]
    details.append(r)

unique_TFs = list(set(metadata_dt.assay))
homer_motifs_dict = helpers.collectMotifFiles(unique_TFs, model_config)
TF_list, tissue_list = [d[0] for d in details], [d[1].replace(' ', '-') for d in details]
motif_files = list(homer_motifs_dict.values())
motif_inputs = helpers.createMotifInputs(homer_motifs_dict, config['homer']['motifs_database'])
motif_outputs = helpers.createMotifOutputs(homer_motifs_dict, HOMERFILES_DIR)

if config['run_enformer'] == True:
    include: workflow/rules/train_enpact.slow.smk
elif config['run_enformer'] == False:
    include: workflow/rules/train_enpact.fast.smk
else:
    print("ERROR - [FATAL] Please specify whether to run ENFORMER or not. Exiting...")
    sys.exit(1)