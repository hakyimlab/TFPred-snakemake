# Description: Given a TF and tissue (or context) , this pipeline trains logistic elastic net models of that TF binding activity in that tissue
# Author: Temi
# Date: Wed Mar 29 2023
# Usage: --> see README

if config["usage"]['offline'] == True:
    if config["usage"]['software'] == 'singularity':
        singularity: config["usage"]['location']
    elif config["usage"]['software'] == 'conda':
        conda: config["usage"]['location']
elif config["usage"]['offline'] == False:
    sys.exit(1)

#"/beagle3/haky/users/temi/projects/build_singularity/def_files/tfxcan_sha256.c9195e1a359cee360f3743d68326208238257e584cbbee36ac6240acba8c0613.sif"
#"library://temi/collection/tfxcan:v1.0" #"/beagle3/haky/users/temi/projects/build_singularity/rb-66e869a3e9e27d2c5140ef22_latest.sif"

import pandas as pd
import os, glob, sys, re, yaml, subprocess
from snakemake.io import glob_wildcards
import numpy as np
import itertools
# from collections import Iterable

sys.path.append('workflow/src')
sys.path.append('workflow/modules')

import helpers

print_progress = False


runname = config['dataset']
rundate = config['date']
run = f'{runname}_{rundate}'

# directories ======

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
STATISTICS_DIR = os.path.join(RUN_DIR, 'statistics') #'output/statistics'


# ====== software ======
RSCRIPT = config['rscript'] if 'rscript' in config.keys() else 'Rscript'
# print(os.system('echo $PATH'))
HOMERDIR = config['homer']['dir'] if 'homer' in config.keys() else '/software/conda_envs/TFXcan-snakemake/share/homer' # if using singularity
# print(os.path.isdir(HOMERDIR))
HOMERSCAN = os.path.join(HOMERDIR, 'bin', 'scanMotifGenomeWide.pl')
HOMERGENOME = os.path.join(HOMERDIR, 'data', 'genomes', 'hg38')
HOMERMOTIFSDATABASE = os.path.join(HOMERDIR, 'motifs') if not 'homer' in config.keys() else config['homer']['motifs_database']
PYTHON3 = '/software/conda_envs/TFXcan-snakemake/bin/python3' if config['usage']['software'] == 'singularity' else 'python3'

#& 'motifs_database' not in config['homer'].keys() 


# try:
#     HOMERDIR = config['homer']['dir'] if 'homer' in config.keys() else '/software/conda_envs/TFXcan-snakemake/share/homer' # if using singularity
#     print(os.path.isdir(HOMERDIR))
#     HOMERSCAN = os.path.join(HOMERDIR, 'bin', 'scanMotifGenomeWide.pl')
#     HOMERGENOME = os.path.join(HOMERDIR, 'data', 'genomes', 'hg38')
#     HOMERMOTIFSDATABASE = os.path.join(HOMERDIR, 'motifs') if 'motifs_database' not in config['homer'].keys() else config['homer']['motifs_database']
# except Exception as e:
#     print("ERROR - [FATAL] Please specify the path to the homer data software. Exiting...")
#     sys.exit(1)

# print(HOMERSCAN)
# print(HOMERGENOME)
# print(HOMERMOTIFSDATABASE)

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
motif_inputs = helpers.createMotifInputs(homer_motifs_dict, HOMERMOTIFSDATABASE)
motif_outputs = helpers.createMotifOutputs(homer_motifs_dict, os.path.join(HOMERFILES_DIR, 'scannedMotifs'))
# check if the motif files are available so you don't have to re-run them with Homer, which takes some time
#motifs_available = [os.path.exists(m) for m in motif_outputs]
# which ones are not 
motifs_unavailable = [m for m in motif_outputs if not os.path.exists(m)]

# print(motifs_available)

# if not all([os.path.exists(m) for m in motif_files]):
#     print(f"ERROR - [FATAL] Please check the motif files. Exiting...")
#     sys.exit(1)

onstart:
    print(f"INFO - Found Rscript at: {RSCRIPT}")
    print(f"INFO - Found Homer at: {HOMERDIR}")
    print(f"INFO - Found python3 at: {PYTHON3}")

    if len(TF_list) < 5:
        print(f'INFO - Assays to train Enpact models for: {TF_list}')
    else:
        print(f"INFO - Verified {len(TF_list)} assays to train Enpact models for.")

    if not motifs_unavailable:
        print("INFO - All scanned motif files are available. Skipping Homer motif scanning...")
    else:
        print("INFO - Some scanned motif files are not available. Running Homer motif scanning...")

    if config['run_enformer'] == True:
        print(f'INFO - Running ENFORMER is set to True. Running the pipeline with ENFORMER...')
        include: 'workflow/rules/train_enpact.slow.smk'
    elif config['run_enformer'] == False:
        print(f'INFO - Running ENFORMER is set to False. Running the pipeline without ENFORMER...')
        include: 'workflow/rules/train_enpact.fast.smk'
    else:
        print(f"ERROR - [FATAL] Please specify whether to run ENFORMER or not. Exiting...")
        sys.exit(1)

rule all:
    input:
        motifs_unavailable,
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt'), tf = set(TF_list)),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'), zip, tf = TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'), zip, tf = TF_list, tissue=tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'{runname}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'train_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'test_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.logistic.rds'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.linear.rds'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.train_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.train_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.test_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.test_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        os.path.join(STATISTICS_DIR, f'{run}.compiled_stats.txt'),
        # os.path.join('reports', f'{run}.report.html')