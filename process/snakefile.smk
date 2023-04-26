# Descriptin: snakefile to process ENFORMER predictions and train elastic net models
# Author: Temi
# Date: Mon Apr 17 2023

import pandas as pd
import os, glob, sys, re, json
from snakemake.io import glob_wildcards

sys.path.append('process/workflow/scripts')
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

# filter for only those available in the metadata
tfs = valid_dt['transcription_factor'].tolist()
tissues = [m.replace(' ', '-') for m in valid_dt['tissue'].tolist()]

# directories
DATA_DIR = 'data'
PREDICTORS_DIR = os.path.join(DATA_DIR, 'predictor_files')
MODELS_DIR = 'output/models'
MODELS_EVAL_DIR = 'output/models_eval'
PREDICTION_PARAMS_DIR = os.path.join(DATA_DIR, 'prediction_parameters')

dts = zip(tfs, tissues)
motif_grouping_dict = module.group_tf_motif_files(dts)
#tissue_grouping_dict = module.group_tf_tissues(TF_tissue_list)

valid_chromosomes = [f'chr{i}' for i in range(1,23)] + ['chrX']

prediction_grouping_paths = {}
for d in zip(tfs, tissues):
    tf, tissue = d[0], d[1]
    p = module.return_prediction_folder(os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{tf}_{tissue}.json'))
    prediction_grouping_paths[f'{tf}_{tissue}'] = p

prediction_run_date = {}
for d in zip(tfs, tissues):
    tf, tissue = d[0], d[1]
    p = module.return_prediction_date(os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{tf}_{tissue}.json'))
    prediction_run_date[f'{tf}_{tissue}'] = p

aggregated_folders = prediction_grouping_paths.values()
tf_tissue_pairs = prediction_grouping_paths.keys() #['_'.join(a) for a in TF_tissue_list]



rule all:
    input:
        expand(os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{{tf_tissue}}.json'), tf_tissue = tf_tissue_pairs),
        expand(os.path.join('{aggfolder}', 'aggregated_predictions', f'{config["dataset"]}_aggByMeanCenter_{{tf_tissue}}.csv'), zip, aggfolder = prediction_grouping_paths.values(), tf_tissue = tf_tissue_pairs),
        expand(os.path.join("{aggfolder}", "aggregated_predictions", f"train_{config['dataset']}_aggByMeanCenter_{{tf_tissue}}.prepared.csv"), zip, aggfolder = prediction_grouping_paths.values(), tf_tissue = tf_tissue_pairs),
        expand(os.path.join("{aggfolder}", "aggregated_predictions", f"test_{config['dataset']}_aggByMeanCenter_{{tf_tissue}}.prepared.csv"), zip, aggfolder = prediction_grouping_paths.values(), tf_tissue = tf_tissue_pairs),
        expand(os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", 'aggByMeanCenter_{tf_tissue}.linear.rds'), zip, date = prediction_run_date.values(), tf_tissue = tf_tissue_pairs),
        expand(os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", 'aggByMeanCenter_{tf_tissue}.logistic.rds'), zip, date = prediction_run_date.values(), tf_tissue = tf_tissue_pairs),
        expand(os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", f'aggByMeanCenter_test_evaluation.linear.rds'), zip, date = prediction_run_date.values(), tf_tissue = tf_tissue_pairs),
        expand(os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", f'aggByMeanCenter_test_evaluation.logistic.rds'), zip, date = prediction_run_date.values(), tf_tissue = tf_tissue_pairs),
        expand(os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", f'aggByMeanCenter_train_evaluation.linear.rds'), zip, date = prediction_run_date.values(), tf_tissue = tf_tissue_pairs),
        expand(os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", f'aggByMeanCenter_train_evaluation.logistic.rds'), zip, date = prediction_run_date.values(), tf_tissue = tf_tissue_pairs)


rule aggregate_predictions:
    input:
        os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{{tf_tissue}}.json')
    output:
        os.path.join('{aggfolder}', 'aggregated_predictions', f"{config['dataset']}_aggByMeanCenter_{{tf_tissue}}.csv")
    message: 
        "working on {wildcards}"
    params:
        jobname = '{tf_tissue}'
    resources:
        mem_mb= 100000
    shell:
        """
            python3 process/workflow/scripts/aggregate.py --metadata_file {input} --agg_types "aggByMeanCenter"
        """

rule prepare_training_data:
    input:
        p1 = rules.aggregate_predictions.output,
        p2 = os.path.join(PREDICTORS_DIR, f'{{tf_tissue}}_ground_truth.txt')
    output:
        p1=os.path.join("{aggfolder}", "aggregated_predictions", f"train_{config['dataset']}_aggByMeanCenter_{{tf_tissue}}.prepared.csv"),
        p2=os.path.join("{aggfolder}", "aggregated_predictions", f"test_{config['dataset']}_aggByMeanCenter_{{tf_tissue}}.prepared.csv")
    message: 
        "preparing {wildcards.tf_tissue} training and test data"
    params:
        jobname = '{tf_tissue}' 
    resources:
        mem_mb= 100000
    shell:
        """
            Rscript process/workflow/scripts/prepare_training_data.R {input.p1} {input.p2} "aggByMeanCenter" {output.p1} {output.p2}
        """

rule train_TFPred_weights:
    input:
        lambda wildcards: os.path.join(prediction_grouping_paths[wildcards.tf_tissue], 'aggregated_predictions', f'train_{config["dataset"]}_aggByMeanCenter_{wildcards.tf_tissue}.prepared.csv')
    output:
        os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", f'aggByMeanCenter_{{tf_tissue}}.logistic.rds'),
        os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", f'aggByMeanCenter_{{tf_tissue}}.linear.rds')
    message:
        "training on {wildcards.tf_tissue} training data"
    params:
        jobname = '{tf_tissue}',
        basename = os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", f'aggByMeanCenter_{{tf_tissue}}')
    resources:
        mem_mb= 100000
    shell:
        """
            Rscript process/workflow/scripts/train_enet.R {input} {params.basename}
        """

rule evaluate_TFPred:
    input:
        test_data = lambda wildcards: os.path.join(prediction_grouping_paths[wildcards.tf_tissue], 'aggregated_predictions', f'test_{config["dataset"]}_aggByMeanCenter_{wildcards.tf_tissue}.prepared.csv'),
        train_data = lambda wildcards: os.path.join(prediction_grouping_paths[wildcards.tf_tissue], 'aggregated_predictions', f'train_{config["dataset"]}_aggByMeanCenter_{wildcards.tf_tissue}.prepared.csv')
    output:        
        os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", 'aggByMeanCenter_test_evaluation.logistic.rds'),
        os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", 'aggByMeanCenter_test_evaluation.linear.rds'),
        os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", 'aggByMeanCenter_train_evaluation.logistic.rds'),
        os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", 'aggByMeanCenter_train_evaluation.linear.rds')
    message:
        "evaluating on {wildcards.tf_tissue} train and test sets"
    params:
        jobname = '{tf_tissue}',
        model_basename = os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}", 'aggByMeanCenter_{tf_tissue}'),
        #data_basename = lambda wildcards: os.path.join(prediction_grouping_paths[wildcards.tf_tissue], 'aggregated_predictions', f'test_{config["dataset"]}_aggByMeanCenter_{wildcards.tf_tissue}'),
        output_dir = os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf_tissue}}_{{date}}")
    resources:
        mem_mb= 100000
    shell:
        """
            Rscript process/workflow/scripts/evaluate_enet.R {params.model_basename} {input.test_data} {params.output_dir} 'test'
            Rscript process/workflow/scripts/evaluate_enet.R {params.model_basename} {input.train_data} {params.output_dir} 'train'
        """
    
    
