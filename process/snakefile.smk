# Usage: snakefile to process ENFORMER predictions
# Author: Temi
# Date: Mon Apr 17 2023

import pandas as pd
import os, glob, sys, re, json
from snakemake.io import glob_wildcards

sys.path.append('process/workflow/scripts')
print(os.getcwd())



cfile = "config/pipeline.yaml"
configfile: cfile

metadata_dt = pd.read_csv(config['metadata'])
rscript = config['rscript']

if not os.path.isdir(config['enformer']['prediction_directives']['metadata_dir']):
    os.makedirs(config['enformer']['prediction_directives']['metadata_dir'])

config_file_path = os.path.abspath(cfile)

# print(config_file_path)
# print(metadata_dt)

# directories
DATA_DIR = 'data'
BEDFILES_DIR = os.path.join(DATA_DIR, 'bed_files')
HOMERFILES_DIR = os.path.join(DATA_DIR, 'homer_files')
PREDICTORS_DIR = os.path.join(DATA_DIR, 'predictor_files')
METADATA_DIR = 'metadata'
MODELS_DIR = 'models'

details = []
for row in metadata_dt.itertuples():
    #print(row)
    r = [row.transcription_factor, row.tissue]
    details.append(r)
print(details)

homer_wildcards = glob_wildcards(os.path.join(HOMERFILES_DIR, '{tfs}/{motif_files}.motif'))

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

# need to read in the predictions folder
def return_prediction_folder(json_file_path):
    with open(json_file_path) as f:
        data = json.load(f)
        p = data['predictions_folder']
    return(p)

def return_prediction_date(json_file_path):
    with open(json_file_path) as f:
        data = json.load(f)
        p = data['run_date']
    return(p)


prediction_grouping_paths = {}
for d in details:
    tf = d[0]
    tissue = d[1]
    p = return_prediction_folder(os.path.join(METADATA_DIR, 'enformer_config', f'aggregation_config_{config["dataset"]}_{tf}_{tissue}.json'))
    print(p)
    prediction_grouping_paths[f'{tf}_{tissue}'] = p
print(prediction_grouping_paths)

prediction_run_date = {}
for d in details:
    tf = d[0]
    tissue = d[1]
    p = return_prediction_date(os.path.join(METADATA_DIR, 'enformer_config', f'aggregation_config_{config["dataset"]}_{tf}_{tissue}.json'))
    print(p)
    prediction_run_date[f'{tf}_{tissue}'] = p
print(prediction_run_date)

aggregated_folders = prediction_grouping_paths.values()
print(aggregated_folders)

print(prediction_run_date)

print(TF_list)
print(tissue_list)

rule all:
    input:
        expand(os.path.join(METADATA_DIR, 'enformer_config', f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join('{aggfolder}', 'aggregated_predictions', f'{config["dataset"]}_aggByMeanCenter_{{tf_tissue}}.csv.gz'), zip, aggfolder = prediction_grouping_paths.values(), tf_tissue = prediction_grouping_paths.keys()),
        expand(os.path.join("{aggfolder}", "aggregated_predictions", f"train_{config['dataset']}_aggByMeanCenter_{{tf}}_{{tissue}}.prepared.csv.gz"), zip, aggfolder = prediction_grouping_paths.values(),  tf = TF_list, tissue = tissue_list),
        expand(os.path.join("{aggfolder}", "aggregated_predictions", f"test_{config['dataset']}_aggByMeanCenter_{{tf}}_{{tissue}}.prepared.csv.gz"), zip, aggfolder = prediction_grouping_paths.values(), tf = TF_list, tissue = tissue_list),
        expand(os.path.join('{aggfolder}', 'aggregated_predictions', f'{config["dataset"]}_aggByMeanCenter_{{tf_tissue}}.csv.gz'), zip, aggfolder = prediction_grouping_paths.values(), tf_tissue = prediction_grouping_paths.keys()),
        expand(os.path.join(MODELS_DIR, '{tf_tissue}_{date}', f'{config["dataset"]}_aggByMeanCenter_{{tf_tissue}}.linear.rds'), zip, date = prediction_run_date.values(), tf_tissue = prediction_run_date.keys()),
        expand(os.path.join(MODELS_DIR, '{tf_tissue}_{date}', f'{config["dataset"]}_aggByMeanCenter_{{tf_tissue}}.logistic.rds'), zip, date = prediction_run_date.values(), tf_tissue = prediction_run_date.keys())


rule aggregate_predictions:
    input:
        os.path.join(METADATA_DIR, 'enformer_config', f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json')
    output:
        os.path.join("{aggfolder}", "aggregated_predictions", f"{config['dataset']}_aggByMeanCenter_{{tf}}_{{tissue}}.csv.gz")
    message: 
        "working on {wildcards}"
    params:
        jobname = '{tf}_{tissue}'
    shell:
        """
            python3 process/workflow/scripts/aggregate.py --metadata_file {input} --agg_types "aggByMeanCenter"
        """

rule prepare_training_data:
    input:
        p1 = rules.aggregate_predictions.output,
        p2 = os.path.join(DATA_DIR, 'predictor_files', f'{{tf}}_{{tissue}}_ground_truth.txt')
    output:
        p1=os.path.join("{aggfolder}", "aggregated_predictions", f"train_{config['dataset']}_aggByMeanCenter_{{tf}}_{{tissue}}.prepared.csv.gz"),
        p2=os.path.join("{aggfolder}", "aggregated_predictions", f"test_{config['dataset']}_aggByMeanCenter_{{tf}}_{{tissue}}.prepared.csv.gz")
    message: 
        "working on {wildcards}"
    params:
        jobname = '{tf}_{tissue}',
        p1 = os.path.join("{aggfolder}", "aggregated_predictions", f"train_{config['dataset']}_aggByMeanCenter_{{tf}}_{{tissue}}.prepared.csv"),
        p2 = os.path.join("{aggfolder}", "aggregated_predictions", f"test_{config['dataset']}_aggByMeanCenter_{{tf}}_{{tissue}}.prepared.csv") 
    shell:
        """
            Rscript process/workflow/scripts/prepare_training_data.R {input.p1} {input.p2} "aggByMeanCenter" {params.p1} {params.p2}
            pigz -kf {params.p1} && pigz -kf {params.p2}
        """

rule train_TFPred:
    input:
        rules.prepare_training_data.output.p1
    output:
        os.path.join(MODELS_DIR, '{tf_tissue}_{date}', f'{config["dataset"]}_aggByMeanCenter_{{tf_tissue}}.logistic.rds'),
        os.path.join(MODELS_DIR, '{tf_tissue}_{date}', f'{config["dataset"]}_aggByMeanCenter_{{tf_tissue}}.linear.rds')
    message:
        "training on {wildcards}"
    params:
        jobname = '{tf_tissue}',
        basename = os.path.join(MODELS_DIR, '{tf_tissue}_{date}', f'{config["dataset"]}_aggByMeanCenter_{{tf_tissue}}')
    shell:
        """
            Rscript process/workflow/scripts/train_TFPred.R {input} {params.basename}
        """