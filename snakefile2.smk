# Description: Given a TF and tissue (or context) , this pipeline trains logistic elastic net models of that TF binding activity in that tissue
# Author: Temi
# Date: Wed Mar 29 2023
# Usage: --> see README

import pandas as pd
import os, glob, sys, re, yaml
from snakemake.io import glob_wildcards
import numpy as np
import itertools
from collections import Iterable

sys.path.append('workflow/src')
sys.path.append('modules')

import module
print_progress = False

runname = config['dataset']
rundate = config['date']

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
MODELS_EVAL_DIR = os.path.join(RUN_DIR, 'models_eval') #'output/models_eval'

# prepare input ======
metadata_dt = pd.read_table(config['models_metadata'], dtype={'assay': 'string', 'context': 'string'})
metadata_dt = metadata_dt.fillna('none')
##print(metadata_dt)
details = []
for row in metadata_dt.itertuples():
    #print(row)
    r = [row.assay, row.context]
    details.append(r)

def get_mem_mb_allocations(wildcards, attempt):
    if attempt > 5:
        return(attempt * 10000)
    else:
        return(attempt * 10000)

def get_cluster_allocation(wildcards, attempt):
    if attempt > 5:
        return('bigmem')
    else:
        return('beagle3')

# # check if the TF is present in the homer and cistrome data_db columns
# data_db = pd.read_table('./info/data_db.txt')
# data_db = data_db.loc[(data_db['homer_db'] == 1) & (data_db['cistrome_db'] == 1)]
# valid_TFs = [d for d in details if d[0] in data_db.TF.tolist()] # filter the lists

# reading YAML file
with open(config['models_config']) as stream:
    try:
        model_config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# write out the list of unavailble TF-tissue pairs
# invalid_TFs = [d for d in details if d[0] not in data_db.TF.tolist()]
# if len(invalid_TFs) > 0:
#     print(f'Some invalid TF-tissue pairs were found. See {os.path.join(METADATA_DIR, "invalid_TFs.txt")} for details.')
#     pd.DataFrame(invalid_TFs).to_csv(os.path.join(METADATA_DIR, 'invalid_TFs.csv'), sep=',', index=False, header=['assay', 'context'])

# if len(valid_TFs) == 0:
#     print('No valid TF-tissue pairs were found. Please check the metadata file.')
#     sys.exit(1)
# else:
#     print(f'{len(valid_TFs)} valid TF-tissue pairs were found.')
#     pd.DataFrame(valid_TFs).to_csv(os.path.join(METADATA_DIR, 'valid_TFs.csv'), sep=',', index=False, header=['assay', 'context'])

#use this to filter for homer files
#print(model_config)
# print(details)


unique_TFs = list(set(metadata_dt.assay))
# print(unique_TFs)

def collectMotifFiles(TFs_list, model_config):
    motif_dict = {tf: model_config[tf]['motifFiles'] for tf in TFs_list}
    return(motif_dict)

def collectPeakFiles(model_details, model_config):
    peak_dict = dict()
    for det in model_details:
        tf = det[0]
        tissue = det[1]
        peak_file = model_config[tf]['peakFiles'][tissue]
        if ' ' in tissue:
            tissue = tissue.replace(' ', '-')
        mname = f'{tf}_{tissue}'
        peak_dict[mname] = peak_file
        return(peak_dict)




homer_motifs_dict = collectMotifFiles(unique_TFs, model_config)
print(f'===========================')
#print(homer_motifs_dict)
# print(homer_motifs.keys())
# print(homer_motifs.values())
# homer_motifs = list(zip(itertools.cycle(homer_motifs.keys()), list(homer_motifs.values())[0]))
# print(f'===========================')
# print(homer_motifs)

print(f'===========================')
# print(collectPeakFiles(details, model_config))

# mapped one-to-one ['AR', 'AR', 'FOXA'], ['Prostate', 'Bone-Marrow', 'Liver']
TF_list, tissue_list = [d[0] for d in details], [d[1].replace(' ', '-') for d in details]

# homer_motifs_wildcards = glob_wildcards(os.path.join(HOMERFILES_DIR, '{tfs}', '{motif_files}.motif'))

motif_files = list(homer_motifs_dict.values())
# motif_files = list(flatten(motif_files))



# print(zip(homer_motifs))

def createMotifInputs(motifsDict):
    res = []
    for key, values in motifsDict.items():
        if isinstance(values, list):
            res.extend(values)
        elif isinstance(values, str):
            res.append(values)
    return(res)

def createMotifOutputs(motifsDict):
    res = []
    for key, values in motifsDict.items():
        if isinstance(values, list):
            res.extend([os.path.join(HOMERFILES_DIR, f'{key}', f'scanMotifsGenomeWide.{v}.txt') for v in values])
        elif isinstance(values, str):
            res.append(os.path.join(HOMERFILES_DIR, f'{key}', f'scanMotifsGenomeWide.{values}.txt'))
    return(res)

motif_inputs = createMotifInputs(homer_motifs_dict)
motif_outputs = createMotifOutputs(homer_motifs_dict)
# print(motif_inputs)
# os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide_{motif_file}.txt')
rule all:
    input:
        motif_outputs,
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt'), tf = set(TF_list)),
#         # expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'), zip, tf=TF_list, tissue=tissue_list),
#         # expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'), zip, tf=TF_list, tissue=tissue_list),
#         # expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.info.txt.gz'), zip, tf=TF_list, tissue=tissue_list),
#         # expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.summary.txt'), zip, tf=TF_list, tissue=tissue_list),
#         # expand(os.path.join(PREDICTION_PARAMS_DIR, f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf=TF_list, tissue=tissue_list),
#         # expand(os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf = TF_list, tissue = tissue_list),
#         # expand(os.path.join(AGGREGATION_DIR, f'{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
#         # expand(os.path.join(AGGREGATION_DIR, f'train_{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
#         # expand(os.path.join(AGGREGATION_DIR, f'test_{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
#         # expand(os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.logistic.rds'), zip, tf = TF_list, tissue = tissue_list),
#         # expand(os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.rds'), zip, tf = TF_list, tissue = tissue_list),
#         # expand(os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.train_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
#         # expand(os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.logistic.train_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
#         # expand(os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.test_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
#         # expand(os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.logistic.test_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list)

# def gatherMotifFiles(wildcards):
#     mfs = glob_wildcards(os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide_{motif_file}.txt')).motif_file
#     return(mfs)

rule find_homer_motifs:
    # input:
    #     os.path.join(config["homer"]['motifs_database'], '{motif_file}')
    output: 
        os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide.{motif_file}.txt')
    params:
        homer_cmd = config['homer']['scanMotifsGenome'],
        genome = config['homer']['genome'],
        mfile = os.path.join(config['homer']['motifs_database'], '{motif_file}'),
        #ofile = lambda output: os.path.join(output, 'scanMotifsGenomeWide_{motif_file}.txt'),
        jobname = '{tf}'
    message: "working on {wildcards}" 
    resources:
        mem_mb = 10000
    shell:
        """
        perl {params.homer_cmd} {params.mfile} {params.genome} > {output}
        """

rule merge_homer_motifs:
    input:
        lambda wildcards: expand(os.path.join(HOMERFILES_DIR, '{tf}','scanMotifsGenomeWide.{motif_file}.txt'), tf = set( set(TF_list)), motif_file = homer_motifs_dict[wildcards.tf])
    output:
        os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt')
    message: "working on {wildcards}" 
    params:
        jobname = '{tf}'
    resources:
        mem_mb = 10000
    run:
        with open(output[0], 'w') as outfile:
            for fname in input:
                with open(fname) as infile:
                    for i, line in enumerate(infile):
                        outfile.write(line)





































# rule create_training_set:
#     input: 
#         rules.merge_homer_motifs.output
#         #os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt')
#     output: 
#         f1=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'),
#         f2=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'),
#         f3=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.info.txt.gz'),
#         f4=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.summary.txt')
#     params:
#         rscript = config['rscript'],
#         bedfiles_dir=config['cistrome_data_dir'],
#         tf_tissue = '{tf}_{tissue}',
#         bedlinks_dir = os.path.join(BEDLINKS_DIR, '{tf}_{tissue}'),
#         sortedbeds_dir = os.path.join(SORTEDBEDS_DIR, '{tf}_{tissue}'),
#         cistrome_mtdt = config['TF_table'],
#         jobname = '{tf}_{tissue}',
#         basename = os.path.join(PREDICTORS_DIR, '{tf}_{tissue}'),
#         train_by_chromosome = config['train_by_chromosome'],
#         dcids = config['dcids']
#         #nfiles = len([name for name in os.listdir(os.path.join(BEDLINKS_DIR, f'{tf}_{tissue}')) if os.path.isfile(name)])
#     message: "working on {wildcards}"
#     resources:
#         partition = "beagle3", #if params.nfiles > 200 else "caslake",
#         #attempt: attempt * 100,
#         mem_cpu = 16, #lambda wildcards, attempt: attempt * 8,
#         nodes = 1,
#         # mem_cpu = 4, #lambda wildcards, attempt: attempt * 8,
#         # nodes = 1,
#         #load= 50 #if resources.partition == 'bigmem' else 1
#     threads: 8
#     shell:
#         """
#         {params.rscript} workflow/src/create_training_sets_bedtools.R --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --predicted_motif_file {input} --sorted_bedfiles_directory {params.sortedbeds_dir} --bedlinks_directory {params.bedlinks_dir} --predictors_file {output.f1} --ground_truth_file {output.f2} --info_file {output.f3} --cistrome_metadata_file {params.cistrome_mtdt} --dcids {params.dcids} --train_by_chromosome {params.train_by_chromosome} --summary_file {output.f4}; sleep 5
#         """

# rule create_enformer_configuration:
#     input: rules.create_training_set.output.f1
#     output: os.path.join(PREDICTION_PARAMS_DIR, f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json')
#     message: "working on {wildcards}"
#     resources:
#         partition="beagle3"
#     params:
#         rscript = config['rscript'],
#         bdirectives = config['enformer']['base_directives'],
#         dset = config['dataset'],
#         model = config['enformer']['model'],
#         fasta_file = config['genome']['fasta'],
#         pdir = config['scratch_dir'], #DATA_DIR,
#         ddate = config['date'],
#         jobname = '{tf}_{tissue}'
#     shell:
#         """
#             {params.rscript} workflow/src/create_enformer_config.R --dataset {params.dset} --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output} --date {params.ddate}
#         """

# rule predict_with_enformer:
#     input:
#         rules.create_enformer_configuration.output
#     output:
#         os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json')
#     resources:
#         partition="beagle3",
#         time="04:00:00",
#         gpu=4,
#         mem_cpu=8,
#         cpu_task=8
#     params:
#         jobname = '{tf}_{tissue}',
#         enformer_predict_script = config['enformer']['predict']
#     message: 
#         "working on {params.jobname}"
#     shell:
#         """
#             python3 {params.enformer_predict_script} --parameters {input}
#         """

# rule aggregate_predictions:
#     input:
#         rules.predict_with_enformer.output
#     output:
#         os.path.join(AGGREGATION_DIR, f'{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz')
#     message: 
#         "working on {wildcards}"
#     resources:
#         partition="beagle3",
#         mem_cpu=8,
#         cpu_task=8,
#         mem_mb=24000
#     params:
#         jobname = '{tf}_{tissue}',
#         aggregation_script = config['enformer']['aggregate'],
#         aggtype = config['enformer']['aggtype'],
#         output_folder = AGGREGATION_DIR,
#         hpc = "beagle3",
#         parsl_executor = "local",
#         delete_enformer_outputs = config["delete_enformer_outputs"]
#     run:
#         if params.delete_enformer_outputs == True:
#             shell("python3 {params.aggregation_script} --metadata_file {input} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor} --delete_enformer_outputs")
#         elif params.delete_enformer_outputs == False: # don't delete the outputs
#             shell("python3 {params.aggregation_script} --metadata_file {input} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor}")

# rule prepare_training_data:
#     input:
#         p1 = rules.aggregate_predictions.output,
#         p2 = rules.create_training_set.output.f2
#     output:
#         p1=os.path.join(AGGREGATION_DIR, f'train_{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.prepared.csv.gz'),
#         p2=os.path.join(AGGREGATION_DIR, f'test_{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.prepared.csv.gz')
#     message: 
#         "preparing {wildcards} training and test data"
#     params:
#         jobname = '{tf}_{tissue}',
#         rscript = config['rscript'],
#         aggtype = config['enformer']['aggtype']
#     resources:
#         mem_mb=24000,
#         partition="beagle3"
#     shell:
#         """
#             {params.rscript} workflow/src/train_test_split.R --data_file {input.p1} --ground_truth_file {input.p2} --aggregation_method {params.aggtype} --train_prepared_file {output.p1} --test_prepared_file {output.p2}
#         """

# rule train_TFPred_weights:
#     input: rules.prepare_training_data.output.p1
#     output: 
#         mlogistic=os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.logistic.rds'),
#         mlinear=os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.rds')
#     message:
#         "training on {wildcards} training data"
#     params:
#         jobname = '{tf}_{tissue}',
#         rscript = config['rscript'],
#         nfolds=5,
#         basename=os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}')
#     resources:
#         #mem_mb=100000,
#         mem_mb = get_mem_mb_allocations,
#         partition = get_cluster_allocation,
#         cpu_task=12,
#         mem_cpu=12,
#         load=50
#     shell:
#         """
#             {params.rscript} workflow/src/train_enet.R --train_data_file {input} --rds_file {params.basename} --nfolds {params.nfolds}; sleep 12
#         """

# rule evaluate_TFPred:
#     input: 
#         linear_model = rules.train_TFPred_weights.output.mlinear,
#         logistic_model = rules.train_TFPred_weights.output.mlogistic,
#         train_data = rules.prepare_training_data.output.p1,
#         test_data = rules.prepare_training_data.output.p2
#     output: 
#         os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.logistic.train_eval.txt.gz'),
#         os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.logistic.test_eval.txt.gz'),
#         os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.train_eval.txt.gz'),
#         os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.test_eval.txt.gz')
#     message:
#         "evaluating on {wildcards} training and test data"
#     params:
#         jobname = '{tf}_{tissue}',
#         rscript = config['rscript'],
#         basename=os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}')
#     resources:
#         mem_mb= 100000,
#         partition="beagle3"
#     shell:
#         """
#             {params.rscript} workflow/src/evaluate_enet.R --linear_model {input.linear_model} --logistic_model {input.logistic_model} --train_data_file {input.train_data} --test_data_file {input.test_data} --eval_output {params.basename}
#         """




# print(TF_list)
# print(tissue_list)

# add check to match all these dictionaries

# homer_data = os.path.join(config['homer']['dir'], 'data/knownTFs/motifs')
# #hpath = os.path.join('/project2/haky/temi/software/homer', 'data/knownTFs/motifs')
# linked_homer_motifs = [module.link_homer_motifs(TF=d[0], tissue=d[1], from_db=homer_data, to_db=os.path.join(HOMERFILES_DIR, d[0])) for d in valid_TFs]
# TF_tissue_motifs_dict = {k: v for elem in linked_homer_motifs if elem is not None for k, v in elem.items()}
# #print(TF_tissue_motifs_dict)

# # link bed files
# #cistrome_mtdt = '/project2/haky/Data/TFXcan/cistrome/raw/human_factor_full_QC.txt' 
# # cistrome_df = pd.read_table(config['TF_table'])
# # linked_bedfiles = [module.link_cistrome_bedfiles(TF=d[0], tissue=d[1], from_db=config['cistrome_data_dir'], to_db=os.path.join(BEDLINKS_DIR, '_'.join(d)).replace(' ', '-'), db_df=cistrome_df) for d in valid_TFs]
# #homer_motifs_wildcards = glob_wildcards(os.path.join(HOMERFILES_DIR, '{tfs}/{motif_files}.motif'))
# dts = zip(homer_motifs_wildcards.tfs, homer_motifs_wildcards.motif_files) # based on what has been linked
# motif_grouping_dict = module.group_tf_motif_files(dts)
# valid_chromosomes = [f'chr{i}' for i in range(1,23)] + ['chrX']

# TF_list = [d[0] for d in valid_TFs]
# tissue_list = [d[1].replace(' ', '-') for d in valid_TFs]

# # print(expand(os.path.join(PREDICTORS_DIR, '{tf_tissue}_predictors.txt'), tf_tissue=TF_tissue_motifs_dict.keys()))

# # def count_number_of_observations(ffile):
# #     with open(ffile, 'r') as fp:
# #         x = len(fp.readlines())
# #         return(x)

# # checkpoint count_number_of_observations:
# #     input: rules.create_training_set.output.f1
# #     output:

# # report: 
# #     f"reports/workflow_{config['date']}.rst"