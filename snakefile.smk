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
print(motif_inputs)
print(motif_outputs)

def catch(func, *args, handle=lambda e : e, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception as e:
        return(None)
        #return handle(e)


def count_number_of_lines(wildcards):
    # tf, tissue, = os.path.join(PREDICTORS_DIR, "{wildcards.tf}_{wildcards.tissue}.predictors.txt")
    # chckd = rules.create_training_set.output.f1 #get(**wildcards)
    chckd = os.path.join(PREDICTORS_DIR, f"{wildcards.tf}_{wildcards.tissue}.predictors.txt")
    nlines = catch(int, subprocess.run(f"wc -l {chckd} | awk '{{print $1}}'", shell=True, capture_output = True).stdout.decode().strip())
    print(f'Checking {nlines}')
    if nlines is None:
        return(0)
    return(nlines)

rule all:
    input:
        motif_outputs,
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt'), tf = set(TF_list)),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'{runname}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'train_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'test_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.logistic.rds'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.linear.rds'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.train_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.train_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.test_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.test_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list)

rule find_homer_motifs:
    input:
        os.path.join(config["homer"]['motifs_database'], '{motif_file}')
    output: 
        os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide.{motif_file}.txt')
    params:
        run = run,
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
        #gatherMotifFiles
        lambda wildcards: expand(os.path.join(HOMERFILES_DIR, wildcards.tf, 'scanMotifsGenomeWide.{motif_file}.txt'), tf = set(homer_motifs_dict.keys()), motif_file = homer_motifs_dict[wildcards.tf])
    output:
        os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt')
    message: "working on {input}" 
    params:
        jobname = '{tf}',
        run = run,
        # ifiles = lambda wildcards: expand(os.path.join(HOMERFILES_DIR, '{tf}','scanMotifsGenomeWide.{motif_file}.txt'), tf = set(TF_list), motif_file = homer_motifs_dict[wildcards.tf]),
        #ifiles = gatherMotifFiles,

    resources:
        mem_mb = 10000
    run:
        with open(output[0], 'w') as outfile:
            for fname in input:
                with open(fname) as infile:
                    for i, line in enumerate(infile):
                        outfile.write(line)


rule create_training_set:
    input: rules.merge_homer_motifs.output
    #     lambda wildcards: os.path.join(HOMERFILES_DIR, wildcards.tf, 'merged_motif_file.txt')
    output:
        #directory(PREDICTORS_DIR)
        f1=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'),
        f2=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'),
        # f3=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.info.txt.gz'),
        # f4=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.summary.txt')
    params:
        run = run,
        rscript = config['rscript'],
        tf_tissue = '{tf}_{tissue}',
        bedfiles_dir = config['peaks']['directory'],
        sortedbeds_dir = os.path.join(SORTEDBEDS_DIR, '{tf}_{tissue}'),
        jobname = '{tf}_{tissue}',
        basename = os.path.join(PREDICTORS_DIR, '{tf}_{tissue}'),
        test_chromosomes = ','.join(config['train_by_chromosome']['test']),
        peaks_files = lambda wildcards: ','.join(model_config[wildcards.tf]["peakFiles"][wildcards.tissue]) if isinstance(model_config[wildcards.tf]["peakFiles"][wildcards.tissue], list) else model_config[wildcards.tf]["peakFiles"][wildcards.tissue],
        #f1=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'),
        # f2=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'),
        f3=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.info.txt.gz'),
        f4=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.summary.txt'),
        mfile = rules.merge_homer_motifs.output
    message: "working on {wildcards}"
    resources:
        partition = "caslake", #if params.nfiles > 200 else "caslake",
        #attempt: attempt * 100,
        mem_cpu = 4, #lambda wildcards, attempt: attempt * 8,
        nodes = 1,
        # mem_cpu = 4, #lambda wildcards, attempt: attempt * 8,
        # nodes = 1,
        #load= 50 #if resources.partition == 'bigmem' else 1
    threads: 8
    shell:
        """
        {params.rscript} workflow/src/create_training_sets_bedtools.R --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --predicted_motif_file {params.mfile} --sorted_bedfiles_directory {params.sortedbeds_dir} --bedlinks_directory {params.bedfiles_dir} --predictors_file {output.f1} --ground_truth_file {output.f2} --info_file {params.f3} --peaks_files {params.peaks_files} --test_chromosomes {params.test_chromosomes} --summary_file {params.f4}; sleep 5
        """

rule aggregate_epigenomes:
    input:
        #rules.create_training_set.output.f1
        #os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt')
        rules.create_training_set.output.f1
        # ifNotRunEnformer = collect_valid_predictor_files()
    output:
        os.path.join(AGGREGATION_DIR, f'{runname}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz')
    message: 
        "working on {wildcards}"
    resources:
        partition="caslake",
        mem_cpu=8,
        cpu_task=8,
        mem_mb=24000
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        enformer_epigenome_directory = config['enformer_epigenome_directory'],
        nlines = count_number_of_lines,
    run:
        if params.nlines < 100:
            shell("touch {output}")
        else:
            shell("python3 workflow/src/aggregate_epigenomes.py --loci_file {input} --reference_epigenome_dir {params.enformer_epigenome_directory} --output_file {output} --use_multiprocessing")


rule prepare_training_data:
    input:
        p1 = rules.aggregate_epigenomes.output,
        p2 = rules.create_training_set.output.f2
    output:
        p1=os.path.join(AGGREGATION_DIR, f'train_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz'),
        p2=os.path.join(AGGREGATION_DIR, f'test_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz')
    message: 
        "preparing {wildcards} training and test data"
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        rscript = config['rscript'],
        aggtype = config['enformer']['aggtype'],
        nlines = count_number_of_lines
    resources:
        mem_mb=24000,
        partition="caslake",
    run:
        if params.nlines < 100:
            shell("touch {output.p1} && touch {output.p2}")
        else:
            shell("{params.rscript} workflow/src/train_test_split.R --data_file {input.p1} --ground_truth_file {input.p2} --aggregation_method {params.aggtype} --train_prepared_file {output.p1} --test_prepared_file {output.p2}")

rule train_TFPred_weights:
    input: rules.prepare_training_data.output.p1
    output: 
        mlogistic=os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.logistic.rds'),
        mlinear=os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.linear.rds')
    message:
        "training on {wildcards} training data"
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        rscript = config['rscript'],
        nfolds=5,
        basename=os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}'),
        nlines = count_number_of_lines
    resources:
        #mem_mb=100000,
        mem_mb = helpers.get_mem_mb_allocations,
        partition = helpers.get_cluster_allocation,
        cpu_task=12,
        mem_cpu=12,
        load=50
    run:
        if params.nlines < 100:
            shell("touch {output.mlogistic} && touch {output.mlinear}")
        else:
            shell("{params.rscript} workflow/src/train_enet.R --train_data_file {input} --rds_file {params.basename} --nfolds {params.nfolds}; sleep 12")

rule evaluate_TFPred:
    input: 
        linear_model = rules.train_TFPred_weights.output.mlinear,
        logistic_model = rules.train_TFPred_weights.output.mlogistic,
        train_data = rules.prepare_training_data.output.p1,
        test_data = rules.prepare_training_data.output.p2
    output: 
        os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.train_eval.txt.gz'),
        os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.test_eval.txt.gz'),
        os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.train_eval.txt.gz'),
        os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.test_eval.txt.gz')
    message:
        "evaluating on {wildcards} training and test data"
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        rscript = config['rscript'],
        basename=os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}'),
        nlines = count_number_of_lines
    resources:
        mem_mb= 100000,
        partition="caslake"
    run:
        if params.nlines < 100:
            shell("touch {output}")
        else:
            shell("{params.rscript} workflow/src/evaluate_enet.R --linear_model {input.linear_model} --logistic_model {input.logistic_model} --train_data_file {input.train_data} --test_data_file {input.test_data} --eval_output {params.basename}")


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


# rule aggregate_predictions:
#     input:
#         ifRunEnformer = rules.predict_with_enformer.output
#     output:
#         os.path.join(AGGREGATION_DIR, f'{runname}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz')
#     message: 
#         "working on {wildcards}"
#     resources:
#         partition="caslake",
#         mem_cpu=8,
#         cpu_task=8,
#         mem_mb=24000
#     params:
#         run = run,
#         jobname = '{tf}_{tissue}',
#         aggregation_script = config['enformer']['aggregate'],
#         aggtype = config['enformer']['aggtype'],
#         output_folder = AGGREGATION_DIR,
#         hpc = "beagle3",
#         parsl_executor = "local",
#         delete_enformer_outputs = config["delete_enformer_outputs"]
#     run:
#         if params.delete_enformer_outputs == True:
#             shell("python3 {params.aggregation_script} --metadata_file {input.ifRunEnformer} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor} --delete_enformer_outputs")
#         elif params.delete_enformer_outputs == False: # don't delete the outputs
#             shell("python3 {params.aggregation_script} --metadata_file {input.ifRunEnformer} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor}")



# rule create_enformer_configuration:
#     input: collect_valid_predictor_files() #checkpoints.create_training_set.get(**wildcards).output.f1 #lambda wildcards: checkpoints.create_training_set.get(tf=wildcards.tf, tissue=wildcards.tissue).output.f1 
#     output: os.path.join(PREDICTION_PARAMS_DIR, f'enformer_config_{run}.{{tf}}_{{tissue}}.json')
#     message: "working on {wildcards}"
#     resources:
#         partition="caslake"
#     params:
#         run = run,
#         rscript = config['rscript'],
#         bdirectives = config['enformer']['base_directives'],
#         dset = runname,
#         model = config['enformer']['model'],
#         fasta_file = config['genome']['fasta'],
#         pdir = config['scratch_dir'], #DATA_DIR,
#         ddate = rundate,
#         jobname = '{tf}_{tissue}',
#     shell:
#         """
#             {params.rscript} workflow/src/create_enformer_config.R --dataset {params.dset} --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output} --date {params.ddate}
#         """
#         #printf 'INFO: This is the input file %s\n' "{input}";

# rule predict_with_enformer:
#     input:
#         rules.create_enformer_configuration.output
#         #lambda wildcards: checkpoints.create_enformer_configuration.get(tf=wildcards.tf, tissue=wildcards.tissue).output
#     output:
#         os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{runname}_{{tf}}_{{tissue}}.json')
#     resources:
#         partition = "caslake", #"beagle3" if config["run_enformer"] is True else 'caslake',
#         time="00:30:00", #06:00:00" if config["run_enformer"] is True else "01:00:00",
#         gpu=0, #config['enformer']['ngpus'] if config["run_enformer"] is True else 0,
#         mem_cpu=8,
#         cpu_task=8
#     params:
#         run = run,
#         jobname = '{tf}_{tissue}',
#         enformer_predict_script = config['enformer']['predict']
#     message: 
#         "working on {params.jobname}"
#     run:
#         if config["run_enformer"] is True:
#             shell("python3 {params.enformer_predict_script} --parameters {input}")
#         elif config["run_enformer"] is False:
#             shell("touch {output}")


# print(f"========================================")
# print(motif_outputs)
# print(f"========================================")
# print(homer_motifs_dict)
#ruleorder: aggregate_epigenomes > aggregate_predictions

# if config["run_enformer"] == True:
#     ruleorder: aggregate_epigenomes > aggregate_predictions
# elif config["run_enformer"] == False:
#     ruleorder: aggregate_predictions > aggregate_epigenomes



# the function counts the lines in a file and returns the number of lines
# def collect_valid_predictor_files():
#     # tf, tissue, = os.path.join(PREDICTORS_DIR, "{wildcards.tf}_{wildcards.tissue}.predictors.txt")
#     #chckd = os.path.dirname(checkpoints.create_training_set.get(**wildcards).output.f1)
#     chckd = PREDICTORS_DIR
#     print(f'Checking {chckd}')
#     AA, BB, = glob_wildcards(os.path.join(chckd, '{tf}_{tissue}.predictors.txt'))
#     mmd = expand(os.path.join(chckd, "{tf}_{tissue}.predictors.txt"), tf = AA, tissue = BB)
#     # count the number of lines of each file
#     numlines = [subprocess.run(f"wc -l {f} | awk '{{print $1}}'", shell=True, capture_output = True).stdout.decode() for f in mmd]
#     numlines = [x.strip() for x in numlines]
#     numlines = [catch(int, x) for x in numlines]
#     assert len(numlines) == len(mmd), 'ERROR - [FATAL] The number of lines in the predictor files do not match the number of files. Exiting...'
#     # FILTER OUT nones and return
#     mmd = [m for m, n in zip(mmd, numlines) if n is not None]
#     return(mmd)

# def collect_valid_predictor_files(wildcards):
#     # tf, tissue, = os.path.join(PREDICTORS_DIR, "{wildcards.tf}_{wildcards.tissue}.predictors.txt")
#     #chckd = os.path.dirname(checkpoints.create_training_set.get(**wildcards).output.f1)
#     # chckd = checkpoints.create_training_set.get(**wildcards).output
#     # print(f'Checking {chckd}')
#     chckd = PREDICTORS_DIR
#     AA, BB, = glob_wildcards(os.path.join(chckd, f'{wildcards.tf}_{wildcards.tissue}.predictors.txt'))
#     mmd = expand(os.path.join(chckd, "{tf}_{tissue}.predictors.txt"), zip, tf = AA, tissue = BB)
#     #count the number of lines of each file
#     numlines = [subprocess.run(f"wc -l {f} | awk '{{print $1}}'", shell=True, capture_output = True).stdout.decode() for f in mmd]
#     numlines = [x.strip() for x in numlines]
#     numlines = [catch(int, x) for x in numlines]
#     assert len(numlines) == len(mmd), 'ERROR - [FATAL] The number of lines in the predictor files do not match the number of files. Exiting...'
#     # FILTER OUT nones and return
#     mmd = [m for m, n in zip(mmd, numlines) if n > 100]
#     print(numlines)
#     return(mmd)




    # if nlines > 0:
    #     return(chckd)
        #raise Exception(f'ERROR - Empty file in {chckd} found. Exiting...')


    # AA, BB, = glob_wildcards(os.path.join(chckd, '{tf}_{tissue}.predictors.txt'))
    # mmd = expand(os.path.join(chckd, "{tf}_{tissue}.predictors.txt"), tf = AA, tissue = BB)
    # # count the number of lines of each file
    # numlines = [subprocess.run(f"wc -l {f} | awk '{{print $1}}'", shell=True, capture_output = True).stdout.decode() for f in mmd]
    # numlines = [x.strip() for x in numlines]
    # numlines = [catch(int, x) for x in numlines]
    # assert len(numlines) == len(mmd), 'ERROR - [FATAL] The number of lines in the predictor files do not match the number of files. Exiting...'
    # # FILTER OUT nones and return
    # mmd = [m for m, n in zip(mmd, numlines) if n is not None]
    # return(mmd)


# def count_lines_and_filter(wildcards):
#     mmd, = glob_wildcards(os.path.join(PREDICTORS_DIR, "{mmd}.predictors.txt")).mmd
#     # count the number of lines
#     mmds_paths = [os.path.join(PREDICTION_PARAMS_DIR, "{md}.predictors.txt") for md in mmd]
#     wcmmds = with open("myfile.txt", "rb") as f:
#             num_lines = sum(1 for _ in f)


# def gatherMotifFiles(wildcards):
#     mfs = glob_wildcards(os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide.{motif_file}.txt')).motif_file
#     mfs = expand(os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide.{m}.txt'), m = mfs)
#     print(mfs)
#     return(mfs)