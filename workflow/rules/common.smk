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

configfile: "config/pipeline.yaml"
metadata_dt = pd.read_csv(config['metadata'])
rscript = config['rscript']
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

    bedfiles_ids = return_cistrome_bedfiles(TF=TF, tissue=tissue, cistrome_db_folder = config['data_dir'], cistrome_db = config['TF_table'])

    for key, values in bedfiles_ids.items():
        dst_folder = os.path.join('data', BEDFILES_DIR, key)
        if not os.path.isdir(dst_folder):
            os.makedirs(dst_folder)
        for value in values:
            dst_file = os.path.join(dst_folder, value)
            if not os.path.isfile(dst_file):
                os.symlink(os.path.join(config['data_dir'], value), dst_file)
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

print(homer_wildcards)


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

print(motif_grouping_dict)
print(tissue_grouping_dict)



# GROUPING_DICT = {}
# ALL_BEDFILES = []
# TF_FOLDERS = []
# for detail in details:
#     bedfiles_ids = return_cistrome_bedfiles(TF=detail[0], tissue=detail[1], cistrome_db_folder=config['data_dir'])
#     bedfiles_inputs = []
#     for key, values in bedfiles_ids.items():
#         for value in values:
#             bedfiles_inputs.append(f"{key}/{value}")
#             ALL_BEDFILES.append(f"{key}/{value}")
    
#     TF_FOLDERS.append('_'.join(detail))
#     nnames = f'{detail[0]}_{detail[1]}'
#     GROUPING_DICT[nnames] = bedfiles_inputs
# #print(TF_FOLDERS)

# #print(snakemake.io.expand(os.path.join('data', BEDFILES_DIR, '{tf_bedfile}'), tf_bedfile = ALL_BEDFILES))

rule all:
    input:
        expand(os.path.join('data', HOMERFILES_DIR, '{tf}', '{motif_file}.motif'), zip, tf = homer_wildcards.tfs, motif_file = homer_wildcards.motif_files),
        expand(os.path.join('data', HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide_{motif_file}.txt'), zip, tf = homer_wildcards.tfs, motif_file = homer_wildcards.motif_files),
        expand(os.path.join('data', HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt'), tf = set(homer_wildcards.tfs)),
        expand(os.path.join('data', 'predictor_files', '{tf}_{tissue}_predictors.txt'), zip, tf = TF_list, tissue = tissue_list)

rule find_homer_motifs:
    input:
        os.path.join('data', HOMERFILES_DIR, '{tf}', '{motif_file}.motif')
    output:
        os.path.join('data', HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide_{motif_file}.txt')
    params:
        homer_cmd = config['homer']['scanMotifsGenome'],
        genome = config['homer']['genome'],
        jobname = '{tf}'
    message: "working on {wildcards}" 
    resources:
        mem_mb = 10000
    shell:
        """
        perl {params.homer_cmd} {input} {params.genome} > {output}
        """

rule merge_homer_motifs:
    input:
        lambda wildcards: expand(os.path.join('data', HOMERFILES_DIR, '{{tf}}','scanMotifsGenomeWide_{motif_file}.txt'), tf = set(motif_grouping_dict.keys()), motif_file = motif_grouping_dict[wildcards.tf])
    output:
        os.path.join('data', HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt')
    message: "working on {wildcards} {input}" 
    params:
        jobname = '{tf}'
    resources:
        mem_mb = 10000
    run:
        #print(input)
        with open(output[0], 'w') as outfile:
            for fname in input:
                with open(fname) as infile:
                    for i, line in enumerate(infile):
                        outfile.write(line)

rule create_training_set:
    output: 
        os.path.join('data', 'predictor_files', '{TF}_{tissue}_predictors.txt')
    params:
        rscript = config['rscript'],
        TF = '{TF}',
        tissue = '{tissue}',
        jobname = '{TF}_{tissue}'
    message: "working on {wildcards}"
    resources:
        mem_mb = 10000
    shell:
        """
        {params.rscript} workflow/rules/create_training_sets.R {params.TF} {params.tissue} data/ {output}
        """

# SRA_FILES_BASENAMES = glob_wildcards(os.path.join("data/raw_samples/raw_fastq", "{SRA_FILES}.fastq.gz")).SRA_FILES

# if not os.path.isfile(metadata_file):
#     raise Exception("ERROR - Metadata sheet file does not exist in `metadata/`")

# # CREATING FILES AND LISTS AND PARAMETERS ===
# grouping_dict = module.grouping_information(metadata_file)
# SAMPLES = [l for l in list(grouping_dict.values()) for l in l]
# # you need to provide the path to the sra files
# SRA_FILES_BASENAMES = glob_wildcards(os.path.join(config["sra_folder"], "{SRA_FILES}.fastq.gz")).SRA_FILES
# # "data/raw_samples/raw_fastq"
# #SRA_FILES_BASENAMES = glob_wildcards(os.path.join("data/raw_samples/raw_fastq", "{SRA_FILES}.fastq.gz")).SRA_FILES
# sra_num_grouping = module.create_sra_num_dict(SRA_FILES_BASENAMES)
# trimming_grouping = {s: module.create_trim_inputs(sra_num_grouping, s) for s in SAMPLES}
# fastq_files = [l for l in list(trimming_grouping.values()) for l in l]
# individuals = list(grouping_dict.keys())
# # read group dictionary ===
# read_groups = module.return_read_group_information(metadata_file)

# # FILES AND DIRECTORIES ===
# project_name = config["project_name"]
# CONDA_YAML_FILE = config['conda_env']
# if not os.path.isfile(CONDA_YAML_FILE):
#     raise Exception("ERROR - conda yaml file cannot be found in `workflow/envs`")
# DATA_DIR = "data/raw_samples"
# GVCF_DIR = "data/gvcf_files"
# REBLOCKED_GVCF_DIR = "data/reblocked_gvcf_files"
# FINAL_VCFS_DIR =  "data/final_vcfs"
# fastqc_dir = "data/fastqc"
# LOG_DIR = os.path.join('logs', 'snakemake_log')
# TMP_DIR = os.path.join('/tmp', project_name)
# #BWA_INDEX = 'resource/alignment/index' # what is this?????????????
# known_variants_files = [f"--known-sites {os.path.join(config['variants_resource']['folder'], v)}" for v in config['variants_resource']['files']]
# genome_file = os.path.join(config['genome_resource']['folder'], config['genome_resource']['files'])
# #print(f'GENOME FILE - {genome_file}')
# #bwa_index_file = f"{config['bwa_resource']['folder']['files']}"
# bwa_index_file = os.path.join(config['bwa_resource']['folder'], config['bwa_resource']['files'][0])
# scratch_folder = os.path.join(config['scratch_folder'], project_name)
# variants_db = config["variants_db"]
# chromosomes = pd.read_csv(config["chromosomes"], header=None).iloc[:, 0].tolist()
# gmap_dir = config['genetic_map']['folder']

