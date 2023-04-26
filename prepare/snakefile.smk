# Description: Common snakefile
# Author: Temi
# Date: Wed Mar 29 2023

import pandas as pd
import os, glob, sys, re
from snakemake.io import glob_wildcards

sys.path.append('workflow/scripts')
sys.path.append('modules')

import module

print(os.getcwd())
print_progress = False

cfile = "config/pipeline.yaml"
configfile: cfile

#metadata_dt = pd.read_csv("/project2/haky/temi/projects/TFPred-snakemake/metadata/metadata.txt")   #config['metadata'])

# if not os.path.isdir(config['enformer']['prediction_directives']['metadata_dir']):
#     os.makedirs(config['enformer']['prediction_directives']['metadata_dir'])

config_file_path = os.path.abspath(cfile)

# print(config_file_path)
# print(metadata_dt)

# directories
DATA_DIR = 'data'
BEDFILES_DIR = os.path.join(DATA_DIR, 'bed_files')
HOMERFILES_DIR = os.path.join(DATA_DIR, 'homer_files')
PREDICTORS_DIR = os.path.join(DATA_DIR, 'predictor_files')
PREDICTION_PARAMS_DIR = os.path.join(DATA_DIR, 'prediction_parameters')
METADATA_DIR = 'metadata'

metadata_dt = pd.read_csv(config['metadata'])
details = []
for row in metadata_dt.itertuples():
    #print(row)
    r = [row.transcription_factor, row.tissue]
    details.append(r)
print(details)

# check if the TF is present in the homer and cistrome data_db columns
data_db = pd.read_table('./info/data_db.txt')
data_db = data_db.loc[(data_db['homer_db'] == 1) & (data_db['cistrome_db'] == 1)]
valid_TFs = [d for d in details if d[0] in data_db.TF.tolist()] # filter the lists

# write out the list of unavailble TF-tissue pairs
invalid_TFs = [d for d in details if d[0] not in data_db.TF.tolist()]
if len(invalid_TFs) > 0:
    print(f'Some invalid TF-tissue pairs were found. See {os.path.join(METADATA_DIR, "invalid_TFs.txt")} for details.')
    pd.DataFrame(invalid_TFs).to_csv(os.path.join(METADATA_DIR, 'invalid_TFs.csv'), sep=',', index=False, header=['transcription_factor', 'tissue'])

if len(valid_TFs) == 0:
    print('No valid TF-tissue pairs were found. Please check the metadata file.')
    sys.exit(1)
else:
    print(f'{len(valid_TFs)} valid TF-tissue pairs were found.')
    pd.DataFrame(valid_TFs).to_csv(os.path.join(METADATA_DIR, 'valid_TFs.csv'), sep=',', index=False, header=['transcription_factor', 'tissue'])

homer_data = os.path.join(config['homer']['dir'], 'data/knownTFs/motifs')
#hpath = os.path.join('/project2/haky/temi/software/homer', 'data/knownTFs/motifs')

linked_homer_motifs = [module.link_homer_motifs(TF=d[0], tissue=d[1], from_db=homer_data, to_db=os.path.join(HOMERFILES_DIR, d[0])) for d in valid_TFs]
TF_tissue_motifs_dict = {k: v for elem in linked_homer_motifs if elem is not None for k, v in elem.items()}
#print(TF_tissue_motifs_dict)

# link bed files
#cistrome_mtdt = '/project2/haky/Data/TFXcan/cistrome/raw/human_factor_full_QC.txt' 
cistrome_df = pd.read_table(config['TF_table'])
linked_bedfiles = [module.link_cistrome_bedfiles(TF=d[0], tissue=d[1], from_db=config['cistrome_data_dir'], to_db=os.path.join(BEDFILES_DIR, '_'.join(d)).replace(' ', '-'), db_df=cistrome_df) for d in valid_TFs]
homer_motifs_wildcards = glob_wildcards(os.path.join(HOMERFILES_DIR, '{tfs}/{motif_files}.motif'))
dts = zip(homer_motifs_wildcards.tfs, homer_motifs_wildcards.motif_files) # based on what has been linked
motif_grouping_dict = module.group_tf_motif_files(dts)
valid_chromosomes = [f'chr{i}' for i in range(1,23)] + ['chrX']

TF_list = [d[0] for d in valid_TFs]
tissue_list = [d[1].replace(' ', '-') for d in valid_TFs]

print(expand(os.path.join(PREDICTORS_DIR, '{tf_tissue}_predictors.txt'), tf_tissue=TF_tissue_motifs_dict.keys()))

rule all:
    input:
        os.path.join('metadata', 'valid_TFs.csv'),
        expand(os.path.join(HOMERFILES_DIR, '{tf}', '{motif_file}.motif'), zip, tf = homer_motifs_wildcards.tfs, motif_file = homer_motifs_wildcards.motif_files),
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide_{motif_file}.txt'), zip, tf = homer_motifs_wildcards.tfs, motif_file = homer_motifs_wildcards.motif_files),
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt'), tf = set(homer_motifs_wildcards.tfs)),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}_predictors.txt'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}_ground_truth.txt'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTION_PARAMS_DIR, f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf=TF_list, tissue=tissue_list)

rule find_homer_motifs:
    input:
        os.path.join(HOMERFILES_DIR, '{tf}', '{motif_file}.motif')
    output:
        os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide_{motif_file}.txt')
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
        lambda wildcards: expand(os.path.join(HOMERFILES_DIR, '{{tf}}','scanMotifsGenomeWide_{motif_file}.txt'), tf = set(motif_grouping_dict.keys()), motif_file = motif_grouping_dict[wildcards.tf])
    output:
        os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt')
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
    input: 
        rules.merge_homer_motifs.output
        #os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt')
    output: 
        #f0=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}_predictor_regions.txt'),
        f1=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}_predictors.txt'),
        f2=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}_ground_truth.txt'),
        f3=os.path.join(PREDICTION_PARAMS_DIR, f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json')
    params:
        rscript = config['rscript'],
        tf_tissue = '{tf}_{tissue}',
        bedfiles_dir = os.path.join(BEDFILES_DIR, '{tf}_{tissue}'),
        config_file = cfile,
        jobname = '{tf}_{tissue}',
        basename = os.path.join(PREDICTORS_DIR, '{tf}_{tissue}')
    message: "working on {wildcards}"
    resources:
        mem_mb= 100000
    threads: 32
    shell:
        """
        {params.rscript} prepare/workflow/scripts/create_training_sets.R {params.tf_tissue} {input} {params.bedfiles_dir} {output.f1} {output.f2} {output.f3} {params.config_file}
        """