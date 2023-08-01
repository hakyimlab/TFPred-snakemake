# Description: Common snakefile
# Author: Temi
# Date: Wed Mar 29 2023

import pandas as pd
import os, glob, sys, re
from snakemake.io import glob_wildcards

sys.path.append('workflow/src')
sys.path.append('modules')

import module
print_progress = False

# directories
DATA_DIR = 'data'
BEDLINKS_DIR = os.path.join(DATA_DIR, 'bed_links')
HOMERFILES_DIR = os.path.join(DATA_DIR, 'homer_files')
PREDICTORS_DIR = os.path.join(DATA_DIR, 'predictor_files')
PREDICTION_PARAMS_DIR = os.path.join(DATA_DIR, 'prediction_parameters')
METADATA_DIR = 'metadata'
PREDICTIONS_DIR = os.path.join(DATA_DIR, 'predictions_folder')
AGGREGATION_DIR = os.path.join(DATA_DIR, 'aggregation_folder')
MODELS_DIR = 'output/models'
MODELS_EVAL_DIR = 'output/models_eval'

metadata_dt = pd.read_csv(config['metadata'])
print(metadata_dt)
details = []
for row in metadata_dt.itertuples():
    #print(row)
    r = [row.assay, row.context]
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
    pd.DataFrame(invalid_TFs).to_csv(os.path.join(METADATA_DIR, 'invalid_TFs.csv'), sep=',', index=False, header=['assay', 'context'])

if len(valid_TFs) == 0:
    print('No valid TF-tissue pairs were found. Please check the metadata file.')
    sys.exit(1)
else:
    print(f'{len(valid_TFs)} valid TF-tissue pairs were found.')
    pd.DataFrame(valid_TFs).to_csv(os.path.join(METADATA_DIR, 'valid_TFs.csv'), sep=',', index=False, header=['assay', 'context'])

#use this to filter for homer files
homer_data = os.path.join(config['homer']['dir'], 'data/knownTFs/motifs')
#hpath = os.path.join('/project2/haky/temi/software/homer', 'data/knownTFs/motifs')

linked_homer_motifs = [module.link_homer_motifs(TF=d[0], tissue=d[1], from_db=homer_data, to_db=os.path.join(HOMERFILES_DIR, d[0])) for d in valid_TFs]
TF_tissue_motifs_dict = {k: v for elem in linked_homer_motifs if elem is not None for k, v in elem.items()}
#print(TF_tissue_motifs_dict)

# link bed files
#cistrome_mtdt = '/project2/haky/Data/TFXcan/cistrome/raw/human_factor_full_QC.txt' 
cistrome_df = pd.read_table(config['TF_table'])
linked_bedfiles = [module.link_cistrome_bedfiles(TF=d[0], tissue=d[1], from_db=config['cistrome_data_dir'], to_db=os.path.join(BEDLINKS_DIR, '_'.join(d)).replace(' ', '-'), db_df=cistrome_df) for d in valid_TFs]
homer_motifs_wildcards = glob_wildcards(os.path.join(HOMERFILES_DIR, '{tfs}/{motif_files}.motif'))
dts = zip(homer_motifs_wildcards.tfs, homer_motifs_wildcards.motif_files) # based on what has been linked
motif_grouping_dict = module.group_tf_motif_files(dts)
valid_chromosomes = [f'chr{i}' for i in range(1,23)] + ['chrX']

TF_list = [d[0] for d in valid_TFs]
tissue_list = [d[1].replace(' ', '-') for d in valid_TFs]

# print(expand(os.path.join(PREDICTORS_DIR, '{tf_tissue}_predictors.txt'), tf_tissue=TF_tissue_motifs_dict.keys()))

rule all:
    input:
        os.path.join('metadata', 'valid_TFs.csv'),
        expand(os.path.join(HOMERFILES_DIR, '{tf}', '{motif_file}.motif'), zip, tf = homer_motifs_wildcards.tfs, motif_file = homer_motifs_wildcards.motif_files),
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide_{motif_file}.txt'), zip, tf = homer_motifs_wildcards.tfs, motif_file = homer_motifs_wildcards.motif_files),
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt'), tf = set(homer_motifs_wildcards.tfs)),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.info.txt.gz'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTION_PARAMS_DIR, f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'train_{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'test_{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.rds'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.rds'), zip, tf = TF_list, tissue = tissue_list)

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
        f1=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'),
        f2=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'),
        f3=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.info.txt.gz')
    params:
        rscript = config['rscript'],
        bedfiles_dir=config['cistrome_data_dir'],
        tf_tissue = '{tf}_{tissue}',
        bedlinks_dir = os.path.join(BEDLINKS_DIR, '{tf}_{tissue}'),
        cistrome_mtdt = config['TF_table'],
        jobname = '{tf}_{tissue}',
        basename = os.path.join(PREDICTORS_DIR, '{tf}_{tissue}')
    message: "working on {wildcards}"
    resources:
        partition = 'beagle3',
        mem_mb= 150000,
        nodes=1
    threads: 32
    shell:
        """
        {params.rscript} workflow/src/create_training_sets.R --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --predicted_motif_file {input} --bedfiles_directory {params.bedfiles_dir} --bedlinks_directory {params.bedlinks_dir} --predictors_file {output.f1} --ground_truth_file {output.f2} --info_file {output.f3} --cistrome_metadata_file {params.cistrome_mtdt} ; sleep 5
        """

rule create_enformer_configuration:
    input: rules.create_training_set.output.f1
    output: os.path.join(PREDICTION_PARAMS_DIR, f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json')
    message: "working on {wildcards}"
    resources:
        partition="beagle3"
    params:
        rscript = config['rscript'],
        bdirectives = config['enformer']['base_directives'],
        dset = config['dataset'],
        model = config['enformer']['model'],
        fasta_file = config['genome']['fasta'],
        pdir = DATA_DIR,
        jobname = '{tf}_{tissue}'
    shell:
        """
            {params.rscript} workflow/src/create_enformer_config.R --dataset {params.dset} --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output}
        """

rule predict_with_enformer:
    input:
        rules.create_enformer_configuration.output
    output:
        os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{config["dataset"]}_{{tf}}_{{tissue}}.json')
    resources:
        slurm="gres=gpu:16",
        partition="beagle3",
        time="04:00:00",
        nodes=4,
        mem_mb = 150000
    params:
        jobname = '{tf}_{tissue}',
        enformer_predict_script = config['enformer']['predict']
    message: 
        "working on {params.jobname}"
    shell:
        """
            python3 {params.enformer_predict_script} --parameters {input}
        """

rule aggregate_predictions:
    input:
        rules.predict_with_enformer.output
    output:
        os.path.join(AGGREGATION_DIR, f'{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz')
    message: 
        "working on {wildcards}"
    resources:
        partition="beagle3"
    params:
        jobname = '{tf}_{tissue}',
        aggregation_script = config['enformer']['aggregate'],
        aggtype = config['enformer']['aggtype'],
        output_folder = AGGREGATION_DIR,
        hpc = "beagle3",
        parsl_executor = "local",
        delete_enformer_outputs = False
    resources:
        mem_mb= 100000
    run:
        if params.delete_enformer_outputs == True:
            shell("python3 {params.aggregation_script} --metadata_file {input} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor} --delete_enformer_outputs")
        elif params.delete_enformer_outputs == False: # don't delete the outputs
            shell("python3 {params.aggregation_script} --metadata_file {input} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor}")

rule prepare_training_data:
    input:
        p1 = rules.aggregate_predictions.output,
        p2 = rules.create_training_set.output.f2
    output:
        p1=os.path.join(AGGREGATION_DIR, f'train_{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.prepared.csv.gz'),
        p2=os.path.join(AGGREGATION_DIR, f'test_{config["dataset"]}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.prepared.csv.gz')
    message: 
        "preparing {wildcards} training and test data"
    params:
        jobname = '{tf}_{tissue}',
        rscript = config['rscript'],
        aggtype = config['enformer']['aggtype']
    resources:
        mem_mb= 100000,
        partition="beagle3"
    shell:
        """
            {params.rscript} workflow/src/train_test_split.R --data_file {input.p1} --ground_truth_file {input.p2} --aggregation_method {params.aggtype} --train_prepared_file {output.p1} --test_prepared_file {output.p2}
        """

rule train_TFPred_weights:
    input: rules.prepare_training_data.output.p1
    output: os.path.join(MODELS_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.rds')
    message:
        "training on {wildcards} training data"
    params:
        jobname = '{tf}_{tissue}',
        rscript = config['rscript'],
        nfolds=5
    resources:
        mem_mb= 100000,
        partition="beagle3"
    shell:
        """
            {params.rscript} workflow/src/train_enet.R --train_data_file {input} --rds_file {output} --nfolds {params.nfolds}; sleep 12
        """

rule evaluate_TFPred:
    input: 
        model_rds = rules.train_TFPred_weights.output,
        train_data = rules.prepare_training_data.output.p1,
        test_data = rules.prepare_training_data.output.p2
    output: os.path.join(MODELS_EVAL_DIR, f"{config['dataset']}_{{tf}}_{{tissue}}_{config['date']}", f'{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.linear.rds')
    message:
        "evaluating on {wildcards} training and test data"
    params:
        jobname = '{tf}_{tissue}',
        rscript = config['rscript']
    resources:
        mem_mb= 100000,
        partition="beagle3"
    shell:
        """
            {params.rscript} workflow/src/evaluate_enet.R --model {input.model_rds} --train_data_file {input.train_data} --test_data_file {input.test_data} --eval_output {output}; sleep 10
        """