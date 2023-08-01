# Author: Temi
# Date: Thursday July 27 2023
# Description: script used to create enformer predict parameters file
# Usage: Rscript create_training_sets.R [OPTIONS]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--transcription_factor", help='A transcription factor e.g. AR'),
    make_option("--tissue", help='A tissue e.g. Breast'),
    make_option("--dataset", help="The name of the dataset"),
    make_option("--base_directives", help='A yaml file containing directives for enformer; will be used to create a json file for the enformer predict pipeline'),
    make_option("--project_directory", help='A project directory for enformer predict'),
	make_option("--predictors_file", help='predictor file containing the intervals to predict on'),
    make_option("--model", help='enformer model'),
	make_option("--fasta_file", help='fasta file, typically hg38'),
    make_option("--parameters_file", help='the json file that will be created')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(yaml)

# read and write the enformer config file
directives <- yaml::yaml.load_file(opt$base_directives)
enformer_parameters_json <- directives$enformer$prediction_directives
# you may change these as appropriate
enformer_parameters_json[['project_dir']] <- normalizePath(opt$project_directory)
enformer_parameters_json[["interval_list_file"]] <- opt$predictors_file
enformer_parameters_json[['prediction_data_name']] <- opt$dataset
enformer_parameters_json[['prediction_id']] <- paste(opt$transcription_factor, opt$tissue, sep='_')
#enformer_parameters_json[['date']] <- Sys.Date()
enformer_parameters_json[['model_path']] <- opt$model
enformer_parameters_json[['fasta_file']] <- opt$fasta_file

# chANGE the metadata dir
enformer_parameters_json[['metadata_dir']] <- dirname(opt$parameters_file)

write(
    jsonlite::toJSON(enformer_parameters_json, na='null', pretty=TRUE, auto_unbox=T),
    file=opt$parameters_file
)


# /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/bin/Rscript prepare/workflow/scripts/create_enformer_config.R --transcription_factor AR --tissue Breast --base_directives config/enformer_base.yaml --project_directory data/predictions_folder --predictors_file data/predictor_files/AR_Breast.predictors.txt --model "/project2/haky/Data/enformer/raw" --fasta_file "/project2/haky/Data/hg_sequences/hg38/Homo_sapiens_assembly38.fasta" --parameters_file data/prediction_parameters/enformer_parameters_cistrome_AR_Breast.json