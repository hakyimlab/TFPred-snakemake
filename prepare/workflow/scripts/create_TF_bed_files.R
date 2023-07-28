

#Rscript -e 'install.packages(c("yaml", "data.table", "glue", "tidyverse"), repos="https://cloud.r-project.org")'

# TF <- 'FOXA1'
# tissue_type <- 'Breast'

args <- commandArgs(trailingOnly = TRUE)

TF <- args[1]
tissue_type <- args[2]
output_dir <- args[3]
config_file <- args[4]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--transcription_factor", help='A list of files to combine'),
    make_option("--tissue", help='A list of files to combine'),
    make_option("--output_directory", help='The final output file.'),
	make_option("--configuration_file", help='The final output file.')
)

opt <- parse_args(OptionParser(option_list=option_list))

# print(TF)
# print(tissue_type)
# print(output_dir)

library(yaml)
library(data.table, quietly=TRUE, warn.conflicts = TRUE)
library(glue)
suppressMessages(library(tidyverse, quietly=TRUE, warn.conflicts = TRUE))

#setwd('/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/pipelines')
cfile <- '/project2/haky/temi/projects/TFPred-snakemake/config/pipeline.yaml'
directives <- yaml::yaml.load_file(cfile)
tf_db <- data.table::fread(directives$TF_table)
# tf_db %>% dplyr::filter(Factor == TF) %>% dplyr::group_by(Cell_line) %>% summarise(n_files=n()) %>% arrange(desc(n_files))

#base::subset(x=tf_db, subset = Factor == TF)

# get the ids 
TF_data <- base::subset(x=tf_db, subset = (Factor == TF) & (Tissue_type == tissue_type))
# 
if( ! (nrow(TF_data) > 1 & ncol(TF_data) > 1)){
    stop('')
} else {
    TF_files <- list.files(glue('{directives$data_dir}/human_factor'), pattern=paste0('^', TF_data$DCid, collapse='_*|'), full.names=T)
}


if(!dir.exists(glue('{output_dir}'))){
    dir.create(glue('{output_dir}'))
}

out <- sapply(TF_files, function(each_file){
    
    output_file <- normalizePath(output_dir)
    bname <- basename(each_file)

    output_file <- file.path(output_file, bname)

    cmd <- glue('ln -s {each_file} {output_file}')
    print(cmd)
    #print(base::normalizePath(output_dir))
    #to_ = glue("{normalizePath(output_dir)}/{basename(each_file)}")
    #base::file.symlink(from=each_file, to=to_)
})

# find homer motifs
# tf_lower <- tolower(TF)
# potential_motif_files <- list.files(glue('{directives$homer$dir}/data/knownTFs/motifs'), glue('^{tf_lower}'), full.names=T)[1]


#glue("{directives$project_dir}/scripts/utilities/scan_for_motifs.sh")

