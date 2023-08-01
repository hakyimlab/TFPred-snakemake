# Author: Temi
# Date: Thursday July 27 2023
# Description: script to evaluate TFPred models on train and test data
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--assay", help='of type `tf`, `dnase` or `histone`'),
    make_option("--context_type", help='of type `tissue`, `cell_type` or `cell_line`')
    make_option("--mt_file", help='a csv file of the assay type, and context'),
    make_option("--info_db", help='test data file'),
    make_option("--assay_db", help='evaluation file in with .rds extension'),
    make_option("--bedfiles_directory", help='path to the bed files')
)

opt <- parse_args(OptionParser(option_list=option_list))

opt <- list()
opt$assay <- 'tf'
opt$context_type <- 'tissue'
opt$mt_file <- '/project2/haky/temi/projects/TFPred-snakemake/metadata/mtfile.csv'
opt$info_db <- '/project2/haky/temi/projects/TFPred-snakemake/info/data_db.txt'
opt$assay_db <- '/project2/haky/Data/TFXcan/cistrome/raw/human_factor_full_QC.txt'
opt$bedfiles_directory <- '/project2/haky/Data/TFXcan/cistrome/raw/human_factor'

library(data.table)
library(tidyverse)
library(glue)

#
if(!(opt$assay) %in% c('tf', 'dnase', 'histone')){
    print(glue('ERROR - assay type must be one of {}'))
}

# read in the mt_file
mt <- data.table::fread(opt$mt_file) %>% dplyr::rename(assay=1, context=2) %>% as.data.frame()
assay_dt <- data.table::fread(opt$assay_db)
if(!is.na(opt$info_db)){
    info_dt <- data.table::fread(opt$info_db)
}

ctext <- switch(as.character(opt$context_type),
    'tissue' = 'Tissue_type',
    'cell_type' = 'Cell_type',
    'cell_line' = 'Cell_line',
    stop('Context type does not exist.')
)

get_dcids <- function(assay_name, context_name, assay_dt, ctext){
    dcids <- assay_dt %>% 
        dplyr::filter(base::get({{ctext}}) == context_name, Factor == assay_name) %>% pull(DCid)

    if(length(dcids) == 0){
        print(glue('INFO - No DCids found for {assay_name} and {context_name}'))
        quit(status = 1)
    } else {
        return(dcids) #bedfiles <- list.files(glue('{opt$bedfiles_directory}'), pattern=paste0('^', TF_data$DCid, collapse='_*|'), full.names=T)
    }
}

check_bedfiles <- function(dcids, bedfiles_directory){
    output_file <- normalizePath(bedfiles_directory)
    bedfiles <- list.files(glue('{bedfiles_directory}'), pattern=paste0('^', dcids, collapse='_*|'), full.names=F)
    out <- sapply(bedfiles, function(bf){return(file.path(bedfiles_directory, bf))})
    valids <- which(file.info(out)$size != 0)

    if(length(valids) == 0){
        return(glue('INFO - No file is valid'))
    }

    bedfiles_collapsed <- paste0(bedfiles, collapse=':')
    return(c(bedfiles_directory, bedfiles_collapsed))

}

completed_mtdt <- apply(mt, 1, function(mm){
    ass <- mm[1] ; ctext_lvl <- mm[2]
    dd <- get_dcids(ass, ctext_lvl, assay_dt, ctext)
    ee <- check_bedfiles(dd, opt$bedfiles_directory)
    nctext_lvl <- gsub(' ', '-', ctext_lvl) # e.g. Bone marrow
    ff <- c(ass, nctext_lvl, ee) |> as.data.frame() |> t()
    return(ff)
}) |> as.data.frame() |> t()

colnames(completed_mtdt) <- c('assay', 'context', 'bed_directory', 'bedfiles')
data.table::fwrite(completed_mtdt, '/project2/haky/temi/projects/TFPred-snakemake/metadata/new_dt.txt', sep='\t', row.names=F)




# Rscript # --mt_file /project2/haky/temi/projects/TFPred-snakemake/metadata/mtfile.csv --assay_db /project2/haky/Data/TFXcan/cistrome/raw/human_factor_full_QC.txt --info_db /project2/haky/temi/projects/TFPred-snakemake/info/data_db.txt --bedfiles_directory