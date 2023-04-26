
args <- commandArgs(trailingOnly = TRUE)

TF_tissue <- args[1] # e.g. AR_Breast
predicted_motifs_file <- args[2] # e.g. data/homer_files/AR/merged_motif_file.txt
bedfiles_dir <- args[3] # e.g. data/bed_files/AR_Breast
output_predictors <- args[4] # e.g. data/predictor_files/AR_Breast_predictors.txt
output_ground_truth <- args[5] #e.g. data/predictor_files/AR_Breast_ground_truth.txt
output_enformer_config <- args[6] # e.g. enformer_config_AR_Breast.json
config_file <- args[7]


# TF <- 'FOXA1' 
# tissue <- 'Breast' 
# predicted_motif_file <- 'data/homer_files/FOXA1/merged_motif_file.txt' 
# data_dir <- 'data/' 
# output_file <- 'data/predictor_files/FOXA1_Breast_predictor_regions.txt' 
# config_file <- 'config/pipeline.yaml'

print(config_file)
valid_chromosomes <- c(paste('chr', 1:22, sep=''), "chrX")

library(yaml)
library(data.table)
library(glue)
library(tidyverse)
library(GenomicRanges)

# if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
# BiocManager::install("GenomicRanges")

print(glue('TF_tissue: {TF_tissue}, predicted_motifs_file: {predicted_motifs_file}'))

print(getwd())
print(file.exists(predicted_motifs_file))

#bname <- output_file %>% str_replace('_predictor_regions.txt', '')
# output_file <- glue('{data_dir}/predictor_files/{TF}_{tissue}_predictors.txt')

genome_wide_predicted_motifs <- data.table::fread(predicted_motifs_file)
genome_wide_predicted_motifs <- genome_wide_predicted_motifs %>% 
    dplyr::select(chr=V2, start=V3, end=V4, strand=V5, score=V6) %>% 
    dplyr::filter(chr %in% valid_chromosomes)

# get threshold
threshold <- genome_wide_predicted_motifs %>% 
    dplyr::pull(score) %>% 
    quantile(0.05)
genome_wide_predicted_motifs <- genome_wide_predicted_motifs %>% 
dplyr::filter(score >= threshold)

# turn into GRanges
tf_motifs_granges <- with(genome_wide_predicted_motifs, GRanges(chr, IRanges(start,end), strand, score))
tf_motifs_granges <- tf_motifs_granges[seqnames(tf_motifs_granges) %in% valid_chromosomes]
tf_motifs_granges <- GenomicRanges::reduce(tf_motifs_granges)


# peak files
peak_files_paths <- list.files(bedfiles_dir, pattern = '*.bed', full.names = TRUE)
peak_files_paths <- peak_files_paths[file.info(peak_files_paths)$size != 0]
peak_files_list <- purrr::map(.x=peak_files_paths, .f=data.table::fread, .progress=T)

pmi_dt_list <- purrr::map(.x=peak_files_paths, function(each_file){
    dt <- data.table::fread(each_file) %>%
        distinct(V1, V2, .keep_all=T) %>%
        dplyr::select(chr=V1, start=V2, end=V3) %>% # select the chr, start and end columns
        with(., GRanges(chr, IRanges(start, end), strand='+', score=0))
    
    dt <- dt[seqnames(dt) %in% valid_chromosomes]

    overlaps <- GenomicRanges::findOverlaps(query=dt, subject=tf_motifs_granges, type='any')

    positive_dt <- tf_motifs_granges[subjectHits(overlaps), ] %>% # because I only want the motifs
        as.data.frame() %>%
        dplyr::select(chr=seqnames, start, end) %>%
        dplyr::mutate(class = 1)
    

    negative_dt <- tf_motifs_granges[-subjectHits(overlaps), ] %>% # because I only want the motifs
        as.data.frame() %>%
        dplyr::select(chr=seqnames, start, end) %>%
        dplyr::mutate(class = 0)

    return(rbind(positive_dt, negative_dt) |> as.data.frame())

}, .progress=T)

pmi_dt_list <- lapply(seq_along(pmi_dt_list), function(i){
    colnames(pmi_dt_list[[i]])[4] <- paste('class_', i, sep='')
    return(data.table::as.data.table(pmi_dt_list[[i]]))
})

dt_merged <- pmi_dt_list %>% purrr::reduce(full_join, by = c('chr', 'start', 'end')) 
dt_merged$binding_counts <- rowSums(dt_merged[, -c(1:3)], na.rm=T)
dt_merged$binding_class <- ifelse(dt_merged$binding_counts > 0, 1, 0)
dt_merged <- dt_merged %>%
    dplyr::relocate(c('binding_class', 'binding_counts'), .after=end)

# shuffle the data
set.seed(2023)
dt_merged <- dt_merged[sample(nrow(dt_merged)), ]
dt_merged$chr <- as.character(dt_merged$chr)

# dt_merged <- data.table::fread('/project2/haky/temi/projects/TFPred-snakemake/data/predictor_files/AR_Breast_predictors.txt')

rg <- range(dt_merged$binding_counts)
positive_set_threshold <- c(rg[1]:rg[2]) |> quantile(0.25)
cistrome_dt_pos <- dt_merged[dt_merged$binding_counts > positive_set_threshold, ][, 1:5]

num_negs <- 1
set.seed(2023)
cistrome_dt_neg <- dplyr::slice_sample(dt_merged[dt_merged$binding_counts == 0, ], n=nrow(cistrome_dt_pos) * num_negs)[, 1:5]
cistrome_dr <- rbind(cistrome_dt_pos, cistrome_dt_neg) %>%
    tidyr::unite('region', c(chr, start, end), remove=T) %>% 
    as.data.table()
cistrome_dr <- cistrome_dr[sample(nrow(cistrome_dr)), ]

print(dim(cistrome_dr))

# write to file
#write.table(dt_merged, file=output_file, sep='\t', quote=F, row.names=F)
# pfile <- glue('{output_predictors')
# gfile <- glue('{output_ground_truth}')
write.table(cistrome_dr[, 1], output_predictors, col.names=F, quote=F, row.names=F)
write.table(cistrome_dr, output_ground_truth, col.names=F, quote=F, row.names=F)


# read and write the enformer config file
directives <- yaml::yaml.load_file(config_file)
enformer_parameters_json <- directives$enformer$prediction_directives
# you may change these as appropriate
enformer_parameters_json[['prediction_data_name']] <- directives$dataset
enformer_parameters_json[['prediction_id']] <- glue('{TF_tissue}')
enformer_parameters_json[['date']] <- Sys.Date()
enformer_parameters_json[['interval_list_file']] <- output_predictors

# chANGE the metadata dir
enformer_parameters_json[['metadata_dir']] <- dirname(output_enformer_config)

write(
    jsonlite::toJSON(enformer_parameters_json, na='null', pretty=TRUE, auto_unbox=T),
    file=output_enformer_config
)
# 
#param_file <- glue('{enformer_parameters_json$metadata_dir}/enformer_parameters_{directives$dataset}_{TF}_{tissue}.json')