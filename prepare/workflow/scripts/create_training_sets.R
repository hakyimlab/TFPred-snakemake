

args <- commandArgs(trailingOnly = TRUE)
TF <- args[1]
tissue <- args[2]
predicted_motif_file <- args[3]
data_dir <- args[4]
output_file <- args[5]
config_file <- args[6]

print(config_file)

valid_chromosomes <- c(paste('chr', 1:22, sep=''), "chrX")

library(yaml)
library(data.table)
library(glue)
library(tidyverse)
library(GenomicRanges)

# if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
# BiocManager::install("GenomicRanges")

print(glue('TF: {TF}, tissue: {tissue}, predicted_motif_file: {predicted_motif_file}, data_dir: {data_dir}, output_file: {output_file}'))

print(getwd())
print(file.exists(predicted_motif_file))

bname <- output_file %>% str_replace('_predictor_regions.txt', '')
# output_file <- glue('{data_dir}/predictor_files/{TF}_{tissue}_predictors.txt')

genome_wide_predicted_motifs <- data.table::fread(predicted_motif_file)
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
peak_files_paths <- list.files(glue('{data_dir}/bed_files/{TF}_{tissue}/'), pattern = '*.bed', full.names = TRUE)
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
cistrome_dt_neg <- dplyr::slice_sample(dt_merged[dt_merged$binding_counts == 0, ], n=nrow(dt_merged) * num_negs)[, 1:5]
cistrome_dr <- rbind(cistrome_dt_pos, cistrome_dt_neg) %>%
    tidyr::unite('region', c(chr, start, end), remove=T) %>% 
    as.data.table()
cistrome_dr <- cistrome_dr[sample(nrow(cistrome_dr)), ]

# write to file
write.table(dt_merged, file=output_file, sep='\t', quote=F, row.names=F)
pfile <- glue('{bname}_predictors.txt')
gfile <- glue('{bname}_ground_truth.txt')
write.table(cistrome_dr[, 1], pfile, col.names=F, quote=F, row.names=F)
write.table(cistrome_dr, gfile, col.names=F, quote=F, row.names=F)

directives <- yaml::yaml.load_file(config_file)
print(directives)
# create tehe enformer prediction config file
enformer_parameters_json <- directives$enformer$prediction_directives
# you may change these as appropriate
enformer_parameters_json[['prediction_data_name']] <- directives$dataset
enformer_parameters_json[['prediction_id']] <- glue('{TF}_{tissue}')
enformer_parameters_json[['date']] <- Sys.Date()
enformer_parameters_json[['interval_list_file']] <- pfile

write(
    jsonlite::toJSON(enformer_parameters_json, na='null', pretty=TRUE, auto_unbox=T),
    file=glue('{enformer_parameters_json$metadata_dir}/enformer_parameters_{directives$dataset}_{TF}_{tissue}.json')
)
# 
#param_file <- glue('{enformer_parameters_json$metadata_dir}/enformer_parameters_{directives$dataset}_{TF}_{tissue}.json')