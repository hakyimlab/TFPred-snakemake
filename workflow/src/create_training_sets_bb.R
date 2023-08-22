# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--transcription_factor", help='A transcription factor e.g. AR'),
    make_option("--tissue", help='A tissue e.g. Breast'),
    make_option("--predicted_motif_file", help="Predicted motif file, particularly from HOMER"),
    make_option("--bedfiles_directory", help='The folder where the cistrome bedfiles are located'),
    make_option("--bedlinks_directory", help='A folder where relevant bedfiles will be linked'),
	make_option("--predictors_file", help='The final predictor file'),
    make_option("--ground_truth_file", help='The final file containing predictors and ground truths'),
    make_option("--info_file", help='The final file containing predictors and extra information'),
    make_option("--cistrome_metadata_file", help='Cistrome metadata file to be read and to find relevant bedfiles'),
    make_option("--train_split", type="double", default=0.8, help='proportion to be used as train')
)

opt <- parse_args(OptionParser(option_list=option_list))


# setwd('./projects/TFPred-snakemake')
# opt <- list()
# opt$transcription_factor <- 'AR'
# opt$tissue <- 'Breast' 
# opt$predicted_motif_file <- 'data/homer_files/AR/merged_motif_file.txt' 
# opt$bedfiles_directory <- '/project2/haky/Data/TFXcan/cistrome/raw/human_factor' 
# opt$bedlinks_directory <- 'data/bed_links/AR_Breast' 
# opt$predictors_file <- 'data/predictor_files/AR_Breast.predictors.txt' 
# opt$ground_truth_file <- 'data/predictor_files/AR_Breast.ground_truth.txt' 
# opt$info_file <- 'data/predictor_files/AR_Breast.info.txt.gz' 
# opt$cistrome_metadata_file <- '/project2/haky/Data/TFXcan/cistrome/raw/human_factor_full_QC.txt'
# opt$train_split <- 0.8

seed <- 2023
set.seed(seed)

library(glue)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(plyranges)


if(!file.exists(opt$cistrome_metadata_file)){
    print(glue('INFO - Cistrome metadata file cannnot be found at {dirname(opt$cistrome_metadata_file)}'))
    quit(status = 1)
} 
if(!file.exists(opt$predicted_motif_file)){
    print(glue('INFO - HOMER predicted motifs file cannnot be found at {dirname(opt$predicted_motif_file)}'))
    quit(status = 1)
} 

valid_chromosomes <- c(paste('chr', 1:22, sep=''), "chrX")

# === prepare the predicted motifs file 
genome_wide_predicted_motifs <- data.table::fread(opt$predicted_motif_file)
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
tf_motifs_granges <- keepSeqlevels(tf_motifs_granges, paste0("chr", c(1:22, 'X')), pruning.mode="coarse")

# === read in the cistrome metadata file and prepare the bed files
mtdt <- data.table::fread(opt$cistrome_metadata_file)
TF_data <- base::subset(x=mtdt, subset = (Factor == opt$transcription_factor) & (Tissue_type == opt$tissue))
if( ! (nrow(TF_data) > 1 & ncol(TF_data) > 1)){
    print(glue('INFO - No DCids found for {opt$transcription_factor} and {opt$tissue}'))
    quit(status = 1)
} else {
    TF_files <- list.files(glue('{opt$bedfiles_directory}'), pattern=paste0('^', TF_data$DCid, collapse='_*|'), full.names=T)
}

out <- sapply(TF_files, function(each_file){
    output_file <- normalizePath(opt$bedlinks_directory)
    bname <- basename(each_file)
    output_file <- file.path(output_file, bname)
    if(!file.exists(output_file)){
        cmd <- glue('ln -s {each_file} {output_file}')
        system(cmd)
    }
    return(output_file)
    #print(base::normalizePath(output_dir))
    #to_ = glue("{normalizePath(output_dir)}/{basename(each_file)}")
    #base::file.symlink(from=each_file, to=to_)
})
peak_files_paths <- out[file.info(out)$size != 0]

pmi_dt_list <- purrr::map(.x=seq_along(peak_files_paths), function(each_file_index){
    each_file <- peak_files_paths[each_file_index]
    dt <- data.table::fread(each_file) %>%
        dplyr::distinct(V1, V2, .keep_all=T) %>%
        dplyr::select(chr=V1, start=V2, end=V3, intensity=V7, peakOffset=V10) %>% # select the chr, start and end columns
        with(., GRanges(chr, IRanges(start, end), strand='+', intensity=intensity, peakOffset=peakOffset))
    
    dt <- dt[seqnames(dt) %in% valid_chromosomes]

    # find overlaps with the predicted motifs
    overlaps <- GenomicRanges::findOverlaps(query=dt, subject=tf_motifs_granges, type='any')

    # those with any overlaps are the positive ones
    positive_mt <- dt[queryHits(overlaps), ]
    positive_dt <- tf_motifs_granges[subjectHits(overlaps), ] %>% # because I only want the motifs
        as.data.frame() %>%
        dplyr::select(chr=seqnames, start, end) %>%
        dplyr::mutate(class = 1, intensity = positive_mt@elementMetadata$intensity, peakOffset = positive_mt@elementMetadata$peakOffset)
    
    # those with no overlaps are negative
    #negative_mt <- dt[-queryHits(overlaps), ]
    negative_dt <- tf_motifs_granges[-subjectHits(overlaps), ] %>% # because I only want the motifs
        as.data.frame() %>%
        dplyr::select(chr=seqnames, start, end) %>%
        dplyr::mutate(class = 0, intensity = 0, peakOffset = 0) %>% 
        dplyr::slice_sample(n=nrow(positive_dt))

    cname <- paste('class_', each_file_index, sep='')
    iname <- paste('intensity_', each_file_index, sep='')
    data <- rbind(positive_dt, negative_dt) %>% 
        data.table::setDT() %>% 
        dplyr::select(-peakOffset) %>%
        dplyr::rename(!!quo_name(cname) := class, !!quo_name(iname) := intensity)
    return(data)

}, .progress=T)

dt_merged <- pmi_dt_list %>% 
    purrr::reduce(full_join, by = c('chr', 'start', 'end')) %>% 
    base::replace(is.na(.), 0) %>%
    dplyr::mutate(binding_counts = rowSums(dplyr::pick(starts_with('class_'))),
                    mean_intensity = rowMeans(dplyr::pick(starts_with('intensity_'))),
                    binding_class = ifelse(binding_counts > 0, 1, 0),
                    chr = as.character(chr)) %>%
    dplyr::select(chr, start, end, binding_class, binding_counts, mean_intensity)

# dt_merged <- data.table::fread('/project2/haky/temi/projects/TFPred-snakemake/data/predictor_files/AR_Breast_predictors.txt')

rg <- range(dt_merged$binding_counts)
positive_set_threshold <- c(rg[1]:rg[2]) |> quantile(0.25)
cistrome_dt_pos <- dt_merged[dt_merged$binding_counts > positive_set_threshold, ]

num_negs <- 1
cistrome_dt_neg <- dplyr::slice_sample(dt_merged[dt_merged$binding_counts == 0, ], n=nrow(cistrome_dt_pos) * num_negs)

cistrome_dr <- rbind(cistrome_dt_pos, cistrome_dt_neg) %>%
    tidyr::unite('locus', c(chr, start, end), remove=T) %>% 
    as.data.table()

tr_size <- ceiling(nrow(cistrome_dr) * opt$train_split)
tr_indices <- sample(1:nrow(cistrome_dr), tr_size)
cistrome_dr$split <- NA
cistrome_dr$split[tr_indices] <- 'train'
cistrome_dr$split[-tr_indices] <- 'test'
cistrome_dr <- cistrome_dr[sample(nrow(cistrome_dr)), ]

print(glue('INFO - Writing out files...'))
data.table::fwrite(as.data.frame(cistrome_dr[, 1]), glue('{opt$predictors_file}'), row.names=F, quote=F, col.names=F, sep='\t')
data.table::fwrite(cistrome_dr, glue('{opt$ground_truth_file}'), row.names=F, quote=F, col.names=T, sep='\t')
data.table::fwrite(cistrome_dr, glue('{opt$info_file}'), row.names=F, quote=F, col.names=T, compress='gzip', sep='\t')