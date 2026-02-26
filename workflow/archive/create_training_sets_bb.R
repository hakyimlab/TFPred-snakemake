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
    make_option("--train_split", type="double", default=0.8, help='proportion to be used as train'),
    make_option("--num_predictors", type="integer", default=40000, help='proportion to be used as train')
)

opt <- parse_args(OptionParser(option_list=option_list))

seed <- 2023
set.seed(seed)

library(glue)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(plyranges)


# setwd('/project2/haky/temi/projects/TFPred-snakemake')
# opt <- list()
# opt$transcription_factor <- 'CTCF'
# opt$tissue <- 'Retina' 
# opt$predicted_motif_file <- glue('data/homer_files/{opt$transcription_factor}/merged_motif_file.txt')
# opt$bedfiles_directory <- '/project2/haky/Data/TFXcan/cistrome/raw/human_factor' 
# opt$bedlinks_directory <- glue('data/bed_links/{opt$transcription_factor}_{opt$tissue}')
# opt$predictors_file <- 'data/predictor_files/{opt$transcription_factor}_{opt$tissue}.predictors.txt' 
# opt$ground_truth_file <- 'data/predictor_files/{opt$transcription_factor}){opt$tissue}.ground_truth.txt' 
# opt$info_file <- 'data/predictor_files/{opt$transcription_factor}){opt$tissue}.info.txt.gz' 
# opt$cistrome_metadata_file <- '/project2/haky/Data/TFXcan/cistrome/raw/human_factor_full_QC.txt'
# opt$train_split <- 0.8
# opt$num_predictors <- 40000

# # a hacky way to deal with the "none" vs 'None'
# if(opt$tissue == 'none'){
#     opt$tissue <- 'None'
# }

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
genome_wide_predicted_motifs <- data.table::fread(opt$predicted_motif_file) %>% 
    dplyr::select(chr=V2, start=V3, end=V4, strand=V5, score=V6) %>% 
    dplyr::filter(chr %in% valid_chromosomes, score > quantile(.$score, 0.25))


# get threshold
# threshold <- genome_wide_predicted_motifs %>% 
#     dplyr::pull(score) %>% 
#     quantile(0.05)
# genome_wide_predicted_motifs <- genome_wide_predicted_motifs %>% 
# dplyr::filter(score >= threshold)

# turn into GRanges
tf_motifs_granges <- with(genome_wide_predicted_motifs, GRanges(chr, IRanges(start,end), strand, score))
tf_motifs_granges <- tf_motifs_granges[seqnames(tf_motifs_granges) %in% valid_chromosomes]
tf_motifs_granges <- GenomicRanges::reduce(tf_motifs_granges)
tf_motifs_granges <- keepSeqlevels(tf_motifs_granges, paste0("chr", c(1:22, 'X')), pruning.mode="coarse")

# === read in the cistrome metadata file and prepare the bed files
# mtdt <- data.table::fread(opt$cistrome_metadata_file)
# ftissue <- gsub('-', ' ', opt$tissue)
# TF_data <- base::subset(x=mtdt, subset = (Factor == opt$transcription_factor) & (Tissue_type == ftissue))
# if( ! (nrow(TF_data) > 1 & ncol(TF_data) > 1)){
#     print(glue('INFO - No DCids found for {opt$transcription_factor} and {opt$tissue}'))
#     quit(status = 1)
# } else {
#     TF_files <- list.files(glue('{opt$bedfiles_directory}'), pattern=paste0('^', TF_data$DCid, collapse='_*|'), full.names=T)
# }



# out <- sapply(TF_files, function(each_file){
#     output_file <- normalizePath(opt$bedlinks_directory)
#     bname <- basename(each_file)
#     output_file <- file.path(output_file, bname)
#     if(!file.exists(output_file)){
#         cmd <- glue('ln -s {each_file} {output_file}')
#         system(cmd)
#     }
#     return(output_file)
#     #print(base::normalizePath(output_dir))
#     #to_ = glue("{normalizePath(output_dir)}/{basename(each_file)}")
#     #base::file.symlink(from=each_file, to=to_)
# })

TF_files <- list.files(glue('{opt$bedlinks_directory}'), full.names=T)
peak_files_paths <- TF_files[file.info(TF_files)$size != 0]

# use just 1000
if(length(peak_files_paths) > 600){
    print(glue('INFO - There are {length(peak_files_paths)} bedfiles for for {opt$transcription_factor} and {opt$tissue}; using a random 600 instead'))
    peak_files_paths <- sample(peak_files_paths, size=600, replace=FALSE)
}

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
        dplyr::mutate(class = 1, intensity = positive_mt@elementMetadata$intensity, peakOffset = positive_mt@elementMetadata$peakOffset) %>%
        dplyr::filter(intensity > quantile(.$intensity, 0.25))
    
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

ptm <- proc.time()
dt_merged <- pmi_dt_list %>% 
    purrr::reduce(full_join, by = c('chr', 'start', 'end')) %>% 
    base::replace(is.na(.), 0) %>%
    dplyr::mutate(binding_counts = rowSums(dplyr::pick(starts_with('class_'))),
                    mean_intensity = rowMeans(dplyr::pick(starts_with('intensity_'))),
                    binding_class = ifelse(binding_counts > 0, 1, 0),
                    chr = as.character(chr)) %>%
    dplyr::select(chr, start, end, binding_class, binding_counts, mean_intensity)

print(glue('INFO - Time to merge files: {as.vector(proc.time() - ptm)[1]} seconds'))

# ss <- sort(unique(dt_pos$binding_counts))
# tt <- ss/sum(ss)
# names(tt) <- ss
# dt_pos <- dt_pos %>%
#     mutate(proba = tt[as.character(binding_counts)])
# dt_neg <- dt_merged[dt_merged$binding_counts == 0, ]
# dt_neg$proba <- 0
# sample(1:10, size=10, replace=F, prob = runif(10, min=0, max=1))


# filtering === use the num of positives to determine the number of negatives if total is less than 40000
dt_pos <- dt_merged[dt_merged$binding_counts > 0, ]

if(nrow(dt_pos) < (opt$num_predictors/2)){
    #npositives <- ceiling(nrow(dt_merged)/2)
    npositives <- nrow(dt_pos)
    nnegatives <- npositives
} else {
    npositives <- ceiling(opt$num_predictors/2)
    nnegatives <- npositives
}

# I need to select the positive sets so that I have all or most of the highly confident regions
if(nrow(dt_pos) >= ceiling(opt$num_predictors/2)){ 
    csum <- cumsum(sort(dt_pos$binding_counts |> table()))
    wsum <- which(csum <= npositives) |> names() |> as.numeric()
    dtp1 <- dt_pos %>% dplyr::filter(binding_counts %in% wsum)
    if(nrow(dtp1 < npositives)){
        dtp2 <- dt_pos %>% dplyr::filter(!binding_counts %in% wsum)
        remainder <- npositives - nrow(dtp1)
        dtp2 <- dplyr::slice_sample(dtp2, n=remainder)
        dt_pos <- rbind(dtp1, dtp2)
    } else {
        dt_pos <- dtp1
    }
    
}

# create probility distribution for the binding counts
# rg <- range(dt_pos$binding_counts)
# positive_set_threshold <- c(rg[1]:rg[2]) |> quantile(0.25)
# dt_pos <- dt_pos[dt_pos$binding_counts > positive_set_threshold, ]#[, 1:5]
dt_neg <- dt_merged[dt_merged$binding_counts == 0, ]
dt_pos <- dt_pos[with(dt_pos, sample(seq_along(binding_counts), size=npositives, replace=F)), ] 
dt_neg <- dt_neg[with(dt_neg, sample(seq_along(binding_counts), size=nnegatives, replace=F)), ]


data.table::fwrite(dt_merged, glue('{opt$info_file}'), row.names=F, quote=F, col.names=T, compress='gzip', sep='\t')

# dt_merged <- data.table::fread('/project2/haky/temi/projects/TFPred-snakemake/data/predictor_files/AR_Breast_predictors.txt')

print(glue('INFO - There are {nrow(dt_merged)} motifs to be trained and tested on before filtering'))
print(glue("INFO - Distribution of negatives and positives before filtering: {paste0(table(dt_merged$binding_class), collapse=' & ')} respectively"))


### merge and split into train and test

dt <- rbind(dt_pos, dt_neg) %>% 
    as.data.frame() %>%
    tidyr::unite('locus', c(chr, start, end), remove=T) %>% 
    as.data.table() %>%
    dplyr::sample_frac(size=1)

# cistrome_dt_pos <- dt[dt$binding_counts > 0, ] %>%
#     dplyr::slice_sample(n=npositives)

# num_negs <- 1
# cistrome_dt_neg <- dplyr::slice_sample(dt_neg, n=nrow(cistrome_dt_pos) * num_negs)

    
print(glue('INFO - There are {nrow(dt)} motifs to be trained and tested on after filtering'))
print(glue("INFO - Distribution of negatives and positives after filtering: {paste0(table(dt$binding_class), collapse=' & ')} respectively"))

tr_size <- ceiling(nrow(dt) * opt$train_split)
tr_indices <- sample(seq_len(nrow(dt)), tr_size)
dt$split <- NA
dt$split[tr_indices] <- 'train'
dt$split[-tr_indices] <- 'test'
dt <- dt[sample(nrow(dt)), ]

print(glue('INFO - Writing out files...'))
data.table::fwrite(as.data.frame(dt[, 1]), glue('{opt$predictors_file}'), row.names=F, quote=F, col.names=F, sep='\t')
data.table::fwrite(dt, glue('{opt$ground_truth_file}'), row.names=F, quote=F, col.names=T, sep='\t')