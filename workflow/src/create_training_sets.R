
# Author: Temi
# Date: Thursday July 27 2023
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
    make_option("--cistrome_metadata_file", help='Cistrome metadata file to be read and to find relevant bedfiles')
)

opt <- parse_args(OptionParser(option_list=option_list))

seed <- 2023

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

peaks_dt <- purrr::map(.x=names(peak_files_paths), .f=function(bfile){
    bt <- data.table::fread(bfile) %>% dplyr::select(chr=V1, start=V2, end=V3, intensity=V7) %>%
        dplyr::filter(chr %in% paste0("chr", c(1:22, 'X'))) %>% dplyr::distinct()
    return(bt)
}, .progress=TRUE)


reduce_combine <- function(x, y){

    if(!'count' %in% colnames(x)){
        x$count <- 0
    }
    if(!'count' %in% colnames(y)){
        y$count <- 0
    }
    
    tryCatch({
        xdt <- with(x, GRanges(chr, IRanges(start, end), strand = '*', score=0)) #%>% keepSeqlevels(., paste0("chr", c(1:22, 'X')), pruning.mode="error")
        ydt <- with(y, GRanges(chr, IRanges(start, end), strand = '*', score=0)) #%>% keepSeqlevels(., paste0("chr", c(1:22, 'X')), pruning.mode="error")
    }, error = function(e){
        message(paste("Error - :"))
        print(x$chr |> unique()) ; print(y$chr |> unique())
    })
    
    overlaps <- GenomicRanges::findOverlaps(query=xdt, subject=ydt, type='any')

    xover_found <- x[queryHits(overlaps), ] 
    xover_nfound <- rbind(x[-queryHits(overlaps), ], y[-subjectHits(overlaps), ]) %>% as.data.frame()
    xover_found$count <- xover_found$count + 1
    xover <- as.data.frame(rbind(xover_found, xover_nfound))

    #print(head(xover_nfound))

    yover_found <- y[subjectHits(overlaps), ] 
    yover_nfound <- rbind(x[-queryHits(overlaps), ], y[-subjectHits(overlaps), ]) %>% as.data.frame()
    yover_found$count <- yover_found$count + 1
    yover <- as.data.frame(rbind(yover_found, yover_nfound))

    newStart <- cbind(xover$start, yover$start) %>% apply(., 1, min) %>% as.numeric()
    newEnd <- cbind(xover$end, yover$end) %>% apply(., 1, max) %>% as.numeric()
    newIntensity <- cbind(xover$intensity, yover$intensity) %>% rowMeans()
    newCount <- cbind(xover$count, yover$count) %>% apply(., 1, max) %>% as.numeric()
    newdt <- data.frame(chr = xover$chr, start=newStart, end=newEnd, intensity = newIntensity, count = newCount) %>% distinct()

    #print(head(newdt))
    #print(dim(newdt))
    return(newdt)
}

print(glue('INFO - Merging files...'))
dt_merged <- peaks_dt %>% purrr::reduce(reduce_combine) %>% dplyr::group_by(chr, start, end, count) %>% dplyr::summarise(intensity = mean(intensity))
dt_granges <- with(dt_merged, GRanges(chr, IRanges(start, end), strand = '*', count=count, intensity=intensity))

# use plyranges::reduce_ranges to reduce and aggregate
dt_granges <- plyranges::reduce_ranges(dt_granges, count = sum(count), intensity = sum(intensity))
motif_overlaps <- GenomicRanges::findOverlaps(query = dt_granges, subject = tf_motifs_granges, type='any')
dt_minus_motifs <- cbind(dt_merged[queryHits(motif_overlaps), ], n_motifs=0, motif=NA) %>% as_tibble()
dt_plus_motifs <- cbind(dt_merged[queryHits(motif_overlaps), ], as.data.frame(tf_motifs_granges[subjectHits(motif_overlaps), ]) %>% 
    dplyr::select(chr_motif=seqnames, start_motif = start, end_motif = end)) %>% 
    dplyr::mutate(motif = paste(chr_motif, start_motif, end_motif, sep='_')) %>% 
    as_tibble() %>% 
    dplyr::group_by(chr, start, end, count, intensity) %>% 
    dplyr::summarise(n_motifs = n(), motif = paste(motif, collapse = ','))

set.seed(seed)
dt_train <- rbind(dt_plus_motifs, dt_minus_motifs) 
dt_train <- dt_train[sample(nrow(dt_train)), ] %>% 
    dplyr::mutate(peakActivityScore = sqrt(intensity * ((count + 0.1) * (n_motifs + 0.1))))
    # dplyr::rename(peakActivityScore = intensity) %>%
    # dplyr::mutate(peakActivityScore = sqrt(intensity * ((count + 0.1) * (n_motifs + 0.1))))
    #dplyr::mutate(peakActivityScore = (intensity * log10((count + 1) * (n_motifs + 1))))

mw <- width(tf_motifs_granges)[1]
all_motifs <- apply(dt_train, 1, function(each_row){
    #print(each_row[['start']])
    if(is.na(each_row[['motif']])){
        sa <- as.numeric(each_row[['start']])
        ea <- as.numeric(each_row[['end']]) 
        a <- ceiling((sa + ea)/2)
        aa <- a - ceiling(mw/2)
        bb <- (a + ceiling(mw/2)) - 1
        return(paste(each_row[['chr']], '_', aa, '_', bb, sep=''))
    } else {
        return(each_row[['motif']])
    }
})

gt <- cbind(motif = all_motifs, peakActivityScore = dt_train$peakActivityScore) %>% as.data.frame() %>% tidyr::separate_longer_delim(motif, delim = ',') %>% dplyr::distinct(motif, .keep_all=TRUE) %>% dplyr::mutate(peakActivityScore = as.numeric(peakActivityScore))

# === write out 
print(glue('INFO - Writing out files...'))
data.table::fwrite(as.data.frame(gt[, 1]), glue('{opt$predictors_file}'), row.names=F, quote=F, col.names=F)
data.table::fwrite(gt, glue('{opt$ground_truth_file}'), row.names=F, quote=F, col.names=F)
data.table::fwrite(dt_train, glue('{opt$info_file}'), row.names=F, quote=F, col.names=T, compress='gzip')


# /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/bin/Rscript prepare/workflow/scripts/create_training_sets.R --transcription_factor AR --tissue Breast --predicted_motif_file data/homer_files/AR/merged_motif_file.txt --bedfiles_directory /project2/haky/Data/TFXcan/cistrome/raw/human_factor --bedlinks_directory data/bed_links/AR_Breast --predictors_file data/predictor_files/AR_Breast.predictors.txt --ground_truth_file data/predictor_files/AR_Breast.ground_truth.txt --info_file data/predictor_files/AR_Breast.info.txt.gz --cistrome_metadata_file /project2/haky/Data/TFXcan/cistrome/raw/human_factor_full_QC.txt

# intensity <- runif(20, 0, 25)
# motifs <- sample(0:5, 20, replace=F)