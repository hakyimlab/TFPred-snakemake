# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--transcription_factor", help='A transcription factor e.g. AR'),
    make_option("--tissue", help='A tissue e.g. Breast'),
    make_option("--predicted_motif_file", help="Predicted motif file, particularly from HOMER"),
    make_option("--sorted_bedfiles_directory", help='The folder where the cistrome bedfiles are located'),
    make_option("--bedlinks_directory", help='A folder where relevant bedfiles will be linked'),
	make_option("--predictors_file", help='The final predictor file'),
    make_option("--ground_truth_file", help='The final file containing predictors and ground truths'),
    make_option("--info_file", help='The final file containing predictors and extra information'),
    make_option("--summary_file", type = "character"),
    make_option("--train_split", type="double", default=0.8, help='proportion to be used as train'),
    make_option("--num_predictors", type="integer", default=40000, help='proportion to be used as train'),
    make_option("--test_chromosomes", type="character", default="chr9,chr22", help='should you randomly split or train by chromosome? num_predictors is ignored if training by chromosome'),
    make_option("--peaks_files", type="character", default=NULL, help='how many bed samples should be used: will use all if NULL, will use the top n by FRiP otherwise'),
    make_option("--sorted_chrom_sizes", type="character", default="/project/haky/users/temi/projects/TFPred-snakemake/info/hg38.chrom.sizes.sorted", help='The sorted chromosome sizes file'),
    make_option("--peaks_counts_threshold", type="integer", default=100, help='...')
)

opt <- parse_args(OptionParser(option_list=option_list))

seed <- 2023
set.seed(seed)

library(glue)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(plyranges)


#  /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/bin/Rscript workflow/src/create_training_sets_bedtools.R --transcription_factor AR --tissue MammaryGland --predicted_motif_file data/homer_instances/AR/merged_motif_file.txt --sorted_bedfiles_directory data/cistrome_2024-06-08/sortedbeds/AR_MammaryGland --bedlinks_directory /project/haky/data/TFXcan/cistrome/raw/human_factor' --predictors_file data/cistrome_2024-06-08/predictor_files/AR_MammaryGland.predictors.txt --ground_truth_file data/cistrome_2024-06-08/predictor_files/AR_MammaryGland.ground_truth.txt --info_file data/cistrome_2024-06-08/predictor_files/AR_MammaryGland.info.txt.gz --peaks_files 36845_sort_peaks.narrowPeak.bed,57275_sort_peaks.narrowPeak.bed,57276_sort_peaks.narrowPeak.bed,57277_sort_peaks.narrowPeak.bed,57278_sort_peaks.narrowPeak.bed,57279_sort_peaks.narrowPeak.bed,57280_sort_peaks.narrowPeak.bed --train_by_chromosome chr9,chr22 --summary_file data/cistrome_2024-06-08/predictor_files/AR_MammaryGland.summary.txt; sleep 5


# /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/bin/Rscript workflow/src/create_training_sets_bedtools.R --transcription_factor MAFF --tissue Colon --predicted_motif_file data/homer_instances/MAFF/merged_motif_file.txt --sorted_bedfiles_directory data/ENPACT_275_2024-06-08/sortedbeds/MAFF_Colon --bedlinks_directory /project/haky/data/TFXcan/cistrome/raw/human_factor --predictors_file data/ENPACT_275_2024-06-08/predictor_files/MAFF_Colon.predictors.txt --ground_truth_file data/ENPACT_275_2024-06-08/predictor_files/MAFF_Colon.ground_truth.txt --info_file data/ENPACT_275_2024-06-08/predictor_files/MAFF_Colon.info.txt.gz --peaks_files 42840_sort_peaks.narrowPeak.bed --test_chromosomes chr9,chr22 --summary_file data/ENPACT_275_2024-06-08/predictor_files/MAFF_Colon.summary.txt;

# setwd('/beagle3/haky/users/temi/projects/TFPred-snakemake')
# opt <- list()
# opt$transcription_factor <- 'MAFF'
# opt$tissue <- 'Colon' 
# opt$predicted_motif_file <- glue('data/homer_instances/{opt$transcription_factor}/merged_motif_file.txt')
# opt$sorted_bedfiles_directory <- glue('data/ENPACT_275_2024-06-08/sortedbeds/{opt$transcription_factor}_{opt$tissue}')
# opt$peaks_files <- '42840_sort_peaks.narrowPeak.bed'
# opt$bedlinks_directory <- '/project/haky/data/TFXcan/cistrome/raw/human_factor'
# opt$predictors_file <- glue('data/predictor_files/{opt$transcription_factor}_{opt$tissue}.predictors.txt')
# opt$ground_truth_file <- glue('data/predictor_files/{opt$transcription_factor}_{opt$tissue}.ground_truth.txt')
# opt$info_file <- glue('data/predictor_files/{opt$transcription_factor}_{opt$tissue}.info.txt.gz')
# opt$summary_file <- glue('data/predictor_files/{opt$transcription_factor}_{opt$tissue}.summary.txt')
# opt$train_split <- 0.8
# opt$num_predictors <- 40000
# opt$test_chromosomes <- "chr9,chr22"



# # a hacky way to deal with the "none" vs 'None'
# if(opt$tissue == 'none'){
#     opt$tissue <- 'None'
# }

if(!dir.exists(opt$sorted_bedfiles_directory)){
    dir.create(opt$sorted_bedfiles_directory, recursive = T)
}


if(!file.exists(opt$predicted_motif_file)){
    print(glue('INFO - HOMER predicted motifs file cannnot be found at {dirname(opt$predicted_motif_file)}'))
    quit(status = 1)
} 

valid_chromosomes <- paste0('chr', 1:22) #, "chrX")

# === prepare the predicted motifs file 
genome_wide_predicted_motifs <- data.table::fread(opt$predicted_motif_file) %>% 
    dplyr::select(chr=V2, start=V3, end=V4, strand=V5, score=V6) %>% 
    dplyr::filter(chr %in% valid_chromosomes, score > quantile(.$score, 0.25))

# turn into GRanges
tf_motifs_granges <- with(genome_wide_predicted_motifs, GRanges(chr, IRanges(start,end), strand, score))
tf_motifs_granges <- tf_motifs_granges[seqnames(tf_motifs_granges) %in% valid_chromosomes]
tf_motifs_granges <- GenomicRanges::reduce(tf_motifs_granges)
tf_motifs_granges <- keepSeqlevels(tf_motifs_granges, paste0("chr", c(1:22)), pruning.mode="coarse")



# check that all of these are available in the db
# TF_files <- list.files(glue('{opt$bedlinks_directory}'), full.names=T)

if(!is.null(opt$peaks_files)){
    peaks_files <- base::strsplit(opt$peaks_files, ',')[[1]] |> sort()
} 

peaks_files <- sapply(peaks_files, function(each_file){
    et <- file.path(opt$bedlinks_directory, each_file)
    if(file.exists(et)){
        return(et)
    }
})

# if(!is.null(opt$dcids)){
#     TF_files <- TF_files[grepl(paste0(dcid_samples, collapse='|'), TF_files)]
# }

peak_files_paths <- peaks_files[file.info(peaks_files)$size != 0]

if(length(peak_files_paths) == 0){
    print(glue('INFO - There are no valid peaks data for {opt$transcription_factor} in {opt$tissue}'))
    quit(status = 1)
}

# use just 1000
#if(length(peak_files_paths) > 600){
print(glue('INFO - There are {length(peak_files_paths)} bedfiles for {opt$transcription_factor} and {opt$tissue}'))
#peak_files_paths <- sample(peak_files_paths, size=600, replace=FALSE)
#}
set.seed(seed)
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

    # creating bedfiles to merge with bedtools
    mf <- positive_dt %>%
        dplyr::select(chr, start, end) %>%
        dplyr::arrange(factor(chr, levels = valid_chromosomes), start)
    data.table::fwrite(mf, file = glue('{opt$sorted_bedfiles_directory}/peaks_{each_file_index}.sorted.bed'), row.names=F, col.names = F, quote=F, sep='\t')
    return(data)

}, .progress=T)

query_bed <- glue('{opt$sorted_bedfiles_directory}/{opt$transcription_factor}_{opt$tissue}.ALL.query.bed')

## get a union of all peaks
set.seed(seed)
purrr::map(pmi_dt_list, function(each_dt){
    each_dt %>%
        filter(.[[4]] == 1) %>%
        dplyr::select(chr, start, end) 
}, .progress = T) %>%
    do.call('rbind', .) %>%
    dplyr::distinct(chr, start, end, .keep_all = FALSE) %>%
    dplyr::filter(chr %in% valid_chromosomes) %>%
    dplyr::arrange(factor(chr, levels = valid_chromosomes), start) %>%
    data.table::fwrite(file = query_bed, row.names=F, col.names = F, quote=F, sep='\t')
#
intersect_bed <- glue("{opt$sorted_bedfiles_directory}/{opt$transcription_factor}_{opt$tissue}.intersect.bed")

# if(!is.null(opt$dcids)){
#     pfiles <- list.files(opt$sorted_bedfiles_directory, pattern = '\\d.sorted.bed', full.names = T)
#     pfiles <- pfiles[grepl(paste0(dcid_samples, collapse = '|'), pfiles)]
# } else {
pfiles <- list.files(opt$sorted_bedfiles_directory, pattern = '^peaks_\\d.sorted.bed', full.names = T)

#pfiles <- list.files(opt$sorted_bedfiles_directory, pattern = '\\d.sorted.bed', full.names = T)
cmd <- glue("bedtools intersect -c -a {query_bed} -g {opt$sorted_chrom_sizes} -sorted -b {paste(pfiles, collapse = ' ')} > {intersect_bed}")
ptm <- proc.time()
system(cmd)
message(glue('INFO - Time to merge files: {as.vector(proc.time() - ptm)[1]} seconds'))


if(!file.exists(intersect_bed) || file.info(intersect_bed)$size <= 0){
    sapply(list(opt$ground_truth_file, opt$info_file,opt$predictors_file, opt$summary_file), file.create)
    message(glue("WARNING - There are no peaks present after bedtools intersect."))
} else if(file.exists(intersect_bed) && file.info(intersect_bed)$size != 0) {

    dt_pos <- data.table::fread(intersect_bed) %>%
        setNames(c('chr', 'start', 'end', 'binding_counts')) %>%
        dplyr::filter(binding_counts > 0) %>%
        dplyr::mutate(binding_class = 1)

    if(nrow(dt_pos) < opt$peaks_counts_threshold){
        message(glue("WARNING - There are less than {opt$peaks_counts_threshold} peaks after intersecting with motif instances. No model will be trained for {opt$transcription_factor}_{opt$tissue}. The peaks that are intersected will still be written out, but there will be no {opt$transcription_factor}_{opt$tissue}.predictors.txt file"))
    }

    # select negatives and merge
    dt_neg <- purrr::map(pmi_dt_list, function(each_dt){
        each_dt %>%
            filter(.[[4]] == 0) %>%
            dplyr::select(chr, start, end) 
    }) %>%
        do.call('rbind', .) %>%
        dplyr::filter(chr %in% valid_chromosomes) %>%
        dplyr::distinct(chr, start, end, .keep_all = FALSE) %>%
        dplyr::mutate(binding_class = 0, binding_counts = 0)

    dt_merged <- dplyr::bind_rows(dt_pos, dt_neg)
    data.table::fwrite(dt_merged, glue('{opt$info_file}'), row.names=F, quote=F, col.names=T, compress='gzip', sep='\t')

    test_chromosomes <- base::strsplit(opt$test_chromosomes, ',')[[1]]

    if(!any(test_chromosomes %in% dt_pos$chr) || !any(test_chromosomes %in% dt_neg$chr)){
        message(glue("WARNING - None of the test chromosomes are present in the positive set. Now trying for other chromosomes"))
        test_chromosomes <- c('chr4', 'chr3')
    }

    train_dt_pos <- dt_pos %>%
        dplyr::filter(!chr %in% test_chromosomes) %>%
        dplyr::arrange(desc(binding_counts)) 

    test_dt_pos <- dt_pos %>%
        dplyr::filter(chr %in% test_chromosomes) %>%
        dplyr::arrange(desc(binding_counts))

    # select training data ===============
    if(nrow(train_dt_pos) < (opt$num_predictors/2)){
        npositives <- nrow(train_dt_pos)
        nnegatives <- npositives
    } else {
        npositives <- ceiling(opt$num_predictors/2)
        nnegatives <- npositives
    }

    # select test data
    if(nrow(test_dt_pos) < 1000){
        tpositives <- nrow(test_dt_pos)
        tnegatives <- tpositives
    } else {
        tpositives <- 1000
        tnegatives <- tpositives
    }

    ### training by splitting into chromosomes ###
    trp_data <- train_dt_pos %>%
        dplyr::slice_head(n=npositives)

    tep_data <- test_dt_pos %>%
        dplyr::slice_head(n=tpositives)

    trn_data <- dt_neg %>%
        dplyr::filter(!chr %in% test_chromosomes) %>%
        dplyr::slice_sample(n=nnegatives)

    ten_data <- dt_neg %>%
        dplyr::filter(chr %in% test_chromosomes) %>%
        dplyr::slice_sample(n=tnegatives)

    set.seed(seed)
    train_data <- rbind(trp_data, trn_data) %>% dplyr::sample_frac(size = 1) %>% dplyr::mutate(split = 'train')
    set.seed(seed)
    test_data <- rbind(tep_data, ten_data) %>% dplyr::sample_frac(size = 1) %>% dplyr::mutate(split = 'test')

    set.seed(seed)
    dt <- rbind(train_data, test_data) %>%
        as.data.frame() %>%
        tidyr::unite('locus', c(chr, start, end), remove=T) %>% 
        as.data.table() %>%
        dplyr::sample_frac(size=1)

    print(glue('INFO - There are {nrow(dt)} motifs to be trained and tested on after filtering'))
    print(glue("INFO - Distribution of negatives and positives after filtering: {paste0(table(dt$binding_class), collapse=' & ')} respectively"))


    ###
    print(glue('INFO - Writing out files...'))
    if(nrow(dt_pos) >= opt$peaks_counts_threshold){
        data.table::fwrite(as.data.frame(dt[, 1]), glue('{opt$predictors_file}'), row.names=F, quote=F, col.names=F, sep='\t')
    } else {
        file.create(opt$predictors_file)
    }
    data.table::fwrite(dt, glue('{opt$ground_truth_file}'), row.names=F, quote=F, col.names=T, sep='\t')


    #### statistics
    train_test_split <- dt$split %>% table()
    train_test_class <- dt$binding_class %>% table()

    train_class <- train_data %>%
        dplyr::group_by(binding_class, chr) %>%
        dplyr::summarise(n_loci = n()) %>%
        tidyr::pivot_wider(id_cols = binding_class, names_from = chr, values_from = n_loci)

    test_class <- test_data %>%
        dplyr::group_by(binding_class, chr) %>%
        dplyr::summarise(n_loci = n()) %>%
        tidyr::pivot_wider(id_cols = binding_class, names_from = chr, values_from = n_loci)


    train_counts <- train_data %>%
        dplyr::group_by(binding_counts, chr) %>%
        dplyr::summarise(n_loci = n()) %>%
        tidyr::pivot_wider(id_cols = binding_counts, names_from = chr, values_from = n_loci)

    test_counts <- test_data %>%
        dplyr::group_by(binding_counts, chr) %>%
        dplyr::summarise(n_loci = n()) %>%
        tidyr::pivot_wider(id_cols = binding_counts, names_from = chr, values_from = n_loci)

    # generate summary statistics ======
    if(!file.exists(opt$summary_file)){
        if(!file.exists(dirname(opt$summary_file))){
            if(!dir.exists(dirname(opt$summary_file))){
                dir.create(dirname(opt$summary_file), recursive = TRUE)
            }
        }

        diagfile <- file(opt$summary_file, open = "a")
        cat("## Training binding sites distribution", file = diagfile, sep = '\n')
        cat("#### train-test split (split) binding sites distribution", file = diagfile, sep = '\n')
        capture.output(train_test_split, file = diagfile, sep = '\n')
        cat("#### train-test split (0s, 1s) binding sites distribution", file = diagfile, sep = '\n')
        capture.output(train_test_class, file = diagfile, sep = '\n')
        cat("#### training binding class distribution", file = diagfile, sep = '\n')
        capture.output(train_class, file = diagfile, sep = '\n')
        cat("#### test binding class distribution", file = diagfile, sep = '\n')
        capture.output(test_class, file = diagfile, sep = '\n')
        cat("#### training binding counts distribution", file = diagfile, sep = '\n')
        capture.output(train_counts, file = diagfile, sep = '\n')
        cat("#### test binding counts distribution", file = diagfile, sep = '\n')
        capture.output(test_counts, file = diagfile, sep = '\n')

    } else {
        diagfile <- NULL
    }

    if(!is.null(opt$diagnostics_file)){
        close(diagfile)
    }
}









# ts <- data.frame(a=c(1:10, 1:10), b = c(LETTERS[1:10], LETTERS[1:10])) %>%
#     dplyr::distinct(a, .keep_all = F)

# ts <- data.table::fread('/project2/haky/temi/projects/TFPred-snakemake/info/hg38.chrom.sizes') %>%
#     dplyr::filter(V1 %in% valid_chromosomes) %>%
#     dplyr::arrange(factor(V1, levels = valid_chromosomes)) %>%
#     data.table::fwrite(file = '/project2/haky/temi/projects/TFPred-snakemake/info/hg38.chrom.sizes.sorted', row.names=F, col.names=F, sep='\t', quote=F)



# if(opt$train_by_chromosome == TRUE){
#     print(glue('INFO - Splitting by chromosome'))
#     csum <- dt %>% 
#     tidyr::separate(locus, into=c('chrom', 'start', 'end'), sep='_') %>%
#     dplyr::pull(chrom) %>%
#     table() %>%
#     sort() 

#     for(i in seq_along(csum)){
#         if(csum[i] <= 100){
#             tpick <- 0
#             next
#         } else {
#             tpick <- i
#             tchrom <- names(csum)[tpick]
#             break
#         }
#     }

#     # if all are less than 100, pick the maximum
#     if(tpick == 0){
#         tpick <- which.max(csum)
#         tchrom <- names(csum)[tpick]
#     }

#     print(glue('INFO - {tchrom} will be used to test'))

#     dt <- dt %>% 
#         tidyr::separate(locus, into=c('chrom', 'start', 'end'), sep='_') %>%
#         dplyr::mutate(split = dplyr::case_when(
#             chrom == tchrom ~ 'test',
#             TRUE ~ 'train'
#         )) %>%
#         tidyr::unite('locus', c(chrom, start, end), remove=T)

# } else {
#     tr_size <- ceiling(nrow(dt) * opt$train_split)
#     tr_indices <- sample(seq_len(nrow(dt)), tr_size)
#     dt$split <- NA
#     dt$split[tr_indices] <- 'train'
#     dt$split[-tr_indices] <- 'test'
#     dt <- dt[sample(nrow(dt)), ]
# }


# dt %>% 
#     tidyr::separate(locus, into=c('chrom', 'start', 'end'), sep='_') %>%
#     dplyr::group_by(chrom, binding_class) %>%
#     dplyr::summarise(cnt = n()) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(prop = (cnt/sum(cnt))*100)

# csum <- cumsum(sort(dt_pos$binding_counts |> table()))
#     wsum <- which(csum <= npositives) |> names() |> as.numeric()
#     dtp1 <- dt_pos %>% dplyr::filter(binding_counts %in% wsum)
    





# I need to select the positive sets so that I have all or most of the highly confident regions
# if(nrow(dt_pos) >= ceiling(opt$num_predictors/2)){ 
#     csum <- cumsum(sort(dt_pos$binding_counts |> table()))
#     wsum <- which(csum <= npositives) |> names() |> as.numeric()
#     dtp1 <- dt_pos %>% dplyr::filter(binding_counts %in% wsum)
#     if(nrow(dtp1 < npositives)){
#         dtp2 <- dt_pos %>% dplyr::filter(!binding_counts %in% wsum)
#         remainder <- npositives - nrow(dtp1)
#         dtp2 <- dplyr::slice_sample(dtp2, n=remainder)
#         dt_pos <- rbind(dtp1, dtp2)
#     } else {
#         dt_pos <- dtp1
#     }
    
# }

# dt_pos <- dt_pos[with(dt_pos, sample(seq_along(binding_counts), size=npositives, replace=F)), ] 
# dt_neg <- dt_neg[with(dt_neg, sample(seq_along(binding_counts), size=nnegatives, replace=F)), ]

# print(glue('INFO - There are {nrow(dt_merged)} motifs to be trained and tested on before filtering'))
# print(glue("INFO - Distribution of negatives and positives before filtering: {paste0(table(dt_merged$binding_class), collapse=' & ')} respectively"))