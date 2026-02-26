# Author: Temi
# Date: Fri Jan 2 2026
# Description: script to create predictors, ground truth and info files on ChIP-seq data only
# Usage: Rscript create_chippeaks_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--transcription_factor", help='A transcription factor e.g. AR'),
    make_option("--tissue", help='A tissue e.g. Breast'),
    make_option("--output_directory", help='The folder where the cistrome bedfiles are located'),
    make_option("--peaks_directory", help='A folder where relevant bedfiles will be linked'),
    make_option("--peaks_files", type="character", default=NULL, help='the peak files to be used, separated by a comma'),
    make_option("--sorted_chrom_sizes", type="character", default=NULL, help='The sorted chromosome sizes file'),
    make_option("--peaks_counts_threshold", type="integer", default=100, help='...'),
    make_option("--test_chromosomes", type="character", default="chr9,chr22", help='should you randomly split or train by chromosome? num_predictors is ignored if training by chromosome'),
    make_option("--train_split", type="double", default=0.8, help='proportion to be used as train'),
    make_option("--num_train_observations", type="integer", default=20000, help='The number of sites to select to train for each group (bound or unbound)'),
    make_option("--num_test_observations", type="integer", default=1000, help='The number of sites to select to test for each group (bound or unbound)'),
    make_option("--predictors_file", help='The final predictor file'),
    make_option("--ground_truth_file", help='The final file containing predictors and ground truths'),
    make_option("--info_file", help='The final file containing predictors and extra information'),
    make_option("--summary_file", type = "A file that summarizes the data collected for training or testing")
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

seed <- 2026
set.seed(seed)

library(glue)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(plyranges)

valid_chromosomes <- paste0('chr', 1:22)

# first compile all peaks across the genome.
if(!is.null(opt$peaks_files)){
    peaks_files <- base::strsplit(opt$peaks_files, ',')[[1]] |> sort()
    peak_files_paths <- sapply(peaks_files, function(each_file){
        et <- file.path(opt$peaks_directory, each_file)
        if(file.exists(et)){
            return(et)
        }
    })

    print(peak_files_paths)
    #[file.info(peaks_files)$size != 0]
} 

# check that the peaks files are not empty; remove empty files
if(length(peak_files_paths) == 0){
    print(glue('INFO - There are no valid peaks data for {opt$transcription_factor} in {opt$tissue}'))
    quit(status = 0)
}

print(glue('INFO - There are {length(peak_files_paths)} bedfiles for {opt$transcription_factor} and {opt$tissue}'))

if(!dir.exists(opt$output_directory)){
    dir.create(opt$output_directory, recursive = TRUE)
}

# I will concatenate all the bed files and sort them
cmd <- glue::glue("cat {paste(peak_files_paths, collapse = ' ')} > {opt$output_directory}/peaks.{opt$transcription_factor}_{opt$tissue}.bed")
# print(cmd)
system(cmd)

system(
    glue::glue("sort -k1,1 -k2,2n {opt$output_directory}/peaks.{opt$transcription_factor}_{opt$tissue}.bed > {opt$output_directory}/peaks.{opt$transcription_factor}_{opt$tissue}.sorted.bed")
)

# then I create shuffle these across the genome to get a null

system(
    glue::glue("bedtools shuffle -i {opt$output_directory}/peaks.{opt$transcription_factor}_{opt$tissue}.sorted.bed -chrom -seed 2026 -chromFirst -noOverlapping -g /beagle3/haky/users/temi/projects/TFPred-snakemake/info/hg38.chrom.sizes.sorted -excl {opt$output_directory}/peaks.{opt$transcription_factor}_{opt$tissue}.sorted.bed > {opt$output_directory}/nullpeaks.{opt$transcription_factor}_{opt$tissue}.sorted.bed")
)

nopeaks.granges <- data.table::fread(glue::glue("{opt$output_directory}/nullpeaks.{opt$transcription_factor}_{opt$tissue}.sorted.bed"), col.names = c('chrom', 'start', 'end', 'id_num', 'width', 'strand', 'score', 'neglog10p', 'neglog10q', 'value')) %>% dplyr::filter(chrom %in% valid_chromosomes) %>%
    with(GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges(start, end), strand = '*')) %>% GenomicRanges::reduce()

peaks.granges <- data.table::fread(glue::glue("{opt$output_directory}/peaks.{opt$transcription_factor}_{opt$tissue}.sorted.bed"), col.names = c('chrom', 'start', 'end', 'id_num', 'width', 'strand', 'score', 'neglog10p', 'neglog10q', 'value')) %>% dplyr::filter(chrom %in% valid_chromosomes) %>%
    with(GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges(start, end), strand = '*')) %>% GenomicRanges::reduce()  

# find and remove overlaps
ovlaps <- findOverlaps(query = peaks.granges, subject = nopeaks.granges) 

if(length(ovlaps) > 0){
    nopeaks.granges <- nopeaks.granges[-subjectHits(ovlaps)]
    peaks.granges <- peaks.granges[-queryHits(ovlaps)]
}


# add metadata column of bound and unbound
mcols(nopeaks.granges)$binding_class <- 0
mcols(peaks.granges)$binding_class <- 1

training_data.granges <- unique(c(peaks.granges, nopeaks.granges)) %>% GenomicRanges::resize(fix = 'center', width = 512)

avl_chroms <- seqnames(training_data.granges)@values %>% as.character() %>% unique()

if(!any(avl_chroms %in% c('chr9', 'chr22'))){
    print(glue('INFO - There are no training or test observations for {opt$transcription_factor} in {opt$tissue}'))

    # create empty files
    columns <- c('chrom', 'start', 'end', 'locus', 'binding_class', 'split')

    # Create an empty data.table with 0 rows but specified columns
    empty_dt <- data.table(matrix(nrow = 0, ncol = length(columns)))
    colnames(empty_dt) <- columns

    empty_dt %>% data.table::fwrite(glue("{opt$info_file}"), sep = '\t', quote = F, col.names = T, row.names = F)

    empty_dt %>%
        dplyr::select(locus, binding_class, split) %>%
        data.table::fwrite(glue("{opt$ground_truth_file}"), sep = '\t', quote = F, col.names = T, row.names = F)

    empty_dt %>%
        dplyr::select(locus) %>% 
        data.table::fwrite(glue("{opt$predictors_file}"), sep = '\t', quote = F, col.names = F, row.names = F)

    quit(status = 0)
}


training_data.dt <- training_data.granges %>% plyranges::mutate(
    chrom = as.character(seqnames)
  ) %>% plyranges::mutate(split =
    case_when(
      !chrom %in% c('chr9', 'chr22') ~ "train",
      chrom %in% c('chr9', 'chr22') ~ "test" # Default case
    )) %>% as.data.table() %>% dplyr::select(-seqnames) %>% dplyr::relocate(chrom)

training_data.dt <- training_data.dt %>% dplyr::filter(chrom %in% valid_chromosomes)

set.seed(seed)
training_set <- training_data.dt %>% 
    dplyr::filter(split == 'train') %>%
    dplyr::group_by(binding_class) %>%
    dplyr::slice_sample(n = opt$num_train_observations) %>%
    dplyr::ungroup() %>%
    dplyr::slice_sample(n = nrow(.)) %>% 
    dplyr::mutate(locus = paste(chrom, start, end, sep = '_'))

test_set <- training_data.dt %>% 
    dplyr::filter(split == 'test') %>%
    dplyr::group_by(binding_class) %>%
    dplyr::slice_sample(n = opt$num_test_observations) %>%
    dplyr::ungroup() %>%
    dplyr::slice_sample(n = nrow(.)) %>% 
    dplyr::mutate(locus = paste(chrom, start, end, sep = '_'))

selected_train <- dplyr::bind_rows(training_set, test_set)

# save 
training_data.dt %>%
    data.table::fwrite(glue("{opt$info_file}"), sep = '\t', quote = F, col.names = T, row.names = F)

selected_train %>%
    dplyr::select(locus, binding_class, split) %>%
    data.table::fwrite(glue("{opt$ground_truth_file}"), sep = '\t', quote = F, col.names = T, row.names = F)

selected_train %>%
    dplyr::select(locus) %>% 
    data.table::fwrite(glue("{opt$predictors_file}"), sep = '\t', quote = F, col.names = F, row.names = F)


# generate summary statistics ======
if(!file.exists(opt$summary_file)){
    if(!file.exists(dirname(opt$summary_file))){
        if(!dir.exists(dirname(opt$summary_file))){
            dir.create(dirname(opt$summary_file), recursive = TRUE)
        }
    }

    diagfile <- file(opt$summary_file, open = "a")
    cat("## Training binding sites distribution", file = diagfile, sep = '\n')

} else {
    diagfile <- NULL
}

if(!is.null(opt$diagnostics_file)){
    close(diagfile)
}