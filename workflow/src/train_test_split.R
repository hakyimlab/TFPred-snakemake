# Author: Temi
# Date: Thursday July 27 2023
# Description: split the aggregated enformer outputs into train and test
# Usage: Rscript train_test_split.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--data_file", help='[Input] A file with the epigenomic features'),
    make_option("--ground_truth_file", help='[Input] A file mapping the loci to the ground truth'),
    make_option("--train_prepared_file", help='[Output] Training file'),
    make_option("--test_prepared_file", help='[Output] Test file')
)

opt <- parse_args(OptionParser(option_list=option_list))

seed <- 2023
set.seed(seed)

library(glue)
library(tidyverse)
library(data.table)

print(opt)

ground_truth <- data.table::fread(opt$ground_truth_file) %>% dplyr::select(locus, binding_class, split)

center_dt <- read.csv(opt$data_file, fill=T)
gt <- ground_truth[ground_truth$locus %in% center_dt$id, ]
gt_dedup <- gt[!duplicated(gt[['locus']]), ]
new_dt <- merge(gt_dedup, center_dt, by.x='locus', by.y='id')
colnames(new_dt) <- c('locus', 'binding_class', 'split', paste('f_', 1:(ncol(new_dt)-3), sep=''))

train <- new_dt %>% 
  dplyr::filter(split == 'train') %>%
  dplyr::select(-split)
test <- new_dt %>%
  dplyr::filter(split == 'test') %>%
  dplyr::select(-split)

data.table::fwrite(x=train, file=opt$train_prepared_file, quote=F, row.names=F, compress='gzip')
# just in case there is a test data
if(nrow(test) > 0){
    data.table::fwrite(x=test, file=opt$test_prepared_file, quote=F, row.names=F, compress = 'gzip')
}