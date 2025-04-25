# Author: Temi
# Date: Thursday July 27 2023
# Description: split the aggregated enformer outputs into train and test
# Usage: Rscript train_test_split.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--data_file", help='A transcription factor e.g. AR'),
    make_option("--ground_truth_file", help='A tissue e.g. Breast'),
    make_option("--aggregation_method", help="Predicted motif file, particularly from HOMER"),
    make_option("--train_prepared_file", help='The folder where the cistrome bedfiles are located'),
    make_option("--test_prepared_file", help='A folder where relevant bedfiles will be linked')
)

opt <- parse_args(OptionParser(option_list=option_list))

seed <- 2023
set.seed(seed)

library(glue)
library(tidyverse)
library(data.table)

# opt <- list()

# opt$data_file <- 'data/ENPACT_734_2025-04-24/aggregation_folder/ENPACT_734_aggByCollect_NEUROG2_FetalLung.csv.gz' 
# opt$ground_truth_file <-'data/ENPACT_734_2025-04-24/predictor_files/NEUROG2_FetalLung.ground_truth.txt'
# opt$aggregation_method <- 'aggByCollect' 

# --train_prepared_file data/ENPACT_734_2025-04-24/aggregation_folder/train_ENPACT_734_2025-04-24_aggByCollect.NEUROG2_FetalLung.prepared.csv.gz --test_prepared_file data/ENPACT_734_2025-04-24/aggregation_folder/test_ENPACT_734_2025-04-24_aggByCollect.NEUROG2_FetalLung.prepared.csv.gz

ground_truth <- data.table::fread(opt$ground_truth_file) %>%
  dplyr::select(locus, binding_class, binding_counts, split) #%>% dplyr::rename(locus=V1, peakActivityScore=V2)

find_duplicates_in_dataframe <- function(dt, col, return_dups=TRUE){
  n_occur <- data.frame(table(dt[[col]]))
  if(return_dups == TRUE){
    return(dt[dt[[col]] %in% n_occur$Var1[n_occur$Freq > 1],])
  } else {
    return(n_occur[n_occur$Freq > 1,])
  }
}
#find_duplicates_in_dataframe(gt, col='V1') # there should be no duplicates

center_dt <- data.table::fread(opt$data_file, fill=T)
gt <- ground_truth[ground_truth$locus %in% center_dt$id, ]
gt_dedup <- gt[!duplicated(gt[['locus']]), ]
new_dt <- merge(gt_dedup, center_dt, by.x='locus', by.y='id')
colnames(new_dt) <- c('locus', 'binding_class', 'binding_counts', 'split', paste('f_', 1:(ncol(new_dt)-4), sep=''))

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
print(glue("INFO - {opt$aggregation_method}'s train (and test, if applicable) data have been saved."))



# purrr::map(.x=names(agg_transform_list), function(each_method){
#     each_dt <- agg_transform_list[[each_method]]
#     tr_size <- ceiling(nrow(each_dt) * 0.8)
#     tr_indices <- sample(1:nrow(each_dt), tr_size)
#     train <- each_dt[tr_indices, ]
#     test <- each_dt[-tr_indices, ]
#     data.table::fwrite(x=train, file=train_prepared_file, quote=F, row.names=F)
#     # just in case there is a test data
#     if(nrow(test) > 0){
#         data.table::fwrite(x=test, file=test_prepared_file, quote=F, row.names=F)
#     }
#     print(glue("INFO - {each_method}'s train (and test, if applicable) data have been saved."))
# })




# purrr::map(.x=names(agg_transform_list), function(each_method){
#     each_dt <- agg_transform_list[[each_method]]
#     tr_size <- ceiling(nrow(each_dt) * 0.8)
#     tr_indices <- sample(1:nrow(each_dt), tr_size)
#     train <- each_dt[tr_indices, ]
#     test <- each_dt[-tr_indices, ]
#     data.table::fwrite(x=train, file=train_prepared_file_without_gz, quote=F, row.names=F)
#     cmd <- glue("pigz -kf {train_prepared_file_without_gz}")
#     system(cmd)
#     # just in case there is a test data
#     if(nrow(test) > 0){
#         data.table::fwrite(x=test, file=test_prepared_file_without_gz, quote=F, row.names=F)
#         cmd <- glue("pigz -kf {train_prepared_file_without_gz}")
#         system(cmd)
#     }
#     print(glue("INFO - {each_method}'s train (and test, if applicable) data have been saved."))
# })

# temp_dt <- data.table::fread(glue('{model_data_dir}/train_aggByPreCenter.csv.gz'))
# temp_dt[1:5, 1:5] ; temp_dt |> dim() ; temp_dt$class |> table() 

# temp_dt <- data.table::fread(glue('{model_data_dir}/test_aggByPreCenter.csv.gz'))
# temp_dt[1:5, 1:5] ; temp_dt |> dim() ; temp_dt$class |> table()

# rm(list = c('temp_dt', 'agg_transform_list'))