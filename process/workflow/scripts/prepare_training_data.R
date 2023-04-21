
# Description:
# Author:
# Usage: 


library(glue)
library(tidyverse)

arguments <- commandArgs(trailingOnly=TRUE)

dt_file <- arguments[1]
gt_file <- arguments[2]
agg_method <- arguments[3]
train_prepared_file <- arguments[4]
test_prepared_file <- arguments[5]

# jfile <- jsonlite::fromJSON(glue('{project_dir}/config_files/aggregation_config_{dataset}_{TF}.json'))
# data_date <- jfile$run_date
# TF <- jfile$transcription_factor

# predictor_file <- Sys.glob(glue('{project_dir}/motif_intervals/{dataset}/intervals_{data_date}/ground_truth/{dataset}_{TF}_*.txt'))
# predictor_file

ground_truth <- data.table::fread(gt_file)
# linearize the binding counts
vbc <- ground_truth$V3
nbc <- (vbc - min(vbc))/(max(vbc) - min(vbc))
ground_truth$norm_bc <- nbc

find_duplicates_in_dataframe <- function(dt, col, return_dups=TRUE){
  n_occur <- data.frame(table(dt[[col]]))
  if(return_dups == TRUE){
    return(dt[dt[[col]] %in% n_occur$Var1[n_occur$Freq > 1],])
  } else {
    return(n_occur[n_occur$Freq > 1,])
  }
}
#find_duplicates_in_dataframe(gt, col='V1') # there should be no duplicates

agg_transform_list <- purrr::map(.x=agg_method, function(each_method){
    center_dt <- data.table::fread(dt_file, fill=T)
    gt <- ground_truth[ground_truth$V1 %in% center_dt$id, ]
    gt_dedup <- gt[!duplicated(gt[['V1']]),]
    new_dt <- merge(gt_dedup, center_dt, by.x='V1', by.y='id')
    colnames(new_dt) <- c('region', 'class', 'binding_counts', 'norm_bc', paste('f_', 1:(ncol(new_dt)-4), sep=''))
    return(new_dt)
}, .progress=T)
names(agg_transform_list) <- agg_method

# train_prepared_file <- gsub('.gz', '', train_prepared_file)
# test_prepared_file <- gsub('.gz', '', test_prepared_file)

set.seed(2023)
purrr::map(.x=names(agg_transform_list), function(each_method){
    each_dt <- agg_transform_list[[each_method]]
    tr_size <- ceiling(nrow(each_dt) * 0.8)
    tr_indices <- sample(1:nrow(each_dt), tr_size)
    train <- each_dt[tr_indices, ]
    test <- each_dt[-tr_indices, ]
    data.table::fwrite(x=train, file=train_prepared_file, quote=F, row.names=F)
    # just in case there is a test data
    if(nrow(test) > 0){
        data.table::fwrite(x=test, file=test_prepared_file, quote=F, row.names=F)
    }
    print(glue("INFO - {each_method}'s train (and test, if applicable) data have been saved."))
})




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