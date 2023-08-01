# Author: Temi
# Date: Thursday July 27 2023
# Description: script to evaluate TFPred models on train and test data
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--model", help='A TFPred model'),
    make_option("--train_data_file", help='training data file'),
    make_option("--test_data_file", help='test data file'),
    make_option("--eval_output", help='evaluation file in with .rds extension')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)

model <- readRDS(glue('{opt$model}'))
# read in the newx data : train or test
newx_list <- purrr::map(.x=c(opt$train_data_file, opt$test_data_file), function(each_data){
    mat_dt <- data.table::fread(each_data)
    newx <- as.matrix(mat_dt[, -c(1:2)])
    # you only need one : link
    link_pred <- predict(model, newx, s = "lambda.1se", type = 'link') |> as.vector()
    df <- mat_dt[, c(1:2)] |> as.data.frame()
    df$prediction <- link_pred
    colnames(df) <- c('locus', 'peakActivityScore', 'predicted')
    return(df)
}, .progress=T)

names(newx_list) <- c('train_eval', 'test_eval')
# save the object to be read later
# print(glue('INFO - Saving `{agg_methods}_{data_type}_evaluation.{model_type}.rds` to {output_dir}'))
# rds_file <- glue('{output_dir}/{agg_methods}_{data_type}_evaluation.{model_type}.rds')
saveRDS(newx_list, file=opt$eval_output)

# for(model_type in c('linear', 'logistic')){
#     # read in the models
#     # load all the models coefficients
#     models_list <- purrr::map(.x=agg_methods, function(each_method){
#         model <- readRDS(glue('{model_basename}.{model_type}.rds'))
#         return(model)
#     }, .progress=T)
#     names(models_list) <- agg_methods

#     # read in the newx data : train or test
#     newx_list <- purrr::map(.x=agg_methods, function(each_method){
#         #mat_file <- glue('{data_basename}.prepared.csv')
#         mat_dt <- data.table::fread(data_file)
#         return(mat_dt)
#     }, .progress=T)
#     names(newx_list) <- agg_methods

#     # predict
#     predictions_list <- parallel::mclapply(agg_methods, function(each_method){

#         newx <- as.matrix(newx_list[[each_method]][, -c(1:4)])
#         if(model_type == 'linear'){

#             # you only need one : link
#             link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
#             df <- newx_list[[each_method]][, c(1:4)] |> as.data.frame()
#             df$prediction_link <- link_pred

#         } else if (model_type == 'logistic'){
#             # return two things

#             link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
#             response_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'response') |> as.vector()

#             df <- newx_list[[each_method]][, c(1:4)] |> as.data.frame()
#             df$prediction_link <- link_pred
#             df$prediction_response <- response_pred
#         }

#         return(df)

#     }, mc.cores=length(agg_methods))

#     names(predictions_list) <- agg_methods

#     # save the object to be read later
#     print(glue('INFO - Saving `{agg_methods}_{data_type}_evaluation.{model_type}.rds` to {output_dir}'))
#     rds_file <- glue('{output_dir}/{agg_methods}_{data_type}_evaluation.{model_type}.rds')
#     saveRDS(predictions_list, file=rds_file)
# }