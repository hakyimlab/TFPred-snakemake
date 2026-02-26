# Author: Temi
# Date: Thursday July 27 2023
# Description: script to evaluate Enpact models on train and test data
# Usage: Rscript evaluate_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--logistic_model", help='An Enpact model'),
    make_option("--train_data_file", help='training data file'),
    make_option("--test_data_file", help='test data file'),
    make_option("--eval_output", help='evaluation file in without an extension; extension will be appended')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)

models <- list()
models[['logistic']] <- readRDS(opt$logistic_model)

for(i in seq_along(models)){
    # read in the newx data : train or test
    newx_list <- purrr::map(.x=c(opt$train_data_file, opt$test_data_file), function(each_data){
        mat_dt <- data.table::fread(each_data)
        newx <- mat_dt |> dplyr::select(starts_with("f_")) |> as.matrix()
        
        # you only need one : link
        link_pred <- predict(models[[i]], newx, s = "lambda.1se", type = 'link') |> as.vector()
        response_pred <- predict(models[[i]], newx, s = "lambda.1se", type = 'response') |> as.vector()
        df <- mat_dt |> dplyr::select(!starts_with("f_")) |> as.data.frame()
        df$TFPred_score <- link_pred
        df$probability <- response_pred
        colnames(df) <- c('locus', 'binding_class', 'Enpact_score', 'probability')
        return(df)
    }, .progress=T)

    data.table::fwrite(newx_list[[1]], glue('{opt$eval_output}.{names(models)[i]}.train_eval.txt.gz'), sep='\t', compress='gzip', quote=F, row.names=F)
    data.table::fwrite(newx_list[[2]], glue('{opt$eval_output}.{names(models)[i]}.test_eval.txt.gz'), sep='\t', compress='gzip', quote=F, row.names=F)
}