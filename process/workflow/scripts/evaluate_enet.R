# Predict on test, train

arguments <- commandArgs(trailingOnly=TRUE)

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)

model_basename <- arguments[1] # ./models
data_file <- arguments[2] # ./data
output_dir <- arguments[3] # ./model_evaluation
data_type <- arguments[4] # train or test

agg_methods <- 'aggByMeanCenter'

for(model_type in c('linear', 'logistic')){
    # read in the models
    # load all the models coefficients
    models_list <- purrr::map(.x=agg_methods, function(each_method){
        model <- readRDS(glue('{model_basename}.{model_type}.rds'))
        return(model)
    }, .progress=T)
    names(models_list) <- agg_methods

    # read in the newx data : train or test
    newx_list <- purrr::map(.x=agg_methods, function(each_method){
        #mat_file <- glue('{data_basename}.prepared.csv')
        mat_dt <- data.table::fread(data_file)
        return(mat_dt)
    }, .progress=T)
    names(newx_list) <- agg_methods

    # predict
    predictions_list <- parallel::mclapply(agg_methods, function(each_method){

        newx <- as.matrix(newx_list[[each_method]][, -c(1:4)])
        if(model_type == 'linear'){

            # you only need one : link
            link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
            df <- newx_list[[each_method]][, c(1:4)] |> as.data.frame()
            df$prediction_link <- link_pred

        } else if (model_type == 'logistic'){
            # return two things

            link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
            response_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'response') |> as.vector()

            df <- newx_list[[each_method]][, c(1:4)] |> as.data.frame()
            df$prediction_link <- link_pred
            df$prediction_response <- response_pred
        }

        return(df)

    }, mc.cores=length(agg_methods))

    names(predictions_list) <- agg_methods

    # save the object to be read later
    print(glue('INFO - Saving `{agg_methods}_{data_type}_evaluation.{model_type}.rds` to {output_dir}'))
    rds_file <- glue('{output_dir}/{agg_methods}_{data_type}_evaluation.{model_type}.rds')
    saveRDS(predictions_list, file=rds_file)
}