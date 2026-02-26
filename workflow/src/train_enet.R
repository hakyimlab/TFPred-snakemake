# Author: Temi
# Date: Thursday July 27 2023
# Description: script to train elastic net Enpact models
# Usage: Rscript train_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--train_data_file", help='data to train with enet'),
    make_option("--rds_file", help='.rds file to be created as the model'),
    make_option("--nfolds", type="integer", default=5L, help='How many cv folds?')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)

seed <- 2023
if(file.exists(opt$train_data_file)){
    print(glue('INFO - Reading train data...'))
    dt_train <- tryCatch({
        read.csv(gzfile(opt$train_data_file))
    }, error = function(err){
        print(err)
    }, warning = function(war){
        print(warnings())
        print(war)
    })

} else {
    stop(glue('ERROR - Training data cannot be found.'))
    quit(-1)
}

print(dt_train[1:5, 1:5])
print(dim(dt_train))

# quit(-1)

# remove missing values
cc <- complete.cases(dt_train)
dt_train <- dt_train[cc, ]

# split the data
X_train <- dt_train |> dplyr::select(starts_with("f_")) |> as.matrix()
y_train <- dt_train |> dplyr::select(!starts_with("f_")) |> as.data.frame()

cl <- 12 #parallel::makeCluster(5)
print(glue('INFO - Found {parallel::detectCores()} cores but using {cl}'))

set.seed(seed)

doParallel::registerDoParallel(cl)
print(glue('INFO - training enet model'))

cv_model <- glmnet::cv.glmnet(x=X_train, y=y_train$binding_class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds, trace.it=F)

print(cv_model)

save_name <-  paste0(opt$rds_file, '.logistic.rds', sep='')
print(glue('INFO - Saving the model to `{save_name}`'))
saveRDS(cv_model, file=save_name)
print(glue('INFO - Finished with model training and saving'))
doParallel::stopImplicitCluster()