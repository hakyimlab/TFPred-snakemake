# Author: Temi
# Date: Thursday July 27 2023
# Description: script to train elastic net TFPred models
# Usage: Rscript train_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--train_data_file", help='data to train with enet'),
    make_option("--rds_file", help='.rds file to be created as the model'),
    make_option("--nfolds", type="integer", default=5L, help='How many cv folds?')
)

opt <- parse_args(OptionParser(option_list=option_list))

# opt <- list()

# bpath <- '/beagle3/haky/users/temi/projects/TFPred-snakemake'
# opt$train_data_file <- file.path(bpath, 'data/ENPACT_5_2024-06-12/aggregation_folder/train_ENPACT_5_2024-06-12_aggByCollect.ARNTL_Bone.prepared.csv.gz')
# opt$rds_file <- file.path(bpath, 'output/models/cistrome_AR_Breast_2023-08-10/aggByCollect_AR_Breast')
# opt$nfolds <- 5

# /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/bin/Rscript workflow/src/train_enet.R --train_data_file data/ENPACT_5_2024-06-12/aggregation_folder/train_ENPACT_5_2024-06-12_aggByCollect.ARNTL_Bone.prepared.csv.gz --rds_file data/ENPACT_5_2024-06-12/models/ARNTL_Bone/ARNTL_Bone_2024-06-12 --nfolds 5

library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)

seed <- 2023
if(file.exists(opt$train_data_file)){
    print(glue('INFO - Reading train data...'))
    dt_train <- data.table::fread(opt$train_data_file)
} else {
    stop(glue('ERROR - Training data cannot be found.'))
}

# remove missing values
cc <- complete.cases(dt_train)
dt_train <- dt_train[cc, ]

# split the data
X_train <- dt_train[, -c(1,2,3)] |> as.matrix()
y_train <- dt_train[, c(1,2,3)] |> as.data.frame()
binding_counts <- y_train$binding_counts
binding_class <- ifelse(y_train$binding_class > 0, 1, 0)
#nbc <- (vbc - min(vbc))/(max(vbc) - min(vbc)) # min-max normalization
# nbc <- log10(1 + y_train$mean_intensity)

cl <- 12 #parallel::makeCluster(5)
print(glue('INFO - Found {parallel::detectCores()} cores but using {cl}'))

set.seed(seed)

doParallel::registerDoParallel(cl)
print(glue('INFO - training enet model'))

train_methods <- c('linear', 'logistic')

parallel::mclapply(train_methods, function(each_method){

    cl <- 6 #parallel::makeCluster(5)
    doParallel::registerDoParallel(cl)

    print(glue('INFO - Starting to build {each_method} enet model'))

    if(each_method == 'linear'){

        cv_model <- tryCatch({
            glmnet::cv.glmnet(x=X_train, y=binding_counts, family = "gaussian", type.measure = "mse", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds)
        }, error = function(e){
            print(glue('ERROR - {e}'))
            return(NULL)
        })
        save_name <- paste0(opt$rds_file, '.linear.rds', sep='') #gsub('.rds', '.linear.rds', opt$rds_file)
    } else if (each_method == 'logistic'){

        cv_model <- tryCatch({
            glmnet::cv.glmnet(x=X_train, y=binding_class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds, trace.it=F)
        }, error = function(e){
            print(glue('ERROR - {e}'))
            return(NULL)
        })
        save_name <-  paste0(opt$rds_file, '.logistic.rds', sep='') #gsub('.rds', '.logistic.rds', opt$rds_file)
    }
    print(cv_model)
    print(glue('INFO - Saving `{save_name}`'))
    #rds_file <- glue('{model_file_basename}.{each_method}.rds')
    saveRDS(cv_model, file=save_name)
    doParallel::stopImplicitCluster()

}, mc.cores=2)

print(glue('INFO - Finished with model training and saving'))
doParallel::stopImplicitCluster()