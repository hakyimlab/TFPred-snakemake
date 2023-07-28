# Author: Temi
# Date: Thursday July 27 2023
# Description: script to train elastic net TFPred models
# Usage: Rscript train_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--train_data_file", help='data to train with enet'),
    make_option("--rds_file", help='.rds file to be created as the model')
)

opt <- parse_args(OptionParser(option_list=option_list))

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

# split the data
X_train <- dt_train[, -c(1,2)] |> as.matrix()
y_train <- dt_train[, c(1,2)] |> as.data.frame()

cl <- 12 #parallel::makeCluster(5)
print(glue('INFO - Found {parallel::detectCores()} cores but using {cl}\n\n'))

set.seed(seed)

doParallel::registerDoParallel(cl)
print(glue('INFO - training enet model\n\n'))

cv_model <- tryCatch({
    glmnet::cv.glmnet(x=X_train, y=y_train[, 2], family = "gaussian", type.measure = "mse", alpha = 0.5, keep=T, parallel=T, nfolds=10)
}, error = function(e){
    print(glue('ERROR - {e}'))
    return(NULL)
})

print(cv_model)
print(glue('INFO - Saving model to `{opt$rds_file}`'))
if(!dir.exists(basename(opt$rds_file))){dir.create(basename(opt$rds_file))}
saveRDS(cv_model, file=opt$rds_file)
doParallel::stopImplicitCluster()
print(glue('INFO - Finished with model training and saving\n\n'))




# train_methods <- c('linear', 'logistic')

# parallel::mclapply(train_methods, function(each_method){

#     cl <- 6 #parallel::makeCluster(5)
#     doParallel::registerDoParallel(cl)

#     print(glue('INFO - Starting to build {each_method} enet model\n\n'))

#     if(each_method == 'linear'){

#         # cv_model <- glmnet::cv.glmnet(x=X_train, y=y_train$norm_bc, family = "gaussian", type.measure = "mse", alpha = 0.5, keep=T, parallel=T, nfolds=5)
#         # print(cv_model)

#         cv_model <- tryCatch({
#             glmnet::cv.glmnet(x=X_train, y=y_train$norm_bc, family = "gaussian", type.measure = "mse", alpha = 0.5, keep=T, parallel=T, nfolds=5)
#         }, error = function(e){
#             print(glue('ERROR - {e}'))
#             return(NULL)
#         })

#     } else if (each_method == 'logistic'){

#         # cv_model <- glmnet::cv.glmnet(x=X_train, y=y_train$class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=5, trace.it=F)
#         # print(cv_model)

#         cv_model <- tryCatch({
#             glmnet::cv.glmnet(x=X_train, y=y_train$class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=5, trace.it=F)
#         }, error = function(e){
#             print(glue('ERROR - {e}'))
#             return(NULL)
#         })
#     }
#     print(cv_model)
#     print(glue('INFO - Saving `{model_file_basename}.{each_method}.rds`'))
#     rds_file <- glue('{model_file_basename}.{each_method}.rds')
#     saveRDS(cv_model, file=rds_file)

#     doParallel::stopImplicitCluster()

# }, mc.cores=2)

# print(glue('INFO - Finished with model training and saving\n\n'))
# doParallel::stopImplicitCluster()



# register a parallel backend
# cl <- 24
# doParallel::registerDoParallel(cl)

#print(glue('[INFO] Registering {foreach::getDoParWorkers()} workers/cores\n'))


#mixing_parameters <- c(0, 0.25, 0.5, 0.75, 1)
# mixing_parameters <- 0.5
# enet_center_binary_models_list <- parallel::mclapply(mixing_parameters, function(mix){

#     cl <- 5 #parallel::makeCluster(5)
#     doParallel::registerDoParallel(cl)
#     #print(glue('[INFO] Registering {len(foreach::getDoParWorkers())} workers/cores for {mix} mixing parameter\n'))

#     model <- glmnet::cv.glmnet(x=X_train, y=y_train$class, family = "binomial", type.measure = "auc", alpha = mix, keep=T, parallel=T, nfolds=5)
#     print(model)
#     doParallel::stopImplicitCluster()
#     #registerDoSEQ()
#     return(model)

# }, mc.cores=48)



#names(enet_center_binary_models_list) <- mixing_parameters



# vbc <- y_train$binding_counts
# nbc <- (vbc - min(vbc))/(max(vbc) - min(vbc))
# enet_center_linear_model <- glmnet::cv.glmnet(x=X_train, y=nbc, family = "gaussian", type.measure = "mse", alpha = 0.5, keep=T, parallel=T, nfolds=5)
# cat('[INFO] Finished with linear model\n')

# enet_center_multinomial_model <- glmnet::cv.glmnet(x=X_train, y=y_train$binding_counts, family = "multinomial", type.multinomial = "ungrouped", alpha=0.5, keep=T, parallel=T, nfolds=5)
# cat('[INFO] Finished with multinomial model\n')



# save the model

#saveRDS(enet_center_linear_model, file=glue('{project_dir}/models/enet_models/{id_data}_center_linear_{unique_id_data}_{run_date}.rds'))
#saveRDS(enet_center_multinomial_model, file=glue('{project_dir}/models/enet_models/{id_data}_center_multinomial_{unique_id_data}_{run_date}.rds'))

#registerDoSEQ()
#cat('[INFO] Job completed')