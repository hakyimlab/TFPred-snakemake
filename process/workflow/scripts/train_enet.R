



arguments <- commandArgs(trailingOnly=TRUE)
#print(arguments)

library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)

data_file <- arguments[1]
model_file_basename <- arguments[2]


if(file.exists(data_file)){
    print(glue('INFO - Training data exists.'))
} else {
    stop(glue('ERROR - Training data cannot be found.'))
}

#id_data <- arguments[2]
#TF <- arguments[3]
#metainfo <- arguments[4] 
#output_dir <- arguments[5]
#training_date <- arguments[6]

# print(glue('id is {id_data}\nTF is {TF}\nmetainfo is {metainfo}\noutput directory is {output_dir}\ntraining date is {training_date}\n\n'))

dt_train <- data.table::fread(data_file)
print(dim(dt_train))

# split the data
X_train <- dt_train[, -c(1,2,3,4)] |> as.matrix()
y_train <- dt_train[, c(1,2,3,4)] |> as.data.frame()
print(head(y_train))

print(glue('INFO - Found {parallel::detectCores()} cores\n\n'))

# build two models
set.seed(2023)
train_methods <- c('linear', 'logistic')

parallel::mclapply(train_methods, function(each_method){

    cl <- 6 #parallel::makeCluster(5)
    doParallel::registerDoParallel(cl)

    print(glue('INFO - Starting to build {each_method} enet model\n\n'))

    if(each_method == 'linear'){

        cv_model <- glmnet::cv.glmnet(x=X_train, y=y_train$norm_bc, family = "gaussian", type.measure = "mse", alpha = 0.5, keep=T, parallel=T, nfolds=5)
        print(cv_model)

    } else if (each_method == 'logistic'){

        cv_model <- glmnet::cv.glmnet(x=X_train, y=y_train$class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=5, trace.it=F)
        print(cv_model)
    }

    print(glue('INFO - Saving `{model_file_basename}.{each_method}.rds`'))
    rds_file <- glue('{model_file_basename}.{each_method}.rds')
    saveRDS(cv_model, file=rds_file)

    doParallel::stopImplicitCluster()

}, mc.cores=2)

print(glue('INFO - Finished with model training and saving\n\n'))
doParallel::stopImplicitCluster()



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