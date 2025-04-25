# Author: Temi
# Date: Friday July 26 2024
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--transcription_factors", help='A transcription factor e.g. AR'),
    make_option("--tissues", help='A tissue e.g. Breast'),
    make_option("--path_pattern", help='pattern to match path e.g. /beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2024-07-26/evaluation/{1}_{2}_2024-07-26.logistic.{3}_eval.txt.gz'),
    make_option("--model_path", help='pattern to match path e.g. /beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2024-07-26/models/{1}_{2}/{1}_{2}_2024-07-26.logistic.rds'),
    make_option("--statistics_file", help='file to write statistics to'),
    make_option("--weights_file_basename", help='file to write weights to'),
    make_option("--training_peaks_directory", help='file to write training statistics to')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
# library(pROC)
library(stringr)
library(foreach)
library(glmnet)
library(glue)

# /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/bin/Rscript workflow/src/compile_statistics.R --transcription_factors AR,AR,AR,AR --tissues Breast,EmbryonicKidney,MammaryGland,Prostate --path_pattern data/ENPACT_AR_2025-02-11/evaluation/{1}_{2}_2025-02-11.linear.{3}_eval.txt.gz --statistics_file data/ENPACT_AR_2025-02-11/statistics/ENPACT_AR_2025-02-11.compiled_stats.txt --model_path data/ENPACT_AR_2025-02-11/models/{1}_{2}/{1}_{2}_2025-02-11.linear.rds --weights_file_basename data/ENPACT_AR_2025-02-11/statistics/ENPACT_AR_2025-02-11.compiled_weights --training_peaks_directory data/ENPACT_AR_2025-02-11/sortedbeds

# NFYA_BoneMarrow
# opt <- list()
# opt$transcription_factors <- 'FOXA1,MYC,YY1,CTCF,GATA4,TAL1,SP1,MYC,MYC,VDR,CTCF,MYC,NEUROG2,CTCF,GATA2,PPARG,FOSL2,CTCF,CTCF,CTCF,CTCF,NANOG,FOXA2,ERG,REST,FOXA1,POU5F1,GATA1,FOXM1,GATA3,GATA3,CTCF,FOXA1,CTCF,CTCF,SPI1,RUNX1,CTCF,HOXB13,AR,CTCF,GATA2,E2F1,ESR1,NR3C1,ESR1,RELA,NR3C1'
# opt$tissues <- 'MammaryGland,Breast,Blood,EmbryonicKidney,Embryo,BoneMarrow,Colon,Blood,Colon,Blood,MammaryGland,BoneMarrow,FetalLung,Liver,BoneMarrow,Adipose,Lung,Brain,Skin,UmbilicalVein,Embryo,Embryo,Skin,Prostate,Blood,Breast,Embryo,BoneMarrow,Breast,MammaryGland,Breast,Lung,Prostate,Colon,BoneMarrow,Blood,Blood,Breast,Prostate,Prostate,Blood,Prostate,Prostate,MammaryGland,Lung,Breast,Blood,Breast'
# opt$path_pattern <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/evaluation/{1}_{2}_2025-04-24.logistic.{3}_eval.txt.gz'
# opt$statistics_file <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/statistics/ENPACT_734_2025-04-24.compiled_stats.txt'
# opt$model_path <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/models/{1}_{2}/{1}_{2}_2025-04-24.logistic.rds'
# opt$training_peaks_directory <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/sortedbeds'

tflist <- strsplit(opt$transcription_factors, ',')[[1]]
tissuelist <- strsplit(opt$tissues, ',')[[1]]

pp <- c("\\{1\\}", "\\{2\\}", "\\{3\\}")

qMatrix <- rbind(
    cbind(tflist, tissuelist, 'train'),
    cbind(tflist, tissuelist, 'test')
) |> as.data.frame() |> setNames(c('transcription_factor', 'tissue', 'type')) %>%
    dplyr::arrange(transcription_factor, tissue, type)

# register cores since I want to parallelize
num_clusters <- 4 #- 5 # 12 - 5
doParallel::registerDoParallel(num_clusters)

cat('INFO - Registering', num_clusters, 'clusters for a parallel run\n')


calculate_logistic_metrics <- function(dt, transcription_factor, tissue, type){

    tryCatch({
        pp <- with(dt, pROC::roc(binding_class, TFPred_score, quiet = TRUE))
        print(pp)
        pp_auc <- pp$auc |> as.numeric() |> round(2)
        pp_var <- pROC::var(pp)|> as.numeric() %>% format(digits = 2)
        pp_ci <- with(dt, pROC::ci.auc(binding_class, TFPred_score, conf.level = 0.95, quiet = T)) |> as.numeric() 
        low <- pp_ci[1] |> round(2)
        upp <- pp_ci[3] |> round(2)

        # t-test
        ttest <- t.test(TFPred_score ~ binding_class, data=dt)
        ttest_pvalue <- ttest$p.value
        # if(ttest$p.value == 0){
        #     ttest_pvalue <- 'p < 2.2e-16'
        # } else {
        #     ttest_pvalue <- paste('p =', format(ttest$p.value, digits = 2))
        # }
        
        return(cbind(transcription_factor, tissue, type, pp_auc, pp_var, low, upp, ttest_pvalue))
    }, error = function(e){
        return(cbind(transcription_factor, tissue, type, NA, NA, NA, NA, NA))
    })
   
}


# calculate_logistic_metrics(dt, 'AR', 'Prostate', 'test')

# dt <- data.table::fread('/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/evaluation/AR_Prostate_2025-04-24.logistic.test_eval.txt.gz')


# tt <- t.test(TFPred_score ~ binding_class, data=dt)
# tt$p.value

# n_training_files <- list.files(file.path(opt$training_peaks_directory, 'AR_Prostate'), '^peaks_.*.bed$', full.names = T)
# n_training_files <- length(n_training_files[file.info(n_training_files)$size > 0])

# 2 * pt(q = abs(tt$statistic), df = tt$parameter, lower=FALSE, log.p = T)

evaluations <- foreach::foreach(i=seq_len(nrow(qMatrix)), .combine='rbind', .inorder=T) %dopar% {
    # do something
    qi <- qMatrix[i, ] |> unlist()
    tfile <- stringr::str_replace_all(opt$path_pattern, purrr::set_names(qi, pp))
    mpath <- stringr::str_replace_all(opt$model_path, purrr::set_names(qi[1:2], pp[1:2]))

    # count number of files in training peaks directory
    n_training_files <- list.files(file.path(opt$training_peaks_directory, paste0(qi[1], '_', qi[2])), '^peaks_.*.bed$', full.names = T)
    # filter if size == 0
    n_training_files <- length(n_training_files[file.info(n_training_files)$size > 0])
    if(file.exists(tfile)){
        dt <- data.table::fread(tfile)
        if(any(nrow(dt) == 0)){
            evaluations <- data.frame(assay = qi[1], context = qi[2], type = qi[3], auc = NA, var = NA, low = NA, upp = NA, ttest_pvalue=NA, num_training_files = NA, model = NA, path = NA, num_0 = NA, num_1 = NA)
        } else if (any(nrow(dt) > 0)) {
            evaluations <- calculate_logistic_metrics(dt, qi[1], qi[2], qi[3])
            evaluations <- evaluations %>% as.data.frame() %>%
                stats::setNames(c('assay', 'context', 'type', 'auc', 'var', 'low', 'upp', 'ttest_pvalue')) %>%
                dplyr::mutate(num_training_files = n_training_files, 
                    model = dplyr::case_when(
                    file.exists(mpath) ~ paste0(qi[1], '_', qi[2]),
                    .default = NA
                ), path = dplyr::case_when(
                    file.exists(mpath) ~ mpath,
                    .default = NA
                ))

            dist <- dt %>%
                dplyr::group_by(binding_class) %>%
                dplyr::summarise(n = dplyr::n()) %>%
                tidyr::pivot_wider(names_from = binding_class, values_from = n) %>%
                dplyr::rename('num_0' = `0`, 'num_1' = `1`)
            evaluations$num_0 <- dist$num_0
            evaluations$num_1 <- dist$num_1

            evaluations <- cbind(evaluations)
        }
    }

    return(evaluations)
}

data.table::fwrite(evaluations, file=opt$statistics_file, sep='\t', row.names=F, col.names=T, quote=F)

# read in the models and extract the weights

models_mtdt <- evaluations %>%
    dplyr::filter(!is.na(model)) %>%
    dplyr::select(model, path) %>%
    dplyr::distinct()

models_weights <- foreach::foreach(i=seq_len(nrow(models_mtdt)), .combine='rbind', .inorder=T) %do% {
    mtdt <- models_mtdt[i, ] |> unlist()
    mpath <- mtdt[2]; model <- mtdt[1]
    model <- readRDS(mpath)
    weights <- coef(model, s = 'lambda.1se')[-1] #|> as.data.frame() |> setNames(model)
    return(weights)
} %>% t() %>% as.data.frame() %>% 
    setNames(models_mtdt$model) %>% 
    dplyr::mutate(feature = paste0('f_', seq_len(nrow(.)))) %>% dplyr::relocate(feature)

data.table::fwrite(models_weights, file=glue("{opt$weights_file}.lambda.1se.txt.gz"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')


models_weights <- foreach::foreach(i=seq_len(nrow(models_mtdt)), .combine='rbind', .inorder=T) %do% {
    mtdt <- models_mtdt[i, ] |> unlist()
    mpath <- mtdt[2]; model <- mtdt[1]
    model <- readRDS(mpath)
    weights <- coef(model, s = 'lambda.min')[-1] #|> as.data.frame() |> setNames(model)
    return(weights)
} %>% t() %>% as.data.frame() %>% 
    setNames(models_mtdt$model) %>% 
    dplyr::mutate(feature = paste0('f_', seq_len(nrow(.)))) %>% dplyr::relocate(feature)

data.table::fwrite(models_weights, file=glue("{opt$weights_file}.lambda.min.txt.gz"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')

doParallel::stopImplicitCluster()
# m <- readRDS('/project/haky/users/temi/Enpact-figures/data/aggByCollect_AR_Prostate.logistic.rds')

# coef(m, s = m$lambda.min)[-1] %>% as.data.frame() %>% setNames('AR_Prostate') %>% dplyr::mutate(feature = paste0('f_', seq_len(nrow(.)))) %>% dplyr::arrange(desc(abs(AR_Prostate))) %>% dplyr::relocate(feature) %>% head(10)













# calculate_linear_metrics <- function(dt, transcription_factor, tissue, type){

#     # remove missing values
#     cc <- complete.cases(dt)
#     dt <- dt[cc, ]

#     tryCatch({
#         # calculate mse
#         #pp <- with(dt, mean((mean_intensity - TFPred_score)**2))
#         sq_err <- with(dt, (binding_counts - TFPred_score)**2)
#         pp_mse <- mean(sq_err) %>% round(2)
#         pp_var <- var(sq_err) %>% format(digits = 2)
#         pp_ci <- with(dt, stats::t.test(binding_counts, TFPred_score)$conf.int) 
#         low <- pp_ci[1] |> round(2)
#         upp <- pp_ci[2] |> round(2)

#         # t-test
#         #ttest <- t.test(TFPred_score ~ binding_class, data=dt)
#         ttest_pvalue <- NA #ttest$p.value
#         # if(ttest$p.value == 0){
#         #     ttest_pvalue <- 'p < 2.2e-16'
#         # } else {
#         #     ttest_pvalue <- paste('p =', format(ttest$p.value, digits = 2))
#         # }
        
#         return(cbind(transcription_factor, tissue, type, pp_mse, pp_var, low, upp, ttest_pvalue))
#     }, error = function(e){
#         return(cbind(transcription_factor, tissue, type, NA, NA, NA, NA, NA))
#     })
   
# }








# calculate_metrics <- function(dt, transcription_factor, tissue, type){

#     tryCatch({
#         pp <- with(dt, pROC::roc(class, TFPred_score, quiet = TRUE))
#         pp_auc <- pp$auc |> as.numeric() |> round(2)
#         pp_var <- pROC::var(pp)|> as.numeric() %>% format(digits = 2)
#         pp_ci <- with(dt, pROC::ci.auc(binding_class, TFPred_score, conf.level = 0.95, quiet = T)) |> as.numeric() 
#         low <- pp_ci[1] |> round(2)
#         upp <- pp_ci[3] |> round(2)

#         # t-test
#         ttest <- t.test(TFPred_score ~ binding_class, data=dt)
#         ttest_pvalue <- ttest$p.value
#         # if(ttest$p.value == 0){
#         #     ttest_pvalue <- 'p < 2.2e-16'
#         # } else {
#         #     ttest_pvalue <- paste('p =', format(ttest$p.value, digits = 2))
#         # }
        
#         return(cbind(transcription_factor, tissue, type, pp_auc, pp_var, low, upp, ttest_pvalue))
#     }, error = function(e){
#         return(cbind(transcription_factor, tissue, type, NA, NA, NA, NA, NA))
#     })
   
# }




# if(file.exists(opt$statistics_file)){
#     write.table(evaluations, file=opt$statistics_file, sep='\t', row.names=F, col.names=F, quote=F, append = T)
# } else if(!file.exists(opt$statistics_file)){
    
# }


# tr_dt <- data.table::fread(opt$train_evaluation)
# te_dt <- data.table::fread(opt$test_evaluation)



# if(any(nrow(tr_dt) == 0) || any(nrow(te_dt) == 0)){
#     evaluations <- data.frame(transcription_factor = opt$transcription_factor, tissue = opt$tissue, type = c('train', 'test'), auc = NA, var = NA, low = NA, upp = NA, num_0 = NA, num_1 = NA)
# } else if (any(nrow(tr_dt) > 0) && any(nrow(te_dt) > 0)){
#     train_evaluations <- calculate_metrics(tr_dt, 'train')
#     test_evaluations <- calculate_metrics(te_dt, 'test')
#     evaluations <- rbind(train_evaluations, test_evaluations) |> as.data.frame() %>%
#         setNames(c('transcription_factor', 'tissue', 'type', 'auc', 'var', 'low', 'upp'))

#     train_dist <- tr_dt %>%
#         dplyr::group_by(binding_class) %>%
#         dplyr::summarise(n = n()) %>%
#         tidyr::pivot_wider(names_from = binding_class, values_from = n) %>%
#         dplyr::rename('num_0' = `0`, 'num_1' = `1`)
#     test_dist <- te_dt %>%
#         dplyr::group_by(binding_class) %>%
#         dplyr::summarise(n = n()) %>%
#         tidyr::pivot_wider(names_from = binding_class, values_from = n) %>%
#         dplyr::rename('num_0' = `0`, 'num_1' = `1`)

#     dist <- rbind(train_dist, test_dist) |> as.data.frame()

#     evaluations <- cbind(evaluations, dist)
# }

