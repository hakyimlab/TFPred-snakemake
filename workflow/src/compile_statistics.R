# Author: Temi
# Date: Friday July 26 2024
# Description: script to compile pvalues of t-tests and wilcoxon tests, training metrics, e.t.c
# Usage: Rscript compile_statistics.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--models_list", help='A two-column tsv file (assay, context) listing the TF-tissue pairs'),
    make_option("--path_pattern", help='pattern to match path to evaluation files e.g. /beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2024-07-26/evaluation/{1}_{2}_2024-07-26.logistic.{3}_eval.txt.gz'),
    make_option("--model_path", help='pattern to match path to models e.g. /beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2024-07-26/models/{1}_{2}/{1}_{2}_2024-07-26.logistic.rds'),
    make_option("--statistics_file", help='file to write statistics to'),
    make_option("--weights_file_basename", help='file to write weights to'),
    make_option("--training_peaks_directory", help='file to write training statistics to')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(stringr)
library(foreach)
library(glmnet)
library(glue)
library(ggplot2)

# NFYA_BoneMarrow
# opt <- list()
# opt$models_list <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/metadata/enpact_models_to_train.minimal.tsv'
# opt$path_pattern <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_MINIMAL_2026-01-02/evaluation/{1}_{2}_2026-01-02.logistic.{3}_eval.txt.gz'
# opt$statistics_file <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/statistics/ENPACT_734_2025-04-24.compiled_stats.txt'
# opt$model_path <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_MINIMAL_2026-01-02/models/{1}_{2}/{1}_{2}_2026-01-02.logistic.rds'
# opt$training_peaks_directory <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_MINIMAL_2026-01-02/sortedbeds'

dt_models <- data.table::fread(opt$models_list)
tflist <- dt_models$assay
tissuelist <- dt_models$context

pp <- c("\\{1\\}", "\\{2\\}", "\\{3\\}")

qMatrix <- rbind(
    cbind(tflist, tissuelist, 'train'),
    cbind(tflist, tissuelist, 'test')
) |> as.data.frame() |> setNames(c('transcription_factor', 'tissue', 'type')) %>%
    dplyr::arrange(transcription_factor, tissue, type)

# register cores since I want to parallelize
num_clusters <- 16 #- 5 # 12 - 5
doParallel::registerDoParallel(num_clusters)

cat('INFO - Registering', num_clusters, 'clusters for a parallel run\n')
calculate_logistic_metrics <- function(dt, transcription_factor, tissue, type){

    tryCatch({
        dt <- dt[complete.cases(dt), ]
        pp <- with(dt, pROC::roc(binding_class, TFPred_score, quiet = TRUE))

        # auroc
        pp_auc <- pp$auc |> as.numeric() |> round(3)
        pp_var <- pROC::var(pp)|> as.numeric() %>% format(digits = 3)
        pp_ci <- with(dt, pROC::ci.auc(binding_class, TFPred_score, conf.level = 0.95, quiet = T)) |> as.numeric() 
        low <- pp_ci[1] |> round(3)
        upp <- pp_ci[3] |> round(3)

        # t-test
        ttest <- t.test(TFPred_score ~ binding_class, data=dt)
        ttest_pvalue <- ttest$p.value

        # wilcoxon-test
        wiltest <- wilcox.test(TFPred_score ~ binding_class, data=dt)
        wiltest_pvalue <- wiltest$p.value

        # auprc
        cl0 <- base::subset(dt, binding_class == 0)$TFPred_score # background
        cl1 <- base::subset(dt, binding_class == 1)$TFPred_score # foreground
        cl0 <- cl0[!is.na(cl0)]
        cl1 <- cl1[!is.na(cl1)]
        pr <- PRROC::pr.curve(scores.class0 = cl1, scores.class1=cl0, curve = FALSE)
        auprc <- pr$auc.integral
        
        return(cbind(transcription_factor, tissue, type, pp_auc, auprc, pp_var, low, upp, ttest_pvalue, wiltest_pvalue))
    }, error = function(e){
        print(e)
        return(cbind(transcription_factor, tissue, type, NA, NA, NA, NA, NA, NA, NA))
    })
   
}

results <- foreach::foreach(i=seq_len(nrow(qMatrix)), .combine='rbind', .inorder=T) %dopar% {

    qi <- qMatrix[i, ] |> unlist()
    tfile <- stringr::str_replace_all(opt$path_pattern, purrr::set_names(qi, pp))
    mpath <- stringr::str_replace_all(opt$model_path, purrr::set_names(qi[1:2], pp[1:2]))
    # do something
    res <- tryCatch({
        # count number of files in training peaks directory
        n_training_files <- list.files(file.path(opt$training_peaks_directory, paste0(qi[1], '_', qi[2])), '^peaks_.*.bed$', full.names = T)
        # filter if size == 0
        n_training_files <- length(n_training_files[file.info(n_training_files)$size > 0])
        if(file.exists(tfile)){
            dt <- data.table::fread(tfile)
            if(any(nrow(dt) == 0)){
                evaluations <- data.frame(assay = qi[1], context = qi[2], type = qi[3], auroc = NA, auprc = NA, var = NA, low = NA, upp = NA, ttest_pvalue=NA, wiltest_pvalue = NA, num_training_files = NA, model = NA, path = NA, num_0 = NA, num_1 = NA)
            } else if (any(nrow(dt) > 0)) {
                evaluations <- calculate_logistic_metrics(dt, qi[1], qi[2], qi[3])
                evaluations <- evaluations %>% as.data.frame() %>%
                    stats::setNames(c('assay', 'context', 'type', 'auroc', 'auprc', 'var', 'low', 'upp', 'ttest_pvalue', 'wiltest_pvalue')) %>%
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
                return(evaluations)
            }
        }
    }, error = function(err){
        evaluations <- data.frame(assay = qi[1], context = qi[2], type = qi[3], auroc = NA, auprc = NA, var = NA, low = NA, upp = NA, ttest_pvalue=NA, wiltest_pvalue = NA, num_training_files = NA, model = NA, path = NA, num_0 = NA, num_1 = NA)

        return(evaluations)
    })

    return(res)
}

data.table::fwrite(results, file=opt$statistics_file, sep='\t', row.names=F, col.names=T, quote=F)

# read in the models and extract the weights

models_mtdt <- results %>%
    dplyr::filter(!is.na(model)) %>%
    dplyr::select(model, path) %>%
    dplyr::distinct()

models_weights <- foreach::foreach(i=seq_len(nrow(models_mtdt)), .combine='cbind', .inorder=T) %do% {
    mtdt <- models_mtdt[i, ] |> unlist()
    mpath <- mtdt[2]; model <- mtdt[1]
    model <- readRDS(mpath)
    weights <- coef(model, s = 'lambda.1se') #|> as.data.frame() |> setNames(model)
    return(weights)
} %>% as.matrix() %>% as.data.frame() %>% 
    setNames(models_mtdt$model) %>% 
   dplyr::mutate(feature = c('intercept', paste0('f_', seq_len(nrow(.)-1)))) %>% dplyr::relocate(feature)

data.table::fwrite(models_weights, file=glue("{opt$weights_file_basename}.lambda.1se.txt.gz"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')


models_weights <- foreach::foreach(i=seq_len(nrow(models_mtdt)), .combine='cbind', .inorder=T) %do% {
    mtdt <- models_mtdt[i, ] |> unlist()
    mpath <- mtdt[2]; model <- mtdt[1]
    model <- readRDS(mpath)
    weights <- coef(model, s = 'lambda.min') #|> as.data.frame() |> setNames(model)
    return(weights)
} %>% as.matrix() %>% as.data.frame() %>% 
    setNames(models_mtdt$model) %>% 
    dplyr::mutate(feature = c('intercept', paste0('f_', seq_len(nrow(.)-1)))) %>% dplyr::relocate(feature)

data.table::fwrite(models_weights, file=glue("{opt$weights_file_basename}.lambda.min.txt.gz"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')

doParallel::stopImplicitCluster()