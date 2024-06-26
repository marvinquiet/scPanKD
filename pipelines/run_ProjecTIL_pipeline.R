# conda activate Seurat
project_dir = "/panfs/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang"
setwd(project_dir)

library(Seurat)
library(ProjecTILs)
library(ggplot2)
library(Matrix)
library(optparse)

# source my own library
source("preprocess/Seurat_functions.R")
source("pipelines/calculate_metrics.R")

result_dir = file.path(project_dir, 'results', 'ProjecTIL')
dir.create(result_dir, showWarnings = FALSE)

## === add argument parser
option_list = list(
    make_option('--experiment', type='character', default=NULL,
                help='Load which dataset')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

ProjecTIL_result_dir = file.path(result_dir, opt$experiment)
dir.create(ProjecTIL_result_dir, showWarnings = FALSE)

# --- run ProjecTIL to HNSC
if (opt$experiment == 'ProjecTILs_each_CD8_to_HNSC_CD8') {
    studies = c('EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464')
    for (study in studies) {
        each_result_dir = file.path(ProjecTIL_result_dir, study)
        dir.create(each_result_dir, showWarnings = FALSE)
        res = load_Seurat_data_for_ProjecTIL_to_HNSC(project_dir, ref_study=study)
        ref = make.reference(ref=res$train_obj, ndim = 20, seed = 1234, recalculate.umap = TRUE,
                                    annotation.column = "curated_celltype")
        ref$curated_celltype = as.factor(ref$curated_celltype) # need to be factorized
        query_obj = res$test_obj
        query_obj$curated_celltype_groundtruth = query_obj$curated_celltype
        predicted_test_obj = Run.ProjecTILs(query_obj, ref=ref, labels.col="curated_celltype", skip.normalize=T)
        write.csv(predicted_test_obj@meta.data, file.path(each_result_dir, 'ProjecTIL_predicted_results.csv'), quote=F)
        na_num = sum(is.na(predicted_test_obj$curated_celltype))
        if (sum(is.na(predicted_test_obj$curated_celltype)) > 0) {
            cat('NA exists in prediction:', na_num, '\n')
            if (na_num < round(0.01*dim(predicted_test_obj)[2])) {
                # If the NAs number is not significant, replace it with a randomly predicted results
                set.seed(1234)
                na_idx = which(is.na(predicted_test_obj$curated_celltype))
                sampled_celltypes = sample(unique(predicted_test_obj$curated_celltype), na_num)
                predicted_test_obj$curated_celltype[na_idx] = sampled_celltypes
            }
            else {
                cat('Too many NAs in the prediction, we cannot compute metrics!!')
                quit("no")
            }
        }
        # compute results
        acc = calculate_accuracy(predicted_test_obj$curated_celltype_groundtruth, predicted_test_obj$curated_celltype)
        F1 = calculate_macroF1(predicted_test_obj$curated_celltype_groundtruth, predicted_test_obj$curated_celltype)
        ARI = calculate_ARI(predicted_test_obj$curated_celltype_groundtruth, predicted_test_obj$curated_celltype)
        result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
        write.csv(result_df, file.path(each_result_dir, 'ProjecTIL_results_metrics.csv'), quote=F)
        }
}  else {
    res = load_Seurat_data(opt$experiment, project_dir)
    if ('multibatch' %in% opt$experiment) {
        ref = make.reference(ref=res$train_obj, ndim = 20, seed = 1234, recalculate.umap = TRUE,
                                    annotation.column = "curated_celltype", assay="integrated")
    } else {
        ref = make.reference(ref=res$train_obj, ndim = 20, seed = 1234, recalculate.umap = TRUE,
                                    annotation.column = "curated_celltype")
    }
    ref$curated_celltype = as.factor(ref$curated_celltype) # need to be factorized
    query_obj = res$test_obj
    query_obj$curated_celltype_groundtruth = query_obj$curated_celltype
    predicted_test_obj = Run.ProjecTILs(query_obj, ref=ref, labels.col="curated_celltype", skip.normalize=T)
    write.csv(predicted_test_obj@meta.data, file.path(ProjecTIL_result_dir, 'ProjecTIL_predicted_results.csv'), quote=F)
    na_num = sum(is.na(predicted_test_obj$curated_celltype))
    if (sum(is.na(predicted_test_obj$curated_celltype)) > 0) {
        cat('NA exists in prediction:', na_num, '\n')
        if (na_num < round(0.01*dim(predicted_test_obj)[2])) {
            # If the NAs number is not significant, replace it with a randomly predicted results
            set.seed(1234)
            na_idx = which(is.na(predicted_test_obj$curated_celltype))
            sampled_celltypes = sample(unique(predicted_test_obj$curated_celltype), na_num)
            predicted_test_obj$curated_celltype[na_idx] = sampled_celltypes
        }
        else {
            cat('Too many NAs in the prediction, we cannot compute metrics!!')
            quit("no")
        }
    }
    # compute results
    acc = calculate_accuracy(predicted_test_obj$curated_celltype_groundtruth, predicted_test_obj$curated_celltype)
    F1 = calculate_macroF1(predicted_test_obj$curated_celltype_groundtruth, predicted_test_obj$curated_celltype)
    ARI = calculate_ARI(predicted_test_obj$curated_celltype_groundtruth, predicted_test_obj$curated_celltype)
    result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
    write.csv(result_df, file.path(ProjecTIL_result_dir, 'ProjecTIL_results_metrics.csv'), quote=F)
}
