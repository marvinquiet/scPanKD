# conda activate ~/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD"
setwd(project_dir)
data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data"

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

seed = 1993
set.seed(seed)
# --- run ProjecTIL to HNSC
if (opt$experiment == 'ProjecTILs_each_CD8_to_HNSC_CD8') {
    studies = c('EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464')
    for (study in studies) {
        each_result_dir = file.path(ProjecTIL_result_dir, study)
        dir.create(each_result_dir, showWarnings = FALSE)
        res = load_Seurat_data_for_ProjecTIL_to_HNSC(data_dir, ref_study=study, ProjecTIL_data_dir=file.path(data_dir, 'ProjecTILs_CD8T/trimmed'))
        ref = make.reference(ref=res$train_obj, ndim = 20, seed = seed, recalculate.umap = TRUE,
                                    annotation.column = "celltype")
        ref$celltype = as.factor(ref$celltype) # need to be factorized
        query_obj = res$test_obj
        query_obj$celltype_groundtruth = query_obj$celltype
        predicted_test_obj = Run.ProjecTILs(query_obj, ref=ref, labels.col="celltype", skip.normalize=T)
        write.csv(predicted_test_obj@meta.data, file.path(each_result_dir, 'ProjecTIL_predicted_results.csv'), quote=F)
        na_num = sum(is.na(predicted_test_obj$celltype))
        if (sum(is.na(predicted_test_obj$celltype)) > 0) {
            cat('NA exists in prediction:', na_num, '\n')
            if (na_num < round(0.01*dim(predicted_test_obj)[2])) {
                # If the NAs number is not significant, replace it with a randomly predicted results
                na_idx = which(is.na(predicted_test_obj$celltype))
                sampled_celltypes = sample(unique(predicted_test_obj$celltype), na_num)
                predicted_test_obj$celltype[na_idx] = sampled_celltypes
            }
            else {
                cat('Too many NAs in the prediction, we cannot compute metrics!!')
                quit("no")
            }
        }
        # compute results
        acc = calculate_accuracy(predicted_test_obj$celltype_groundtruth, predicted_test_obj$celltype)
        F1 = calculate_macroF1(predicted_test_obj$celltype_groundtruth, predicted_test_obj$celltype)
        ARI = calculate_ARI(predicted_test_obj$celltype_groundtruth, predicted_test_obj$celltype)
        result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
        write.csv(result_df, file.path(each_result_dir, 'ProjecTIL_results_metrics.csv'), quote=F)
        }
} else if (opt$experiment == 'ProjecTILs_each_CD8_cancertrimmed_to_HNSC_CD8') {
    studies = c('HNSCC', 'Melanoma', 'SCC', 'Lung', 'Endometrial', 'Renal', 'Breast')
    for (study in studies) {
        each_result_dir = file.path(ProjecTIL_result_dir, study)
        dir.create(each_result_dir, showWarnings = FALSE)
        res = load_Seurat_data_for_ProjecTIL_to_HNSC(data_dir, ref_study=study, ProjecTIL_data_dir=file.path(data_dir, 'ProjecTILs_CD8T/Tissue_trimmed'))
        ref = make.reference(ref=res$train_obj, ndim = 20, seed = seed, recalculate.umap = TRUE,
                                    annotation.column = "celltype")
        ref$celltype = as.factor(ref$celltype) # need to be factorized
        query_obj = res$test_obj
        query_obj$celltype_groundtruth = query_obj$celltype
        predicted_test_obj = Run.ProjecTILs(query_obj, ref=ref, labels.col="celltype", skip.normalize=T)
        write.csv(predicted_test_obj@meta.data, file.path(each_result_dir, 'ProjecTIL_predicted_results.csv'), quote=F)
        na_num = sum(is.na(predicted_test_obj$celltype))
        if (sum(is.na(predicted_test_obj$celltype)) > 0) {
            cat('NA exists in prediction:', na_num, '\n')
            if (na_num < round(0.01*dim(predicted_test_obj)[2])) {
                # If the NAs number is not significant, replace it with a randomly predicted results
                na_idx = which(is.na(predicted_test_obj$celltype))
                sampled_celltypes = sample(unique(predicted_test_obj$celltype), na_num)
                predicted_test_obj$celltype[na_idx] = sampled_celltypes
            }
            else {
                cat('Too many NAs in the prediction, we cannot compute metrics!!')
                quit("no")
            }
        }
        # compute results
        acc = calculate_accuracy(predicted_test_obj$celltype_groundtruth, predicted_test_obj$celltype)
        F1 = calculate_macroF1(predicted_test_obj$celltype_groundtruth, predicted_test_obj$celltype)
        ARI = calculate_ARI(predicted_test_obj$celltype_groundtruth, predicted_test_obj$celltype)
        result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
        write.csv(result_df, file.path(each_result_dir, 'ProjecTIL_results_metrics.csv'), quote=F)
    }

} else {
    res = load_Seurat_data(data_dir, opt$experiment)
    train_obj = res$train_obj
    query_obj = res$test_obj
    if (opt$experiment == 'Chu_CD4_multibatch_to_Zheng_CD4') {
        common_genes = intersect(rownames(train_obj), rownames(query_obj))
        train_obj = subset(train_obj, features=common_genes)
        query_obj = subset(query_obj, features=common_genes)
    }
    if (!'celltype' %in% colnames(train_obj@meta.data)) {
        train_obj$celltype = train_obj$curated_celltype
    }
    ref = make.reference(ref=train_obj, ndim = 20, seed = 1234, recalculate.umap = TRUE,
                                    annotation.column = "celltype", assay="integrated")
    # if ('multibatch' %in% opt$experiment) {
    #     ref = make.reference(ref=res$train_obj, ndim = 20, seed = 1234, recalculate.umap = TRUE,
    #                                 annotation.column = "curated_celltype", assay="integrated")
    # } else {
    #     ref = make.reference(ref=res$train_obj, ndim = 20, seed = 1234, recalculate.umap = TRUE,
    #                                 annotation.column = "curated_celltype")
    # }
    ref$celltype = as.factor(ref$celltype) # need to be factorized
    
    if (!'celltype' %in% colnames(query_obj@meta.data)) {
        query_obj$celltype = query_obj$curated_celltype
    }
    query_obj$curated_celltype_groundtruth = query_obj$celltype
    predicted_test_obj = Run.ProjecTILs(query_obj, ref=ref, labels.col="celltype", skip.normalize=T)
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
    acc = calculate_accuracy(predicted_test_obj$curated_celltype_groundtruth, predicted_test_obj$celltype)
    F1 = calculate_macroF1(predicted_test_obj$curated_celltype_groundtruth, predicted_test_obj$celltype)
    ARI = calculate_ARI(predicted_test_obj$curated_celltype_groundtruth, predicted_test_obj$celltype)
    result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
    write.csv(result_df, file.path(ProjecTIL_result_dir, 'ProjecTIL_results_metrics.csv'), quote=F)
}
