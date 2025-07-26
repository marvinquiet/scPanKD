# conda activate ~/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD"
setwd(project_dir)
data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data"

library(Seurat)
library(ggplot2)
library(Matrix)
library(optparse)

# source my own library
source("preprocess/Seurat_functions.R")
source("pipelines/calculate_metrics.R")

result_dir = file.path(project_dir, 'results', 'Seurat')
dir.create(result_dir, showWarnings = FALSE)

## === add argument parser
option_list = list(
    make_option('--experiment', type='character', default=NULL,
                help='Load which dataset')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

Seurat_result_dir = file.path(result_dir, opt$experiment)
dir.create(Seurat_result_dir, showWarnings = FALSE)

set.seed(1993)
# --- run ProjecTIL to HNSC
if (opt$experiment == 'ProjecTILs_each_CD8_to_HNSC_CD8') {
    studies = c('EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464')
    for (study in studies) {
        each_result_dir = file.path(Seurat_result_dir, study)
        dir.create(each_result_dir, showWarnings = FALSE)
        res = load_Seurat_data_for_ProjecTIL_to_HNSC(data_dir, ref_study=study)
        predicted_test_obj = predict_Seurat(res$train_obj, res$test_obj, ref_celltype_ind='celltype')
        write.csv(predicted_test_obj@meta.data, file.path(each_result_dir, 'Seurat_predicted_results.csv'), quote=F)
        # compute results
        acc = calculate_accuracy(predicted_test_obj$curated_celltype, predicted_test_obj$predicted.id)
        F1 = calculate_macroF1(predicted_test_obj$curated_celltype, predicted_test_obj$predicted.id)
        ARI = calculate_ARI(predicted_test_obj$curated_celltype, predicted_test_obj$predicted.id)
        result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
        write.csv(result_df, file.path(each_result_dir, 'Seurat_results_metrics.csv'), quote=F)
    }
} else { 
    # --- run other experiments
    res = load_Seurat_data(data_dir, opt$experiment)
    # # if multibatch, then plot integrated plots
    # if (opt$experiment == 'Chu_CD4T_multibatch_validation' || opt$experiment == 'Chu_CD8T_multibatch_validation') {
    #     train_obj = res$train_obj
    #     train_obj = RunUMAP(train_obj, reduction = "pca", dims = 1:30)
    #     p1 = DimPlot(train_obj, reduction = "umap", group.by = "CancerType")
    #     ggsave(file.path(Seurat_result_dir, 'Seurat_integrated_reference_batch.png'))
    #     p2 = DimPlot(train_obj, reduction = "umap", group.by = "curated_celltype", repel = TRUE)
    #     ggsave(file.path(Seurat_result_dir, 'Seurat_integrated_reference_celltype.png'))
    # }
    # if (opt$experiment == 'Chu_GSE179994_CD8T_multibatch_to_HNSC_CD8') {
    #     train_obj = res$train_obj
    #     train_obj = RunUMAP(train_obj, reduction = "pca", dims = 1:30)
    #     p1 = DimPlot(train_obj, reduction = "umap", group.by = "dataset")
    #     ggsave(file.path(Seurat_result_dir, 'Seurat_integrated_reference_batch.png'))
    #     p2 = DimPlot(train_obj, reduction = "umap", group.by = "curated_celltype", repel = TRUE)
    #     ggsave(file.path(Seurat_result_dir, 'Seurat_integrated_reference_celltype.png'))
    # }
    # if ('ProjecTILs_CD8_multibatch' %in% opt$experiment) {
    #     train_obj = res$train_obj
    #     train_obj = RunUMAP(train_obj, reduction = "pca", dims = 1:30)
    #     p1 = DimPlot(train_obj, reduction = "umap", group.by = "Cohort")
    #     ggsave(file.path(Seurat_result_dir, 'Seurat_integrated_reference_batch.png'))
    #     p2 = DimPlot(train_obj, reduction = "umap", group.by = "curated_celltype", repel = TRUE)
    #     ggsave(file.path(Seurat_result_dir, 'Seurat_integrated_reference_celltype.png'))
    # }
    # if (opt$experiment == 'ProjecTILs_CD8_multibatch_trimmed_to_HNSC_CD8') {
    #     predicted_test_obj = predict_Seurat(res$train_obj, res$test_obj, ref_celltype_ind='celltype')
    # }
    # else {
        
    # }
    # ProjecTILs_CD8_multibatch_trimmed_to_HNSC_CD8 / Chu_CD4_multibatch_validation
    predicted_test_obj = predict_Seurat(res$train_obj, res$test_obj, ref_celltype_ind='celltype')

    write.csv(predicted_test_obj@meta.data, file.path(Seurat_result_dir, 'Seurat_predicted_results.csv'), quote=F)
    if (!'curated_celltype' %in% colnames(predicted_test_obj@meta.data)) {
        predicted_test_obj$curated_celltype = predicted_test_obj$celltype
    }
    # compute results
    acc = calculate_accuracy(predicted_test_obj$curated_celltype, predicted_test_obj$predicted.id)
    F1 = calculate_macroF1(predicted_test_obj$curated_celltype, predicted_test_obj$predicted.id)
    ARI = calculate_ARI(predicted_test_obj$curated_celltype, predicted_test_obj$predicted.id)
    result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
    write.csv(result_df, file.path(Seurat_result_dir, 'Seurat_results_metrics.csv'), quote=F)
}
