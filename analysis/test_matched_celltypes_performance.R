# conda activate ~/envs/SpaDOT
.libPaths("/net/mulan/home/wenjinma/Rlib")
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD"
setwd(project_dir)

library(aricode)
library(networkD3)
library(dplyr)
library(ggplot2)
library(Matrix)
source("/net/mulan/home/wenjinma/projects/scPanKD/pipelines/calculate_metrics.R")


data_dir = '/net/zootopia/disk1/wenjinma/data/scPanKD_data'
# --- ProjecTILs_CD8_multibatch_to_Chu_CD8, Exp7
true_df = read.csv(file.path(data_dir, 'Chu_pancancer_CD8T', 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
scPanKD_pred_df = read.csv(file.path(project_dir, '/results/scPanKD/ProjecTILs_CD8_multibatch_to_Chu_CD8/scPanKD_predicted_celltypes.csv'), header=T, row.names=1)
Seurat_multibatch_pred_df = read.csv(file.path(project_dir, 'results/Seurat/ProjecTILs_CD8_multibatch_to_Chu_CD8T/Seurat_predicted_results.csv'), header=T, row.names=1)
Seurat_pred_df = read.csv(file.path(project_dir, 'results/Seurat/ProjecTILs_CD8_to_Chu_CD8T/Seurat_predicted_results.csv'), header=T, row.names=1)

explore_true_vs_pred = function(true_df, pred_df, pred_col) {
    cat('Cell names are same:', all(rownames(true_df)==rownames(pred_df)), '\n')
    true_df = true_df[true_df$curated_celltype %in% c('CD8.Naivelike', 'CD8.CM', 'CD8.Exh'), ] # 4399
    cat('GT number of cells:', nrow(true_df), '\n')
    pred_df = pred_df[pred_df[[pred_col]] %in% c('CD8.NaiveLike', 'CD8.CM', 'CD8.TEX'), ] # 4882
    cat('Pred number of cells:', nrow(pred_df), '\n')
    common_cells = intersect(rownames(true_df), rownames(pred_df)) # 3548
    cat('Common number of cells:', length(common_cells), '\n')
    filtered_true_df = true_df[common_cells, ]
    filtered_pred_df = pred_df[common_cells, ]
    cat('ARI:', calculate_ARI(filtered_true_df$curated_celltype, filtered_pred_df[[pred_col]]), '\n') # 0.914
}

explore_true_vs_pred(true_df, scPanKD_pred_df, pred_col='pred_celltype')
explore_true_vs_pred(true_df, scPanKD_pred_df, pred_col='firstround_pred_celltype')
explore_true_vs_pred(true_df, Seurat_multibatch_pred_df, pred_col='predicted.id')
explore_true_vs_pred(true_df, Seurat_pred_df, pred_col='predicted.id')


# --- ProjecTILs_CD8_to_HNSC_CD8, Exp5
true_df = read.csv(file.path(data_dir, 'HNSC_CD8T', 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
scPanKD_pred_df = read.csv(file.path(project_dir, 'results/scPanKD/ProjecTILs_CD8_multibatch_to_HNSC_CD8/scPanKD_predicted_celltypes.csv'), header=T, row.names=1)
Seurat_multibatch_pred_df = read.csv(file.path(project_dir, 'results/Seurat/ProjecTILs_CD8_multibatch_to_HNSC_CD8/Seurat_predicted_results.csv'), header=T, row.names=1)
Seurat_pred_df = read.csv(file.path(project_dir, 'results/Seurat/ProjecTILs_CD8_to_HNSC_CD8/Seurat_predicted_results.csv'), header=T, row.names=1)

explore_true_vs_pred = function(true_df, pred_df, pred_col) {
    cat('Cell names are same:', all(rownames(true_df)==rownames(pred_df)), '\n')
    cat('GT number of cells:', nrow(true_df), '\n')
    pred_df = pred_df[pred_df[[pred_col]] %in% c('CD8.NaiveLike', 'CD8.EM', 'CD8.TEMRA', 'CD8.TEX'), ] # 4882
    cat('Pred number of cells:', nrow(pred_df), '\n')
    common_cells = intersect(rownames(true_df), rownames(pred_df)) # 3548
    cat('Common number of cells:', length(common_cells), '\n')
    filtered_true_df = true_df[common_cells, ]
    filtered_pred_df = pred_df[common_cells, ]
    cat('ARI:', calculate_ARI(filtered_true_df$curated_celltype, filtered_pred_df[[pred_col]]), '\n') # 0.914
    excluded_acc = calculate_accuracy(filtered_true_df$curated_celltype, filtered_pred_df[[pred_col]])
    excluded_cov = length(common_cells)/nrow(true_df)
    excluded_f1 = 2*excluded_acc*excluded_cov/(excluded_acc+excluded_cov)
    cat('Excluded Accuracy:', excluded_acc, '\n') # 0.914
    cat('Excluded Coverage:', excluded_cov, '\n')
    cat('Excluded F1:', excluded_f1, '\n')

}

explore_true_vs_pred(true_df, scPanKD_pred_df, pred_col='pred_celltype')
explore_true_vs_pred(true_df, scPanKD_pred_df, pred_col='firstround_pred_celltype')
explore_true_vs_pred(true_df, Seurat_multibatch_pred_df, pred_col='predicted.id')
explore_true_vs_pred(true_df, Seurat_pred_df, pred_col='predicted.id')
