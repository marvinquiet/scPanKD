# conda activate ~/envs/SpaDOT
.libPaths("/net/mulan/home/wenjinma/Rlib")
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD"
setwd(project_dir)

library(ggpubr)
library(ggplot2)
library(dplyr)
source("/net/mulan/home/wenjinma/projects/scPanKD/pipelines/calculate_metrics.R")

experiment = "ProjecTILs_each_CD8_to_HNSC_CD8"
studies = c('EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464')
methods = c('Cellcano', 'Seurat', 'ProjecTIL', 'CellTypist', 'scType')
metrics = c('Acc', 'macroF1', 'ARI')

metadata = read.csv(file.path(project_dir, 'data', 'HNSC_CD8T', 'HNSC_CD8T_metadata.csv'),
                         header=T, row.names=1)
metric_df = data.frame(matrix(ncol = 4, nrow = 0))
colnames(metric_df) = c('study', 'method', 'metric', 'value')
for (method in methods) {
    result_dir = file.path(project_dir, 'results', method, experiment)
    for (study in studies) {
        method_result_dir = file.path(result_dir, study)
        if (method == 'CellTypist') {
            res = read.csv(file.path(method_result_dir, 'predicted_labels.csv'), row.names=1, header=T)
        } else if (method == 'Cellcano') {
            res = read.csv(file.path(method_result_dir, paste0(method, '_predicted_celltypes.csv')), row.names=1, header=T)
        } else {
            res = read.csv(file.path(method_result_dir, paste0(method, '_predicted_results.csv')), row.names=1, header=T)
        }
        common_cells = intersect(rownames(metadata), rownames(res))
        cat("No. Ground Truth cells:", nrow(metadata), 'No. Predicted cells:', nrow(res), 'No. Common cells:', length(common_cells), '\n')
        res = res[common_cells, ,drop=F]
        res_metadata = metadata[common_cells, ]
        if (method == 'scType' || method == 'CellTypist') {
            pred_res = data.frame(source=res_metadata$curated_celltype, target=res[, 1]) 
        }
        if (method == 'Cellcano') {
            pred_res = data.frame(source=res_metadata$curated_celltype, target=res$pred_celltype)
        }
        if (method == 'Seurat') {
            pred_res = data.frame(source=res_metadata$curated_celltype, target=res$predicted.id)
        }
        if (method == 'ProjecTIL') {
            pred_res = data.frame(source=res_metadata$curated_celltype, target=res$curated_celltype)
        }
        acc = calculate_accuracy(pred_res$source, pred_res$target)
        F1 = calculate_macroF1(pred_res$source, pred_res$target)
        ARI = calculate_ARI(pred_res$source, pred_res$target)
        metric_df[nrow(metric_df)+1, ] = c(study, method, 'Acc', acc)
        metric_df[nrow(metric_df)+1, ] = c(study, method, 'macroF1', F1)
        metric_df[nrow(metric_df)+1, ] = c(study, method, 'ARI', ARI)
    }
}
# add scType results -> this only has one result, so I separate it out
method = "scType"
result_dir = file.path(project_dir, 'results', method, experiment)
res = read.csv(file.path(result_dir, paste0(method, '_predicted_results.csv')), row.names=1, header=T)
common_cells = intersect(rownames(metadata), rownames(res))
cat("No. Ground Truth cells:", nrow(metadata), 'No. Predicted cells:', nrow(res), 'No. Common cells:', length(common_cells), '\n')
res = res[common_cells, ,drop=F]
res_metadata = res_metadata[common_cells, ]
pred_res = data.frame(source=res_metadata$curated_celltype, target=res[, 1]) 
acc = calculate_accuracy(pred_res$source, pred_res$target)
F1 = calculate_macroF1(pred_res$source, pred_res$target)
ARI = calculate_ARI(pred_res$source, pred_res$target)
metric_df[nrow(metric_df)+1, ] = c('nostudy', method, 'Acc', acc)
metric_df[nrow(metric_df)+1, ] = c('nostudy', method, 'macroF1', F1)
metric_df[nrow(metric_df)+1, ] = c('nostudy', method, 'ARI', ARI)
write.csv(metric_df, file.path(project_dir, 'results', 'summarized_results.csv'), quote=F)

# --- plot a boxplot
metric_df$value = as.numeric(metric_df$value)
for (metric in metrics) {
    sub_metric_df = metric_df[metric_df$metric == metric, ]
    # calculate group mean
    group_means = sub_metric_df %>%
        group_by(method) %>%
        summarize(mean_value = mean(value)) %>%
        arrange(desc(mean_value))
    sub_metric_df$method = factor(sub_metric_df$method, levels = group_means$method)
    ggboxplot(sub_metric_df, x='method', y='value') + 
        stat_summary(fun.y=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.title =element_blank())
    ggsave(file.path(project_dir, 'results', paste0(metric, '_summarized_results.png')), width=3, height=3)
}


