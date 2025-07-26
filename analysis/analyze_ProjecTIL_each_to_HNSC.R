# conda activate ~/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD"
setwd(project_dir)
data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data"

library(ggpubr)
library(ggplot2)
library(cowplot)
library(dplyr)
source("/net/mulan/home/wenjinma/projects/scPanKD/pipelines/calculate_metrics.R")

experiment = "ProjecTILs_each_CD8_to_HNSC_CD8"
studies = c('EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464')
methods = c('scPanKD', 'Seurat', 'ProjecTIL', 'CellTypist')
metrics = c('Acc', 'macroF1', 'ARI')

metadata = read.csv(file.path(data_dir, 'HNSC_CD8T', 'HNSC_CD8T_metadata.csv'),
                         header=T, row.names=1)
metric_df = data.frame(matrix(ncol = 4, nrow = 0))
colnames(metric_df) = c('study', 'method', 'metric', 'value')
for (method in methods) {
    result_dir = file.path(project_dir, 'results', method, experiment)
    for (study in studies) {
        method_result_dir = file.path(result_dir, study)
        if (method == 'CellTypist') {
            res = read.csv(file.path(method_result_dir, 'predicted_labels.csv'), row.names=1, header=T)
        } else if (method == 'scPanKD') {
            res = read.csv(file.path(method_result_dir, paste0(method, '_predicted_celltypes.csv')), row.names=1, header=T)
        } else {
            res = read.csv(file.path(method_result_dir, paste0(method, '_predicted_results.csv')), row.names=1, header=T)
        }
        common_cells = intersect(rownames(metadata), rownames(res))
        cat("No. Ground Truth cells:", nrow(metadata), 'No. Predicted cells:', nrow(res), 'No. Common cells:', length(common_cells), '\n')
        if (ncol(res) == 1) {
            res = res[common_cells, ,drop=F]
        } else {
            res = res[common_cells, ]
        } 
        res_metadata = metadata[common_cells, ]
        if (method == 'CellTypist') {
            pred_res = data.frame(source=res_metadata$curated_celltype, target=res[, 1]) 
        }
        if (method == 'scPanKD') {
            pred_res = data.frame(source=res_metadata$curated_celltype, target=res$pred_celltype)
        }
        if (method == 'Seurat') {
            pred_res = data.frame(source=res_metadata$curated_celltype, target=res$predicted.id)
        }
        if (method == 'ProjecTIL') {
            pred_res = data.frame(source=res_metadata$curated_celltype, target=res$celltype)
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
    # Separate scType data from other methods
    scType_data = sub_metric_df[sub_metric_df$method == 'scType', ]
    boxplot_data = sub_metric_df[sub_metric_df$method != 'scType', ]
    # calculate group mean
    group_means = boxplot_data %>%
        group_by(method) %>%
        summarize(mean_value = mean(value)) %>%
        arrange(desc(mean_value))
    boxplot_data$method = factor(boxplot_data$method, levels = group_means$method)
    p = ggboxplot(boxplot_data, x='method', y='value', color='method') + 
        stat_summary(fun.y=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
        scale_color_brewer(palette = 'Set1') + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.title =element_blank())
        # Add scType as horizontal line with better positioning
    legend = get_plot_component(p, 'guide-box', return_all = TRUE)[[4]]
    ggsave(file.path(project_dir, 'results', 'shared_legend.pdf'), 
       plot = legend, width = 6, height = 1)
    if (nrow(scType_data) > 0) {
        scType_value = scType_data$value[1]
        y_range = layer_scales(p)$y$range$range
        p = p + 
            geom_hline(yintercept = scType_value, color = "blue", linetype = "dashed", size = 1) +
            # Position label at the right edge
            annotate("text", x = Inf, y = -Inf, 
                     label = paste0("scType: ", round(scType_value, 3)), 
                     hjust = 1.1, vjust = -0.5, color = "blue", size = 3.5, fontface = "bold") +
            # Extend plot area slightly to accommodate label
            coord_cartesian(clip = "off") +
            theme(plot.margin = margin(5.5, 40, 5.5, 5.5, "pt"))
    }
    p = p + theme(legend.position = "none")
    ggsave(file.path(project_dir, 'results', paste0(metric, '_summarized_results.pdf')), width=3, height=3)
}



