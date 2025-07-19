# conda activate ~/envs/SpaDOT
.libPaths("/net/mulan/home/wenjinma/Rlib")
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD"
setwd(project_dir)

library(aricode)
library(networkD3)
library(dplyr)
library(ggplot2)
library(Matrix)
library(optparse)
source("/net/mulan/home/wenjinma/projects/scPanKD/pipelines/calculate_metrics.R")

## === add argument parser
option_list = list(
    make_option('--experiment', type='character', default=NULL,
                help='Load which dataset'),
    make_option('--method', type='character', default=NULL,
                help='Load which dataset')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

result_dir = file.path(project_dir, 'results', opt$method)
method_result_dir = file.path(result_dir, opt$experiment)
data_dir = '/net/zootopia/disk1/wenjinma/data/scPanKD_data'
## load ground truth label
if (opt$experiment == 'Chu_CD4T_validation' || opt$experiment == 'Zheng_CD4_to_Chu_CD4' || 
    opt$experiment == 'Chu_CD4T_multibatch_validation' || opt$experiment == 'Zheng_CD4_multibatch_to_Chu_CD4') {
    metadata = read.csv(file.path(data_dir, 'Chu_pancancer_CD4T', 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
}
if (opt$experiment == 'Chu_CD8T_validation' || opt$experiment == 'GSE179994_CD8_to_Chu_CD8T' || opt$experiment == 'ProjecTILs_CD8_to_Chu_CD8T' ||
    opt$experiment == 'Chu_CD8T_multibatch_validation' || opt$experiment == 'ProjecTILs_CD8_multibatch_to_Chu_CD8' || opt$experiment == 'GSE179994_CD8T_multibatch_to_Chu_CD8') {
    metadata = read.csv(file.path(data_dir, 'Chu_pancancer_CD8T', 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
}
if (opt$experiment == 'GSE179994_CD8_to_HNSC_CD8' || opt$experiment == 'ProjecTILs_CD8_to_HNSC_CD8' ||
    opt$experiment == 'GSE179994_CD8T_multibatch_to_HNSC_CD8' || opt$experiment == 'ProjecTILs_CD8_multibatch_to_HNSC_CD8') {
    metadata = read.csv(file.path(data_dir, 'HNSC_CD8T', 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
    # metadata$curated_celltype = metadata$Label
}

## load predicted data
if (opt$method == 'CellTypist' || opt$method == 'Cellcano' || opt$method == 'scPanKD') {
    res = read.csv(file.path(method_result_dir, paste0(opt$method, '_predicted_celltypes.csv')), row.names=1, header=T)
} else {
    res = read.csv(file.path(method_result_dir, paste0(opt$method, '_predicted_results.csv')), row.names=1, header=T)
}
common_cells = intersect(rownames(metadata), rownames(res))
cat("No. Ground Truth cells:", nrow(metadata), 'No. Predicted cells:', nrow(res), 'No. Common cells:', length(common_cells), '\n')
res = res[common_cells, ,drop=F]
metadata = metadata[common_cells, ]
if (opt$method == 'Cellcano' || opt$method == 'scPanKD') {
    pred_res = data.frame(source=metadata$curated_celltype, target=res$pred_celltype)
    acc = calculate_accuracy(pred_res[, 1], pred_res[, 2])
    F1 = calculate_macroF1(pred_res[, 1], pred_res[, 2])
    ARI = calculate_ARI(pred_res[, 1], pred_res[, 2])
    result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
    pred_res = data.frame(source=metadata$curated_celltype, target=res$firstround_pred_celltype)
    acc = calculate_accuracy(pred_res[, 1], pred_res[, 2])
    F1 = calculate_macroF1(pred_res[, 1], pred_res[, 2])
    ARI = calculate_ARI(pred_res[, 1], pred_res[, 2])
    result_df = rbind(result_df, data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI))
    rownames(result_df) = c('Second round', 'First round')
    print(result_df)
    write.csv(result_df, file.path(method_result_dir, paste0(opt$method, '_results_metrics.csv')), quote=F)
}

if (opt$method == 'scType' || opt$method == 'CellTypist') {
    # pred_res = data.frame(source=metadata$curated_celltype, target=res[, 1]) 
    pred_res = data.frame(source=metadata$curated_celltype, target=res$predicted.id) 
    acc = calculate_accuracy(pred_res[, 1], pred_res[, 2])
    F1 = calculate_macroF1(pred_res[, 1], pred_res[, 2])
    ARI = calculate_ARI(pred_res[, 1], pred_res[, 2])
    result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
    write.csv(result_df, file.path(method_result_dir, paste0(opt$method, '_results_metrics.csv')), quote=F)
}

if (opt$method == 'Seurat') {
    pred_res = data.frame(source=metadata$curated_celltype, target=res$predicted.id)
}
if (opt$method == 'ProjecTIL') {
    pred_res = data.frame(source=metadata$curated_celltype, target=res$curated_celltype)
}

## Create Sankey link
# links = pred_res %>%
#       group_by_all() %>%
#       summarise(value = n())
# links = as.data.frame(links)
# links$source = paste0('source_', links$source)
# links$target = paste0('target_', links$target)
# nodes = data.frame(name=c(as.character(links$source), 
#                            as.character(links$target)) %>% unique()
#                       )
# links$IDsource <- match(links$source, nodes$name)-1 
# links$IDtarget <- match(links$target, nodes$name)-1

# # Create Sankey plot
# p = sankeyNetwork(Links = links, Nodes = nodes,
#                 Source = "IDsource", Target = "IDtarget",
#                 Value = "value", NodeID = "name", 
#                 sinksRight=FALSE, units = "TWh", fontSize = 12, nodeWidth = 30)
# require(htmlwidgets)
# saveWidget(p, file=file.path(method_result_dir, paste0(opt$method, '_sankey_plot.html')))
# require(webshot)
# webshot(file.path(method_result_dir, paste0(opt$method, '_sankey_plot.html')), file.path(method_result_dir, paste0(opt$method, '_sankey_plot.pdf')))
