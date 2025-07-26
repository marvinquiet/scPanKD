# conda activate ~/envs/SpaDOT
.libPaths("/net/mulan/home/wenjinma/Rlib")
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD"
setwd(project_dir)
data_dir = '/net/zootopia/disk1/wenjinma/data/scPanKD_data'

library(aricode)
library(networkD3)
library(dplyr)
library(ggplot2)
library(Matrix)
library(optparse)
source("pipelines/calculate_metrics.R")

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

# --- load ground truth label
if (opt$experiment == 'ProjecTILs_each_CD8_to_HNSC_CD8') {
    metadata = read.csv(file.path(data_dir, 'HNSC_CD8T', 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
    studies = c('EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464')
    for (study in studies) {
        each_result_dir = file.path(method_result_dir, study)
        ## load predicted data
        if (opt$method == 'CellTypist' || opt$method == 'scPanKD') {
            res = read.csv(file.path(method_result_dir, study, paste0(opt$method, '_predicted_celltypes.csv')), row.names=1, header=T)
        } else {
            res = read.csv(file.path(method_result_dir, paste0(opt$method, '_predicted_results.csv')), row.names=1, header=T)
        }
        # --- calculate metrics
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
            write.csv(result_df, file.path(each_result_dir, paste0(opt$method, '_results_metrics.csv')), quote=F)
        }
        if (opt$method == 'scType' || opt$method == 'CellTypist') {
            # pred_res = data.frame(source=metadata$curated_celltype, target=res[, 1]) 
            pred_res = data.frame(source=metadata$curated_celltype, target=res$predicted.id) 
            acc = calculate_accuracy(pred_res[, 1], pred_res[, 2])
            F1 = calculate_macroF1(pred_res[, 1], pred_res[, 2])
            ARI = calculate_ARI(pred_res[, 1], pred_res[, 2])
            result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
            write.csv(result_df, file.path(each_result_dir, paste0(opt$method, '_results_metrics.csv')), quote=F)
        }
        if (opt$method == 'Seurat') {
            pred_res = data.frame(source=metadata$curated_celltype, target=res$predicted.id)
        }
        if (opt$method == 'ProjecTIL') {
            pred_res = data.frame(source=metadata$curated_celltype, target=res$curated_celltype)
        }
    }
} else {
    if (grepl('to_HNSC_CD8', opt$experiment)) {
        metadata = read.csv(file.path(data_dir, 'HNSC_CD8T', 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
    } else if (grepl('to_GSE179994_CD8', opt$experiment)) {
        metadata = read.csv(file.path(data_dir, 'GSE179994_CD8T', 'GSE179994_CD8T_ref_metadata.csv'), header=T, row.names=1)
    } else if (opt$experiment == 'Chu_CD4_multibatch_validation' || grepl('to_Chu_CD4', opt$experiment)) {
        metadata = read.csv(file.path(data_dir, 'Chu_pancancer_CD4T', 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
    } else if (opt$experiment == 'Chu_CD8_multibatch_validation' || grepl('to_Chu_CD8', opt$experiment)) {
        metadata = read.csv(file.path(data_dir, 'Chu_pancancer_CD8T', 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
    } else if (grepl('to_ProjecTILs_CD8', opt$experiment)) {
        metadata = read.csv(file.path(data_dir, 'ProjecTILs_CD8T', 'ProjecTIL_CD8T_metadata.csv'), header=T, row.names=1)
    } else if (grepl('to_Zheng_CD4', opt$experiment)) {
        metadata = read.csv(file.path(data_dir, 'Zheng_CD4T', 'Zheng_CD4T_metadata.csv'), header=T, row.names=1)
    } else {
        stop('Unknown experiment:', opt$experiment)
    }

    ## load predicted data
    if (opt$method == 'CellTypist' || opt$method == 'scPanKD') {
        res = read.csv(file.path(method_result_dir, paste0(opt$method, '_predicted_celltypes.csv')), row.names=1, header=T)
    } else {
        res = read.csv(file.path(method_result_dir, paste0(opt$method, '_predicted_results.csv')), row.names=1, header=T)
    }
    common_cells = intersect(rownames(metadata), rownames(res))
    cat("No. Ground Truth cells:", nrow(metadata), 'No. Predicted cells:', nrow(res), 'No. Common cells:', length(common_cells), '\n')
    res = res[common_cells, ,drop=F]
    metadata = metadata[common_cells, ]
    if ('curated_celltype' %in% colnames(metadata)) {
        metadata$curated_celltype = as.factor(metadata$curated_celltype)
    } else {
        metadata$curated_celltype = as.factor(metadata$celltype)
    }
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
        pred_res = data.frame(source=res$curated_celltype_groundtruth, target=res$celltype)
        na_num = sum(is.na(pred_res$target))
        cat('NA exists in prediction:', na_num, '\n')
        if (na_num < round(0.01*dim(pred_res)[1])) {
            # If the NAs number is not significant, replace it with a randomly predicted results
            set.seed(1234)
            na_idx = which(is.na(pred_res$target))
            sampled_celltypes = sample(unique(pred_res$target), na_num)
            pred_res$target[na_idx] = sampled_celltypes
        }
        else {
            cat('Too many NAs in the prediction, we cannot compute metrics!!')
            quit("no")
        }
        acc = calculate_accuracy(pred_res[, 1], pred_res[, 2])
        F1 = calculate_macroF1(pred_res[, 1], pred_res[, 2])
        ARI = calculate_ARI(pred_res[, 1], pred_res[, 2])
        result_df = data.frame('Acc'=acc, 'macroF1'=F1, 'ARI'=ARI)
        print(result_df)
        write.csv(result_df, file.path(method_result_dir, paste0(opt$method, '_results_metrics.csv')), quote=F)
    }
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
