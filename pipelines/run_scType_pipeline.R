# conda activate Seurat
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# get cell-type-specific gene sets from our in-built database (DB)
#gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

project_dir = "/panfs/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang"
setwd(project_dir)
marker_file = 'data/scType_marker_XYu.xlsx'
workbook = loadWorkbook(marker_file)
sheet_names = getSheetNames(marker_file)
marker_df = read.xlsx(marker_file, sheet = sheet_names[1], rowNames=T)
celltypes = rownames(marker_df)
marker_gs_list = list()
marker_gs_list[['gs_positive']] = list()
marker_gs_list[['gs_negative']] = list()
for (celltype in celltypes) {
    marker_gs_list[['gs_positive']][[celltype]] = unlist(strsplit(marker_df[celltype, 'Positive'], split=','))
    marker_gs_list[['gs_negative']][[celltype]] = unlist(strsplit(marker_df[celltype, 'Negative'], split=','))
}

library(ggplot2)
library(Matrix)
library(optparse)

# source my own library
source("preprocess/Seurat_functions.R")
source("pipelines/calculate_metrics.R")

result_dir = file.path(project_dir, 'results', 'scType')
dir.create(result_dir, showWarnings = FALSE)

## === add argument parser
option_list = list(
    make_option('--experiment', type='character', default=NULL,
                help='Load which dataset')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

scType_result_dir = file.path(result_dir, opt$experiment)
dir.create(scType_result_dir, showWarnings = FALSE)

if (opt$experiment == 'ProjecTILs_each_CD8_to_HNSC_CD8') {
    # since we only have one test dataset, the results will be the same
    studies = c('EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464')
    res = load_Seurat_data_for_ProjecTIL_to_HNSC(project_dir, ref_study=studies[1])
    test_obj = res$test_obj
    test_obj = ScaleData(test_obj, features=rownames(test_obj))
    test_scaled_data = test_obj@assays$RNA@scale.data
    scType_res = sctype_score(scRNAseqData = test_scaled_data, scaled = TRUE, gs = marker_gs_list$gs_positive, gs2 = marker_gs_list$gs_negative)
    write.csv(t(scType_res), file.path(scType_result_dir, 'scType_predicted_scores.csv'), quote=F)
    scType_predict = apply(t(scType_res), 1, which.max)
    scType_predict_celltypes = colnames(t(scType_res))[unname(scType_predict)]
    names(scType_predict_celltypes) = names(scType_predict)
    write.csv(as.data.frame(scType_predict_celltypes), file.path(scType_result_dir, 'scType_predicted_results.csv'), quote=F)
} else {
    res = load_Seurat_data(opt$experiment, project_dir)
    test_obj = res$test_obj
    test_obj = ScaleData(test_obj, features=rownames(test_obj))
    test_scaled_data = test_obj@assays$RNA@scale.data
    scType_res = sctype_score(scRNAseqData = test_scaled_data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    write.csv(t(scType_res), file.path(scType_result_dir, 'scType_predicted_scores.csv'), quote=F)
    scType_predict = apply(t(scType_res), 1, which.max)
    scType_predict_celltypes = colnames(t(scType_res))[unname(scType_predict)]
    names(scType_predict_celltypes) = names(scType_predict)
    write.csv(as.data.frame(scType_predict_celltypes), file.path(scType_result_dir, 'scType_predicted_results.csv'), quote=F)
}


