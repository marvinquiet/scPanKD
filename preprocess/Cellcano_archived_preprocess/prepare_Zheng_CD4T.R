.libPaths("/net/mulan/home/wenjinma/Rlib")
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)

split_by = '80:20' # 80:20 | tissue | 80_enlarge:20

project_dir = "/net/mulan/home/wenjinma/minicluster/projects/test_Cellcano_XuefengWang"
result_dir = file.path(project_dir, 'ZhengCD4T_results')
dir.create(result_dir, showWarnings=F)

data_obj = readRDS(file.path(project_dir, 'Zheng.PancanerT', 'CD4.rds'))
data_obj$curated_celltype = data_obj$Annotation

# filter genes
integrated_counts = GetAssayData(data_obj, slot="data", assay="integrated")   
selected_genes = rownames(integrated_counts)
counts = GetAssayData(data_obj, slot="counts", assay="RNA")
counts = counts[selected_genes,] ## 800 genes x 12631 cells

pancancerT_obj = CreateSeuratObject(counts, meta.data=data_obj@meta.data)

# perform data analysis
pancancerT_obj = NormalizeData(pancancerT_obj)
pancancerT_obj = FindVariableFeatures(pancancerT_obj, selection.method = "vst", nfeatures = 2000)
pancancerT_obj = ScaleData(pancancerT_obj, features = rownames(pancancerT_obj))
pancancerT_obj = RunPCA(pancancerT_obj, features = VariableFeatures(object = pancancerT_obj))
pancancerT_obj = FindNeighbors(pancancerT_obj, dims=1:20)
pancancerT_obj = FindClusters(pancancerT_obj, resolution=0.5)
pancancerT_obj = RunUMAP(pancancerT_obj, dims=1:20)
gg1 = DimPlot(pancancerT_obj, reduction='pca', group.by='curated_celltype', label=T)
gg2 = DimPlot(pancancerT_obj, reduction='pca', group.by='cancerType', label=T)
gg1 + gg2
ggsave(file.path(result_dir, 'PCA_analysis.png'), width=10, height=5)
gg1 = DimPlot(pancancerT_obj, reduction='umap', group.by='curated_celltype', label=T)
gg2 = DimPlot(pancancerT_obj, reduction='umap', group.by='cancerType', label=T)
gg1 + gg2
ggsave(file.path(result_dir, 'UMAP_analysis.png'), width=10, height=5)


cal_macroF1 <- function(true_celltypes, pred_celltypes) {
    ## calculate macroF1
    union_celltypes = union(true_celltypes, pred_celltypes)
    cm = matrix(0, nrow=length(union_celltypes), ncol=length(union_celltypes))
    rownames(cm) = union_celltypes
    colnames(cm) = union_celltypes

    ## correct levels of target object
    predict_tb = table(true_celltypes, pred_celltypes)
    cm[rownames(predict_tb), colnames(predict_tb)] = predict_tb

    diag = diag(cm)
    precision = diag / colSums(cm) 
    precision[is.nan(precision)] = 0
    recall = diag / rowSums(cm) 
    recall[is.nan(recall)] = 0
    f1 = 2 * precision * recall / (precision + recall) 
    f1[is.nan(f1)] = 0
    macroF1 = mean(f1)
    return(macroF1)
}



# run data split
if (split_by == 'tissue') { # split by cancer tissue
    # random sample
    set.seed(1234)
    train_tissues = sample(selected_tissues, 5)
    test_tissues = sample(setdiff(selected_tissues, train_tissues), 5)
    train_obj = subset(pancancerT_obj, subset=CancerType %in% train_tissues)
    test_obj = subset(pancancerT_obj, subset=CancerType %in% test_tissues)
    writeMM(train_obj@assays$RNA@counts, file.path(result_dir, 'train_pancanT.mtx'))
    write(rownames(train_obj@assays$RNA@counts), file.path(result_dir, 'train_pancanT_genes.tsv'))
    write(colnames(train_obj@assays$RNA@counts), file.path(result_dir, 'train_pancanT_barcodes.tsv'))
    write.csv(train_obj@meta.data, file.path(result_dir, 'train_pancanT_metadata.csv'), quote=F)
    ## run Cellcano
    #Cellcano train -i train_pancanT -m train_pancanT_metadata.csv --prefix pancanT_subtypes_trained
    #---change the column names of metadata
    #Cellcano train -i train_pancanT -m train_pancanT_metadata.csv --prefix pancanT_combinedtypes_trained
    writeMM(test_obj@assays$RNA@counts, file.path(result_dir, 'test_pancanT.mtx'))
    write(rownames(test_obj@assays$RNA@counts), file.path(result_dir, 'test_pancanT_genes.tsv'))
    write(colnames(test_obj@assays$RNA@counts), file.path(result_dir, 'test_pancanT_barcodes.tsv'))
    write.csv(test_obj@meta.data, file.path(result_dir, 'test_pancanT_metadata.csv'), quote=F)
    ## run Cellcano predict
    #Cellcano predict -i test_pancanT --trained_model output/pancanT_subtypes_trainedMLP_model -o output --prefix predict_pancanT_subtypes
    #Cellcano predict -i test_pancanT --trained_model output/pancanT_combinedtypes_trainedMLP_model -o output --prefix predict_pancanT_combinedtypes
    ## analyze results
    # subcelltype
    tgt_metadata = read.csv(file.path(result_dir, 'test_pancanT_metadata.csv'), header=T, row.names=1)
    subcelltype_pred_result = read.csv(file.path(result_dir, 'output', 'predict_pancanT_subtypescelltypes.csv'), header=T, row.names=1)
    sum(tgt_metadata$cell.type == subcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype']) / nrow(tgt_metadata)
    sum(tgt_metadata$cell.type == subcelltype_pred_result[rownames(tgt_metadata), 'firstround_pred_celltype']) / nrow(tgt_metadata)
    table(tgt_metadata$cell.type, subcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype'])
    # combined cell type
    combinedcelltype_pred_result = read.csv(file.path(result_dir, 'output', 'predict_pancanT_combinedtypescelltypes.csv'), header=T, row.names=1)
    sum(tgt_metadata$curated_celltype == combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype']) / nrow(tgt_metadata)
    sum(tgt_metadata$curated_celltype == combinedcelltype_pred_result[rownames(tgt_metadata), 'firstround_pred_celltype']) / nrow(tgt_metadata)
    table(tgt_metadata$curated_celltype, combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype'])
}

if (split_by == "80:20") { # split by proportion
    set.seed(1234)
    train_idx = sample(1:ncol(pancancerT_obj), round(0.8*ncol(pancancerT_obj)))
    test_idx = setdiff(1:ncol(pancancerT_obj), train_idx)
    train_obj = pancancerT_obj[, train_idx]
    test_obj = pancancerT_obj[, test_idx]
    writeMM(train_obj@assays$RNA@counts, file.path(result_dir, 'train_80_pancanT.mtx'))
    write(rownames(train_obj@assays$RNA@counts), file.path(result_dir, 'train_80_pancanT_genes.tsv'))
    write(colnames(train_obj@assays$RNA@counts), file.path(result_dir, 'train_80_pancanT_barcodes.tsv'))
    write.csv(train_obj@meta.data, file.path(result_dir, 'train_80_pancanT_metadata.csv'), quote=F)
    #---change the column names of metadata and run Cellcano
    #Cellcano train -i train_80_pancanT -m train_80_pancanT_metadata.csv --prefix pancanT_80_combinedtypes_trained
    writeMM(test_obj@assays$RNA@counts, file.path(result_dir, 'test_20_pancanT.mtx'))
    write(rownames(test_obj@assays$RNA@counts), file.path(result_dir, 'test_20_pancanT_genes.tsv'))
    write(colnames(test_obj@assays$RNA@counts), file.path(result_dir, 'test_20_pancanT_barcodes.tsv'))
    write.csv(test_obj@meta.data, file.path(result_dir, 'test_20_pancanT_metadata.csv'), quote=F)
    ## run Cellcano predict
    #Cellcano predict -i test_20_pancanT --trained_model output/pancanT_80_combinedtypes_trainedMLP_model -o output --prefix predict_20_pancanT_combinedtypes
    ## analyze results
    tgt_metadata = read.csv(file.path(result_dir, 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
    combinedcelltype_pred_result = read.csv(file.path(result_dir, 'output', 'predict_20_pancanT_combinedtypescelltypes.csv'), header=T, row.names=1)
    sum(tgt_metadata$curated_celltype == combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype']) / nrow(tgt_metadata)
    sum(tgt_metadata$curated_celltype == combinedcelltype_pred_result[rownames(tgt_metadata), 'firstround_pred_celltype']) / nrow(tgt_metadata)
    table(tgt_metadata$curated_celltype, combinedcelltype_pred_result[rownames(tgt_metadata), 'firstround_pred_celltype'])
    table(tgt_metadata$curated_celltype, combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype'])
    cal_macroF1(tgt_metadata$curated_celltype, combinedcelltype_pred_result[rownames(tgt_metadata), 'firstround_pred_celltype'])
    cal_macroF1(tgt_metadata$curated_celltype, combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype'])
    ## run Cellcano predict on pancancer T
    #Cellcano predict -i ../Chu.PancanerCD4T_results/test_20_pancanT --trained_model output/pancanT_80_combinedtypes_trainedMLP_model -o output --prefix predict_Chu_pancanT_combinedtypes
    tgt_metadata = read.csv(file.path(dirname(result_dir), 'Chu.PancanerCD4T_results', 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
    combinedcelltype_pred_result = read.csv(file.path(result_dir, 'output', 'predict_Chu_pancanT_combinedtypescelltypes.csv'), header=T, row.names=1)
    sum(tgt_metadata$curated_celltype == combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype']) / nrow(tgt_metadata)
    sum(tgt_metadata$curated_celltype == combinedcelltype_pred_result[rownames(tgt_metadata), 'firstround_pred_celltype']) / nrow(tgt_metadata)
    table(tgt_metadata$curated_celltype, combinedcelltype_pred_result[rownames(tgt_metadata), 'firstround_pred_celltype'])
    table(tgt_metadata$curated_celltype, combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype'])

}

if (split_by == "80_enlarge:20") { ## make sure cell type proportion is the same
    set.seed(1234)
    train_idx = sample(1:ncol(pancancerT_obj), round(0.8*ncol(pancancerT_obj)))
    test_idx = setdiff(1:ncol(pancancerT_obj), train_idx)
    train_obj = pancancerT_obj[, train_idx]
    # sample per cell type
    sampled_idx = c()
    for (celltype in levels(train_obj$curated_celltype)) {
        celltype_idx = which(train_obj$curated_celltype == celltype)
        sampled_num = round(ncol(train_obj)/length(levels(train_obj$curated_celltype)))
        sampled_idx = c(sampled_idx, sample(celltype_idx, size=sampled_num, replace=T))
    }
    train_sampled_counts = train_obj@assays$RNA@counts[, sampled_idx]
    writeMM(train_sampled_counts, file.path(result_dir, 'train_80_enlarge_pancanT.mtx'))
    write(rownames(train_sampled_counts), file.path(result_dir, 'train_80_enlarge_pancanT_genes.tsv'))
    write(colnames(train_sampled_counts), file.path(result_dir, 'train_80_enlarge_pancanT_barcodes.tsv'))
    write.csv(train_obj@meta.data[sampled_idx, ], file.path(result_dir, 'train_80_enlarge_pancanT_metadata.csv'), quote=F)
    #---change the column names of metadata and run Cellcano
    #Cellcano train -i train_80_enlarge_pancanT -m train_80_enlarge_pancanT_metadata.csv --prefix pancanT_80_enlarge_combinedtypes_trained
    ## run Cellcano predict
    #Cellcano predict -i test_20_pancanT --trained_model output/pancanT_80_enlarge_combinedtypes_trainedMLP_model -o output --prefix predict_20_enlarge_pancanT_combinedtypes
    tgt_metadata = read.csv(file.path(result_dir, 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
    combinedcelltype_pred_result = read.csv(file.path(result_dir, 'output', 'predict_20_enlarge_pancanT_combinedtypescelltypes.csv'), header=T, row.names=1)
    sum(tgt_metadata$curated_celltype == combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype']) / nrow(tgt_metadata)
    sum(tgt_metadata$curated_celltype == combinedcelltype_pred_result[rownames(tgt_metadata), 'firstround_pred_celltype']) / nrow(tgt_metadata)
    table(tgt_metadata$curated_celltype, combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype'])
}


