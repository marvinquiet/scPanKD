.libPaths("/net/mulan/home/wenjinma/Rlib")
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)

split_by = '80_enlarge:20' # 80:20 | tissue | 80_enlarge:20

project_dir = "/net/mulan/home/wenjinma/minicluster/projects/test_Cellcano_XuefengWang"
data_dir = "/net/mulan/home/wenjinma/minicluster/data/pan-cancer-T"
result_dir = file.path(project_dir, 'pancancerT_ref_data')
dir.create(result_dir, showWarnings=F)

pancancerT_obj = readRDS(file.path(data_dir, 'CD8'))
pancancerT_obj$curated_celltype = as.character(pancancerT_obj$cell.type)
idx = pancancerT_obj$cell.type %in% c('CD8_c3_Tn', 'CD8_c13_Tn_TCF7')
pancancerT_obj@meta.data[idx, 'curated_celltype'] = 'naive'
idx = pancancerT_obj$cell.type %in% c('CD8_c0_t-Teff', 'CD8_c2_Teff', 'CD8_c8_Teff_KLRG1', 'CD8_c10_Teff_CD244', 'CD8_c11_Teff_SEMA4A')
pancancerT_obj@meta.data[idx, 'curated_celltype'] = 'effector'
pancancerT_obj$curated_celltype = as.factor(pancancerT_obj$curated_celltype)
selected_tissues = names(which(table(pancancerT_obj$CancerType) > 1000))
filtered_obj = subset(pancancerT_obj, subset=CancerType %in% selected_tissues)

# perform data analysis
filtered_obj = NormalizeData(filtered_obj)
filtered_obj = FindVariableFeatures(filtered_obj, selection.method = "vst", nfeatures = 2000)
filtered_obj = ScaleData(filtered_obj, features = rownames(filtered_obj))
filtered_obj = RunPCA(filtered_obj, features = VariableFeatures(object = filtered_obj))
filtered_obj = FindNeighbors(filtered_obj, dims=1:20)
filtered_obj = FindClusters(filtered_obj, resolution=0.5)
filtered_obj = RunUMAP(filtered_obj, dims=1:20)
gg1 = DimPlot(filtered_obj, reduction='pca', group.by='curated_celltype', label=T)
gg2 = DimPlot(filtered_obj, reduction='pca', group.by='CancerType', label=T)
gg1 + gg2
ggsave(file.path(result_dir, 'PCA_analysis.png'), width=10, height=5)
gg1 = DimPlot(filtered_obj, reduction='umap', group.by='curated_celltype', label=T)
gg1 + gg2
ggsave(file.path(result_dir, 'UMAP_analysis.png'), width=10, height=5)

# run data split
if (split_by == 'tissue') { # split by cancer tissue
    # random sample
    set.seed(1234)
    train_tissues = sample(selected_tissues, 5)
    test_tissues = sample(setdiff(selected_tissues, train_tissues), 5)
    train_obj = subset(filtered_obj, subset=CancerType %in% train_tissues)
    test_obj = subset(filtered_obj, subset=CancerType %in% test_tissues)
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
    train_idx = sample(1:ncol(filtered_obj), round(0.8*ncol(filtered_obj)))
    test_idx = setdiff(1:ncol(filtered_obj), train_idx)
    train_obj = filtered_obj[, train_idx]
    test_obj = filtered_obj[, test_idx]
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
    table(tgt_metadata$curated_celltype, combinedcelltype_pred_result[rownames(tgt_metadata), 'pred_celltype'])
}

if (split_by == "80_enlarge:20") { ## make sure cell type proportion is the same
    set.seed(1234)
    train_idx = sample(1:ncol(filtered_obj), round(0.8*ncol(filtered_obj)))
    test_idx = setdiff(1:ncol(filtered_obj), train_idx)
    train_obj = filtered_obj[, train_idx]
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


