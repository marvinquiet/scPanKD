.libPaths("/net/mulan/home/wenjinma/Rlib")
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)

project_dir = "/net/mulan/home/wenjinma/minicluster/projects/test_Cellcano_XuefengWang"
data_dir = "/net/mulan/home/wenjinma/minicluster/data/pan-cancer-T"
result_dir = file.path(project_dir, 'ZhengCD4T_to_Chu_results_curated')
dir.create(result_dir, showWarnings=F)

# read Zheng data
data_obj = readRDS(file.path(project_dir, 'Zheng.PancanerT', 'CD4.rds'))
selected_idx = data_obj$Annotation %in% c('CD4.EOMES+', 'CD4.Naivelike', 'Th17', 'Tfh', 'Treg')
data_obj = data_obj[,selected_idx] ## selected cell types, 10911 cells remained
data_obj$curated_celltype = data_obj$Annotation

# filter genes
integrated_counts = GetAssayData(data_obj, slot="data", assay="integrated")   
selected_genes = rownames(integrated_counts)
counts = GetAssayData(data_obj, slot="counts", assay="RNA")
counts = counts[selected_genes,] ## 800 genes x 10911 cells

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

writeMM(pancancerT_obj@assays$RNA@counts, file.path(result_dir, 'Zheng_train_pancanT.mtx'))
write(rownames(pancancerT_obj@assays$RNA@counts), file.path(result_dir, 'Zheng_train_pancanT_genes.tsv'))
write(colnames(pancancerT_obj@assays$RNA@counts), file.path(result_dir, 'Zheng_train_pancanT_barcodes.tsv'))
write.csv(pancancerT_obj@meta.data, file.path(result_dir, 'Zheng_train_pancanT_metadata.csv'), quote=F)
## run Cellcano
#Cellcano train -i Zheng_train_pancanT -m Zheng_train_pancanT_metadata.csv --prefix Zheng_pancanT_trained
 
# read Chu data
load(file.path(project_dir, 'Chu.PancanerT', 'CD4.keep.RData'))
pancancerT_obj = readRDS(file.path(data_dir, 'CD4'))
CD4_keep_cells = rownames(CD4.keep)
common_cells = intersect(colnames(pancancerT_obj), CD4_keep_cells) # 94571 cells
pancancerT_obj = pancancerT_obj[, common_cells]
pancancerT_obj$curated_celltype = CD4.keep[common_cells, 'Annotation']
idx = pancancerT_obj$curated_celltype == 'D4.Naivelike'
pancancerT_obj$curated_celltype[idx] = 'CD4.Naivelike'
selected_idx = pancancerT_obj$curated_celltype %in% c('CD4.CTL', 'CD4.Naivelike', 'Tfh', 'Th17', 'Treg') # 54079 cells
pancancerT_obj = pancancerT_obj[, selected_idx]

# filter genes
counts = GetAssayData(pancancerT_obj, slot="counts", assay="RNA")   
genes_percent = rowMeans(counts > 0)
genes_selected = genes_percent > 0.01 
pancancerT_obj = pancancerT_obj[genes_selected] # remaining 12647 genes

# perform data analysis
pancancerT_obj = NormalizeData(pancancerT_obj)
pancancerT_obj = FindVariableFeatures(pancancerT_obj, selection.method = "vst", nfeatures = 2000)
pancancerT_obj = ScaleData(pancancerT_obj, features = rownames(pancancerT_obj))
pancancerT_obj = RunPCA(pancancerT_obj, features = VariableFeatures(object = pancancerT_obj))
pancancerT_obj = FindNeighbors(pancancerT_obj, dims=1:20)
pancancerT_obj = FindClusters(pancancerT_obj, resolution=0.5)
pancancerT_obj = RunUMAP(pancancerT_obj, dims=1:20)
gg1 = DimPlot(pancancerT_obj, reduction='pca', group.by='curated_celltype', label=T)
gg2 = DimPlot(pancancerT_obj, reduction='pca', group.by='CancerType', label=T)
gg1 + gg2
ggsave(file.path(result_dir, 'Chu_PCA_analysis.png'), width=10, height=5)
gg1 = DimPlot(pancancerT_obj, reduction='umap', group.by='curated_celltype', label=T)
gg2 = DimPlot(pancancerT_obj, reduction='umap', group.by='CancerType', label=T)
gg1 + gg2
ggsave(file.path(result_dir, 'Chu_UMAP_analysis.png'), width=10, height=5)

writeMM(pancancerT_obj@assays$RNA@counts, file.path(result_dir, 'Chu_pancanT.mtx'))
write(rownames(pancancerT_obj@assays$RNA@counts), file.path(result_dir, 'Chu_pancanT_genes.tsv'))
write(colnames(pancancerT_obj@assays$RNA@counts), file.path(result_dir, 'Chu_pancanT_barcodes.tsv'))
write.csv(pancancerT_obj@meta.data, file.path(result_dir, 'Chu_pancanT_metadata.csv'), quote=F)
#Cellcano predict -i Chu_pancanT --trained_model output/Zheng_pancanT_trainedMLP_model -o output --prefix predict_Chu_pancanT
 

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

tgt_metadata = read.csv(file.path(result_dir, 'Chu_pancanT_metadata.csv'), header=T, row.names=1)
idx = tgt_metadata$curated_celltype == 'CD4.CTL'
tgt_metadata$curated_celltype[idx] = 'CD4.EOMES+'
pred_result = read.csv(file.path(result_dir, 'output', 'predict_Chu_pancanTcelltypes.csv'), header=T, row.names=1)
sum(tgt_metadata$curated_celltype == pred_result[rownames(tgt_metadata), 'pred_celltype']) / nrow(tgt_metadata)
sum(tgt_metadata$curated_celltype == pred_result[rownames(tgt_metadata), 'firstround_pred_celltype']) / nrow(tgt_metadata)
table(tgt_metadata$curated_celltype, pred_result[rownames(tgt_metadata), 'pred_celltype'])

cal_macroF1(tgt_metadata$curated_celltype, pred_result[rownames(tgt_metadata), 'firstround_pred_celltype'])
cal_macroF1(tgt_metadata$curated_celltype, pred_result[rownames(tgt_metadata), 'pred_celltype'])

