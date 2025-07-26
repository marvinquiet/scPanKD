## conda activate ~/envs/Cellcano
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)

project_dir = "/net/mulan/home/wenjinma/minicluster/projects/test_Cellcano_XuefengWang"
result_dir = file.path(project_dir, 'CD8T_ref_data')
dir.create(result_dir, showWarnings=F)
load(file.path(project_dir, 'CD8T.ref.rawcount.RData'))
# analyze data
rownames(CD8T.ref.label) = CD8T.ref.label$CellID
colnames(CD8T.ref.label) = c('CellID', 'celltype')
CD8T.ref.label$batch = unlist(sapply(strsplit(CD8T.ref.label$CellID, split='_'), '[', 1))

#' Seurat analysis pipeline to check batch information
obj = CreateSeuratObject(counts=CD8T.ref.rawcount, meta.data=CD8T.ref.label)
obj = NormalizeData(obj)
obj = FindVariableFeatures(obj)
obj = ScaleData(obj)
obj = RunPCA(obj)
obj = FindNeighbors(obj, dims=1:20)
obj = FindClusters(obj, resolution=0.5)
obj = RunUMAP(obj, dims=1:20)
gg1 = DimPlot(obj, reduction='pca', group.by='batch', label=T)
gg2 = DimPlot(obj, reduction='pca', group.by='celltype', label=T)
gg1 + gg2
ggsave(file.path(result_dir, 'PCA_analysis.png'), width=10, height=5)
gg1 = DimPlot(obj, reduction='umap', group.by='batch', label=T)
gg2 = DimPlot(obj, reduction='umap', group.by='celltype', label=T)
gg1 + gg2
ggsave(file.path(result_dir, 'UMAP_analysis.png'), width=10, height=5)

#' split the data into different batches
for (batch in unique(CD8T.ref.label$batch)) {
    batch_df = CD8T.ref.label[CD8T.ref.label$batch == batch, ]
    batch_data = CD8T.ref.rawcount[, rownames(batch_df)]
    batch_mat = Matrix(batch_data, sparse=T)
    writeMM(batch_mat, file.path(result_dir, paste0(batch, '.mtx')))
    write(colnames(batch_mat), file.path(result_dir, paste0(batch, '_barcodes.tsv')))
    write(rownames(batch_mat), file.path(result_dir, paste0(batch, '_genes.tsv')))
}


#load(file.path(project_dir, 'GSE179994.CD8T.ref.rawcount.RData'))
#ref_data = Matrix(CD8T.ref.mat, sparse=T)
#ref_data = Matrix(CD8T.ref.rawcount, sparse=T)
ref_data = Matrix(GSE179994.CD8T.ref.rawcount.mat, sparse=T)
#ref_data = round(ref_data, 3)
writeMM(ref_data, file.path(project_dir, 'GSE179994_CD8T_ref.mtx'))
write(rownames(ref_data), file.path(project_dir, 'GSE179994_CD8T_ref_genes.tsv'))
write(colnames(ref_data), file.path(project_dir, 'GSE179994_CD8T_ref_barcodes.tsv'))

#rownames(CD8T.ref.label) = CD8T.ref.label$CellID
#colnames(CD8T.ref.label) = c('CellID', 'celltype')
#write.csv(CD8T.ref.label, file.path(project_dir, 'CD8T_ref_metadata.csv'), quote=F)

rownames(GSE179994.CD8T.ref.label) = GSE179994.CD8T.ref.label$CellID
colnames(GSE179994.CD8T.ref.label) = c('CellID', 'celltype')
write.csv(GSE179994.CD8T.ref.label, file.path(project_dir, 'GSE179994_CD8T_ref_metadata.csv'), quote=F)


## run Cellcano
#Cellcano train -i CD8T_ref -m CD8T_ref_metadata.csv --prefix CD8T_trained
#Cellcano train -i GSE179994_CD8T_ref -m GSE179994_CD8T_ref_metadata.csv --prefix GSE179994_CD8T_trained

load(file.path(project_dir, 'HNSC.CD8T.rawcount.RData'))
#tgt_data = Matrix(HNSC.CD8T.mat, sparse=T)
tgt_data = Matrix(HNSC.CD8T.rawcount, sparse=T)
#tgt_data = round(tgt_data, 3)
writeMM(tgt_data, file.path(project_dir, 'HNSC_CD8T.mtx'))
write(rownames(tgt_data), file.path(project_dir, 'HNSC_CD8T_genes.tsv'))
write(colnames(tgt_data), file.path(project_dir, 'HNSC_CD8T_barcodes.tsv'))

rownames(HNSC.CD8T.label) = HNSC.CD8T.label$CellID
write.csv(HNSC.CD8T.label, file.path(project_dir, 'HNSC_CD8T_metadata.csv'), quote=F)

## run Cellcano predict
#Cellcano predict -i HNSC_CD8T --trained_model output/CD8T_trainedMLP_model -o output --prefix predict_HNSC_CD8
#Cellcano predict -i HNSC_CD8T --trained_model output/GSE179994_CD8T_trainedMLP_model -o output --prefix predict_HNSC_CD8

## analyze results
tgt_metadata = read.csv(file.path(project_dir, 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
pred_result = read.csv(file.path(project_dir, 'output', 'predict_HNSC_CD8celltypes.csv'), header=T, row.names=1)
sum(tgt_metadata$Label == pred_result[rownames(tgt_metadata), 'pred_celltype']) / nrow(tgt_metadata)

