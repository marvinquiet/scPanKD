# conda activate ~/envs/spatialATAC

.libPaths("/net/mulan/home/wenjinma/Rlib")
data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data"
require(Seurat)
require(Matrix)
require(ggplot2)

# --- generate Chu CD4T integrated data
Chu_CD4_data_dir = file.path(data_dir, 'Chu_pancancer_CD4T')
Chu_CD4_train_80_counts = readMM(file.path(Chu_CD4_data_dir, 'train_80_pancanT.mtx.gz'))
Chu_CD4_train_genes = scan(file.path(Chu_CD4_data_dir, 'train_80_pancanT_genes.tsv'), what=character())
Chu_CD4_train_barcodes = scan(file.path(Chu_CD4_data_dir, 'train_80_pancanT_barcodes.tsv'), what=character())
rownames(Chu_CD4_train_80_counts) = Chu_CD4_train_genes
colnames(Chu_CD4_train_80_counts) = Chu_CD4_train_barcodes
Chu_CD4_train_metadata = read.csv(file.path(Chu_CD4_data_dir, 'CancerType', 'train_80_pancanT_batch_metadata.csv'), header=T, row.names=1)
Chu_CD4_train_object = CreateSeuratObject(counts=Chu_CD4_train_80_counts[, rownames(Chu_CD4_train_metadata)], meta.data=Chu_CD4_train_metadata)
Chu_CD4_train_object_list = SplitObject(Chu_CD4_train_object, split.by = "CancerType")
Chu_CD4_train_object_list = lapply(X = Chu_CD4_train_object_list, FUN = function(x) {
                    x = NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
                    x = FindVariableFeatures(x)
                    x = ScaleData(x)
                    })
features = SelectIntegrationFeatures(object.list=Chu_CD4_train_object_list)
anchors = FindIntegrationAnchors(object.list = Chu_CD4_train_object_list, anchor.features = features)
Chu_CD4_train_object_integrated = IntegrateData(anchorset = anchors)
# integrate trainining object using tissue as batch
Chu_CD4_train_object_integrated = ScaleData(Chu_CD4_train_object_integrated)
Chu_CD4_train_object_integrated = RunPCA(Chu_CD4_train_object_integrated)
Chu_CD4_train_object = Chu_CD4_train_object_integrated
saveRDS(Chu_CD4_train_object, file=file.path(Chu_CD4_data_dir, 'Chu_CD4T_train_80_integrated.RDS'))
# --- plot the integrated object
Chu_CD4_train_object = RunUMAP(Chu_CD4_train_object, reduction = "pca", dims = 1:30)
p1 = DimPlot(Chu_CD4_train_object, reduction = "umap", group.by = "CancerType")
ggsave(file.path(Chu_CD4_data_dir, 'Seurat_integrated_Chu_CD4T_CancerType.png'))
p2 = DimPlot(Chu_CD4_train_object, reduction = "umap", group.by = "celltype", repel = TRUE)
ggsave(file.path(Chu_CD4_data_dir, 'Seurat_integrated_Chu_CD4T_celltype.png'))

# --- generate Chu CD8T integrated data
Chu_CD8_data_dir = file.path(data_dir, 'Chu_pancancer_CD8T')
Chu_CD8_train_80_counts = readMM(file.path(Chu_CD8_data_dir, 'train_80_pancanT.mtx.gz'))
Chu_CD8_train_genes = scan(file.path(Chu_CD8_data_dir, 'train_80_pancanT_genes.tsv'), what=character())
Chu_CD8_train_barcodes = scan(file.path(Chu_CD8_data_dir, 'train_80_pancanT_barcodes.tsv'), what=character())
rownames(Chu_CD8_train_80_counts) = Chu_CD8_train_genes
colnames(Chu_CD8_train_80_counts) = Chu_CD8_train_barcodes
Chu_CD8_train_metadata = read.csv(file.path(Chu_CD8_data_dir, 'CancerType', 'train_batch_metadata.csv'), header=T, row.names=1)
Chu_CD8_train_object = CreateSeuratObject(counts=Chu_CD8_train_80_counts[, rownames(Chu_CD8_train_metadata)], meta.data=Chu_CD8_train_metadata)
Chu_CD8_train_object_list = SplitObject(Chu_CD8_train_object, split.by = "CancerType")
Chu_CD8_train_object_list = lapply(X = Chu_CD8_train_object_list, FUN = function(x) {
                    x = NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
                    x = FindVariableFeatures(x)
                    x = ScaleData(x)
                    })
features = SelectIntegrationFeatures(object.list=Chu_CD8_train_object_list)
anchors = FindIntegrationAnchors(object.list = Chu_CD8_train_object_list, anchor.features = features)
Chu_CD8_train_object_integrated = IntegrateData(anchorset = anchors)
# integrate trainining object using tissue as batch
Chu_CD8_train_object_integrated = ScaleData(Chu_CD8_train_object_integrated)
Chu_CD8_train_object_integrated = RunPCA(Chu_CD8_train_object_integrated)
Chu_CD8_train_object = Chu_CD8_train_object_integrated
saveRDS(Chu_CD8_train_object, file=file.path(Chu_CD8_data_dir, 'Chu_CD8T_train_80_integrated.RDS'))
# --- plot the integrated object
Chu_CD8_train_object = RunUMAP(Chu_CD8_train_object, reduction = "pca", dims = 1:30)
p1 = DimPlot(Chu_CD8_train_object, reduction = "umap", group.by = "CancerType")
ggsave(file.path(Chu_CD8_data_dir, 'Seurat_integrated_Chu_CD8T_CancerType.png'))
p2 = DimPlot(Chu_CD8_train_object, reduction = "umap", group.by = "celltype", repel = TRUE)
ggsave(file.path(Chu_CD8_data_dir, 'Seurat_integrated_Chu_CD8T_celltype.png'))


# --- generate ProjecTILs CD8T integrated data - no trimmed, based on cancer type
ProjecTILs_CD8_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
ProjecTILs_counts = readMM(file.path(ProjecTILs_CD8_data_dir, 'ProjecTIL_CD8T.mtx.gz'))
ProjecTILs_genes = scan(file.path(ProjecTILs_CD8_data_dir, 'ProjecTIL_CD8T_genes.tsv'), what=character())
ProjecTILs_barcodes = scan(file.path(ProjecTILs_CD8_data_dir, 'ProjecTIL_CD8T_barcodes.tsv'), what=character())
rownames(ProjecTILs_counts) = ProjecTILs_genes
colnames(ProjecTILs_counts) = ProjecTILs_barcodes
ProjecTILs_metadata = read.csv(file.path(ProjecTILs_CD8_data_dir, 'ProjecTIL_CD8T_metadata.csv'), header=T, row.names=1)
ProjecTILs_object = CreateSeuratObject(counts=ProjecTILs_counts, meta.data=ProjecTILs_metadata[colnames(ProjecTILs_counts), ])
ProjecTILs_object_list = SplitObject(ProjecTILs_object, split.by = "Tissue")
ProjecTILs_object_list = lapply(X = ProjecTILs_object_list, FUN = function(x) {
                    x = NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
                    x = FindVariableFeatures(x)
                    x = ScaleData(x)
                    })
features = SelectIntegrationFeatures(object.list=ProjecTILs_object_list)
anchors = FindIntegrationAnchors(object.list = ProjecTILs_object_list, anchor.features = features)
ProjecTILs_object_integrated = IntegrateData(anchorset = anchors)
# integrate trainining object using tissue as batch
ProjecTILs_object_integrated = ScaleData(ProjecTILs_object_integrated)
ProjecTILs_object_integrated = RunPCA(ProjecTILs_object_integrated)
ProjecTILs_object = ProjecTILs_object_integrated
saveRDS(ProjecTILs_object, file=file.path(ProjecTILs_CD8_data_dir, 'ProjecTILs_CD8T_Tissue_integrated.RDS'))
# --- plot the integrated object
ProjecTILs_object = RunUMAP(ProjecTILs_object, reduction = "pca", dims = 1:30)
p1 = DimPlot(ProjecTILs_object, reduction = "umap", group.by = "Tissue")
ggsave(file.path(ProjecTILs_CD8_data_dir, 'Seurat_integrated_ProjecTILs_CD8T_Tissue.png'))
p2 = DimPlot(ProjecTILs_object, reduction = "umap", group.by = "celltype", repel = TRUE)
ggsave(file.path(ProjecTILs_CD8_data_dir, 'Seurat_integrated_ProjecTILs_CD8T_celltype.png'))


# --- generate ProjecTILs CD8T integrated data - trimmed, based on cancer type
ProjecTILs_CD8_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
studies = c('EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464')
ProjecTILs_object_list = list()
for (study in studies) {
    study_counts = readMM(file.path(ProjecTILs_CD8_data_dir, 'trimmed', paste0('ProjecTIL_', study, '_CD8T_trimmed.mtx.gz')))
    study_genes = scan(file.path(ProjecTILs_CD8_data_dir, 'trimmed', paste0('ProjecTIL_', study, '_CD8T_trimmed_genes.tsv')), what=character())
    study_barcodes = scan(file.path(ProjecTILs_CD8_data_dir, 'trimmed', paste0('ProjecTIL_', study, '_CD8T_trimmed_barcodes.tsv')), what=character())
    rownames(study_counts) = study_genes
    colnames(study_counts) = study_barcodes
    study_metadata = read.csv(file.path(ProjecTILs_CD8_data_dir, 'trimmed', paste0('ProjecTIL_', study, '_CD8T_trimmed_metadata.csv')), header=T, row.names=1)
    study_object = CreateSeuratObject(counts=study_counts, meta.data=study_metadata[colnames(study_counts), ])
    study_object = NormalizeData(study_object, normalization.method = "LogNormalize", scale.factor = 10000)
    study_object = FindVariableFeatures(study_object)
    study_object = ScaleData(study_object)
    ProjecTILs_object_list[[study]] = study_object
}
features = SelectIntegrationFeatures(object.list=ProjecTILs_object_list)
anchors = FindIntegrationAnchors(object.list = ProjecTILs_object_list, anchor.features = features)
ProjecTILs_object_integrated = IntegrateData(anchorset = anchors)
# integrate trainining object using tissue as batch
ProjecTILs_object_integrated = ScaleData(ProjecTILs_object_integrated)
ProjecTILs_object_integrated = RunPCA(ProjecTILs_object_integrated)
ProjecTILs_object = ProjecTILs_object_integrated
saveRDS(ProjecTILs_object, file=file.path(ProjecTILs_CD8_data_dir, 'ProjecTILs_CD8T_trimmed_Cohort_integrated.RDS'))
# --- plot the integrated object
ProjecTILs_object = RunUMAP(ProjecTILs_object, reduction = "pca", dims = 1:30)
p1 = DimPlot(ProjecTILs_object, reduction = "umap", group.by = "Cohort")
ggsave(file.path(ProjecTILs_CD8_data_dir, 'Seurat_integrated_ProjecTILs_CD8T_trimmed_Cohort.png'))
p2 = DimPlot(ProjecTILs_object, reduction = "umap", group.by = "celltype", repel = TRUE)
ggsave(file.path(ProjecTILs_CD8_data_dir, 'Seurat_integrated_ProjecTILs_CD8T_trimmed_celltype.png'))


# --- generate Zheng CD4T integrated data
Zheng_CD4_data_dir = file.path(data_dir, 'Zheng_CD4T')
Zheng_CD4_counts = readMM(file.path(Zheng_CD4_data_dir, 'Zheng_CD4T.mtx.gz'))
Zheng_CD4_genes = scan(file.path(Zheng_CD4_data_dir, 'Zheng_CD4T_genes.tsv'), what=character())
Zheng_CD4_barcodes = scan(file.path(Zheng_CD4_data_dir, 'Zheng_CD4T_barcodes.tsv'), what=character())
rownames(Zheng_CD4_counts) = Zheng_CD4_genes
colnames(Zheng_CD4_counts) = Zheng_CD4_barcodes
Zheng_CD4_metadata = read.csv(file.path(Zheng_CD4_data_dir, 'Zheng_CD4T_metadata.csv'), header=T, row.names=1)
Zheng_CD4_object = CreateSeuratObject(counts=Zheng_CD4_counts, meta.data=Zheng_CD4_metadata[colnames(Zheng_CD4_counts), ])
Zheng_CD4_object_list = SplitObject(Zheng_CD4_object, split.by = "cancerType")
Zheng_CD4_object_list = lapply(X = Zheng_CD4_object_list, FUN = function(x) {
                    x = NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
                    x = FindVariableFeatures(x)
                    x = ScaleData(x)
                    })
features = SelectIntegrationFeatures(object.list=Zheng_CD4_object_list)
anchors = FindIntegrationAnchors(object.list = Zheng_CD4_object_list, anchor.features = features)
Zheng_CD4_object_integrated = IntegrateData(anchorset = anchors)
# integrate trainining object using tissue as batch
Zheng_CD4_object_integrated = ScaleData(Zheng_CD4_object_integrated)
Zheng_CD4_object_integrated = RunPCA(Zheng_CD4_object_integrated)
Zheng_CD4_object = Zheng_CD4_object_integrated
saveRDS(Zheng_CD4_object, file=file.path(Zheng_CD4_data_dir, 'Zheng_CD4T_cancerType_integrated.RDS'))
# --- plot the integrated object
Zheng_CD4_object = RunUMAP(Zheng_CD4_object, reduction = "pca", dims = 1:30)
p1 = DimPlot(Zheng_CD4_object, reduction = "umap", group.by = "cancerType")
ggsave(file.path(Zheng_CD4_data_dir, 'Seurat_integrated_Zheng_CD4T_cancerType.png'))
p2 = DimPlot(Zheng_CD4_object, reduction = "umap", group.by = "celltype", repel = TRUE)
ggsave(file.path(Zheng_CD4_data_dir, 'Seurat_integrated_Zheng_CD4T_celltype.png'))