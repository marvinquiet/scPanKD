# write out embeddings obtained from Seurat
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(Seurat)

Seurat_result_dir = "/net/mulan/home/wenjinma/projects/scPanKD/results/Seurat/ProjecTILs_CD8_multibatch_cancertrimmed_to_HNSC_CD8"
Seurat_obj = readRDS(file.path(Seurat_result_dir, 'Seurat_predicted_ref_obj.rds'))

# --- write out embeddings
write.csv(Seurat_obj@meta.data, file.path(Seurat_result_dir, 'Seurat_predicted_ref_obj_metadata.csv'), quote=F)
embeddings = Seurat_obj@assays$integrated@data
write.csv(t(embeddings), file.path(Seurat_result_dir, 'Seurat_predicted_ref_obj_embeddings.csv'), quote=F)
umap_embeddings = Embeddings(Seurat_obj, 'umap')
write.csv(umap_embeddings, file.path(Seurat_result_dir, 'Seurat_predicted_ref_obj_umap_embeddings.csv'), quote=F)