# conda activate ~/envs/scGPT_bak/ # corrupted, need reinstall
# conda activate ~/envs/spatialATAC

import os
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

scPanKD_result_dir = '/net/mulan/home/wenjinma/projects/scPanKD/results/scPanKD/ProjecTILs_CD8_multibatch_cancertrimmed_to_HNSC_CD8'
scPanKD_adata = anndata.read_h5ad(os.path.join(scPanKD_result_dir, 'scPanKD_predicted_reference_features.h5ad'))
Seurat_result_dir = "/net/mulan/home/wenjinma/projects/scPanKD/results/Seurat/ProjecTILs_CD8_multibatch_cancertrimmed_to_HNSC_CD8"
Seurat_integrated = pd.read_csv(Seurat_result_dir+os.sep+'Seurat_predicted_ref_obj_embeddings.csv', index_col=0)
Seurat_umap = pd.read_csv(Seurat_result_dir+os.sep+'Seurat_predicted_ref_obj_umap_embeddings.csv', index_col=0)
# --- load normalized counts ---
mtx_prefix='/net/zootopia/disk1/wenjinma/data/scPanKD_data/ProjecTILs_CD8T/ProjecTIL_CD8T'
adata = anndata.read_mtx(mtx_prefix+'.mtx.gz').T
cells = pd.read_csv(mtx_prefix+'_barcodes.tsv', header=None, sep='\t')
adata.obs["barcode"] = cells[0].values
adata.obs_names = adata.obs['barcode']
adata.obs_names_make_unique(join="-")
adata.obs.index.name = None
adata = adata[scPanKD_adata.obs_names, :]
adata.obs = scPanKD_adata.obs.loc[adata.obs_names, :]
adata.obsm['scPanKD_embedding'] = scPanKD_adata.X
adata.obsm['scPanKD_umap'] = scPanKD_adata.obsm['X_umap']
adata.obsm['Seurat_embedding'] = Seurat_integrated.loc[adata.obs_names, :].values
adata.obsm['Seurat_umap'] = Seurat_umap.loc[adata.obs_names, :].values
sc.pp.normalize_per_cell(adata)

# --- load Seurat results ---
bm = Benchmarker(
    adata,
    batch_key="Tissue",
    label_key="celltype",
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    embedding_obsm_keys=["scPanKD_embedding", "scPanKD_umap", "Seurat_embedding", "Seurat_umap"],
    n_jobs=1,
)
bm.benchmark()
scPanKD_Tissue_res = bm.get_results(min_max_scale=False)
scPanKD_Tissue_res.to_csv(os.path.join(scPanKD_result_dir, 'scPanKD_Tissue_benchmark_results.csv'))

bm = Benchmarker(
    adata,
    batch_key="Cohort",
    label_key="celltype",
    embedding_obsm_keys=["scPanKD_embedding", "scPanKD_umap", "Seurat_embedding", "Seurat_umap"],
    n_jobs=1,
)
bm.benchmark()
scPanKD_Cohort_res = bm.get_results(min_max_scale=False)
scPanKD_Cohort_res.to_csv(os.path.join(scPanKD_result_dir, 'scPanKD_Cohort_benchmark_results.csv'))

