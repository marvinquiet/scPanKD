# conda activate /net/zootopia/disk1/wenjinma/envs/SCimilarity/

import os, sys 
import argparse                                                                                                                                                                                                                                                                                                    
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=ResourceWarning)
import anndata
import scipy
import scanpy as sc
import numpy as np
import pandas as pd

import pipelines.load_utils as load_utils
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset, lognorm_counts

import matplotlib.pyplot as plt
plt.style.context('default')
from sklearn.cluster import KMeans
# calculate ARI
from sklearn.metrics import adjusted_rand_score

data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data/"
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD/"
result_dir = project_dir+os.sep+'results'+os.sep+'SCimilarity'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# add argparse
parser = argparse.ArgumentParser(description='Run scGPT pipeline for single-slice data.')
parser.add_argument('--experiment', type=str,
                    help='Load experiment related dataset.')
args = parser.parse_args()
experiment = args.experiment

SCimilarity_result_dir = result_dir+os.sep+experiment
if not os.path.exists(SCimilarity_result_dir):
    os.makedirs(SCimilarity_result_dir)

# --- SCimilarity model directory ---
SCimilarity_dir = "/net/zootopia/disk1/wenjinma/LLM_models/SCimilarity"
SCimilarity_model_dir = SCimilarity_dir+os.sep+'annotation_model_v1'
ca = CellAnnotation(model_path=SCimilarity_model_dir)

# Aim: use train embedding to train a SVM/XGBoost classifier and then predict on test embedding
if experiment == 'Chu_CD4T_validation':
    # --- curated_celltype
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'
    train_adata = load_utils.load_mtx(Chu_CD4T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata = load_utils.load_mtx(Chu_CD4T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x')) 
    celltype_col = 'curated_celltype' 

if experiment == 'Chu_CD8T_validation':
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    train_adata = load_utils.load_mtx(Chu_CD8T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata = load_utils.load_mtx(Chu_CD8T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))  
    test_adata.obs['celltype'] = test_adata.obs['curated_celltype'].values
    celltype_col = 'celltype' 
    
if experiment == 'Zheng_CD4_to_Chu_CD4':
    Zheng_CD4_dir = data_dir+os.sep+'Zheng_CD4T'
    train_adata = load_utils.load_mtx(Zheng_CD4_dir+os.sep+'Zheng_CD4T')
    train_metadata = pd.read_csv(Zheng_CD4_dir+os.sep+'Zheng_CD4T_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'
    test_adata = load_utils.load_mtx(Chu_CD4T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'
    
if experiment == 'GSE179994_CD8_to_HNSC_CD8':
    GSE179994_CD8_dir = data_dir+os.sep+'GSE179994_CD8T'
    train_adata = load_utils.load_mtx(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref')
    train_metadata = pd.read_csv(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    HNSC_CD8_dir = data_dir+os.sep+'HNSC_CD8T'
    test_adata = load_utils.load_mtx(HNSC_CD8_dir+os.sep+'HNSC_CD8T')
    test_metadata = pd.read_csv(HNSC_CD8_dir+os.sep+'HNSC_CD8T_metadata_curated.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'

if experiment == 'ProjecTILs_CD8_to_HNSC_CD8':
    ProjecTILs_CD8_dir = data_dir+os.sep+'ProjecTILs_CD8T'
    train_adata = load_utils.load_mtx(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T')
    train_metadata = pd.read_csv(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    HNSC_CD8_dir = data_dir+os.sep+'HNSC_CD8T'
    test_adata = load_utils.load_mtx(HNSC_CD8_dir+os.sep+'HNSC_CD8T')
    test_metadata = pd.read_csv(HNSC_CD8_dir+os.sep+'HNSC_CD8T_metadata_curated.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'


if experiment == 'GSE179994_CD8_to_Chu_CD8T':
    GSE179994_CD8_dir = data_dir+os.sep+'GSE179994_CD8T'
    train_adata = load_utils.load_mtx(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref')
    train_metadata = pd.read_csv(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    test_adata = load_utils.load_mtx(Chu_CD8T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata.obs['celltype'] = test_adata.obs['curated_celltype'].values
    celltype_col = 'celltype' 

if experiment == 'ProjecTILs_CD8_to_Chu_CD8T':
    ProjecTILs_CD8_dir = data_dir+os.sep+'ProjecTILs_CD8T'
    train_adata = load_utils.load_mtx(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T')
    train_metadata = pd.read_csv(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    test_adata = load_utils.load_mtx(Chu_CD8T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata.obs['celltype'] = test_adata.obs['curated_celltype'].values
    celltype_col = 'celltype'


def impute_genes(adata):
    reference_genes = list(ca.gene_order)
    missing_genes = set(reference_genes) - set(adata.var_names)
    n_cells = adata.n_obs
    zero_data = np.zeros((n_cells, len(missing_genes))) # impute zero to avoid limited overlapping panel
    X_missing = scipy.sparse.csr_matrix(zero_data)
    missing_adata = anndata.AnnData(
        X=X_missing,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=list(missing_genes))
    )
    adata_impute = anndata.concat([adata, missing_adata], axis=1)
    adata_impute.obs = adata.obs
    adata_impute.var['gene_name'] = adata_impute.var_names.astype(str).values
    return adata_impute


common_genes = train_adata.var_names.intersection(ca.gene_order)
if len(common_genes) < 5000:
    train_adata = impute_genes(train_adata).copy()
train_adata.layers['counts'] = train_adata.X.copy()  # save the original counts
train_adata = align_dataset(train_adata, ca.gene_order)
train_adata = lognorm_counts(train_adata)
train_adata.obsm["X_scimilarity"] = ca.get_embeddings(train_adata.X)

common_genes = test_adata.var_names.intersection(ca.gene_order)
if len(common_genes) < 5000:
    test_adata= impute_genes(test_adata).copy()
test_adata.layers['counts'] = test_adata.X.copy()  # save the original counts
test_adata = align_dataset(test_adata, ca.gene_order)
test_adata = lognorm_counts(test_adata)
test_adata.obsm["X_scimilarity"] = ca.get_embeddings(test_adata.X)

# train a XGBoost classifier on train_adata.obsm['X_scimilarity'] and train_adata.obs[celltype_col]
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
# Convert celltype to categorical and get the mapping
celltype_categorical = train_adata.obs[celltype_col].astype('category')
celltype_codes = celltype_categorical.cat.codes
celltype_categories = celltype_categorical.cat.categories

# Split the data into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(
    train_adata.obsm['X_scimilarity'],
    celltype_codes,  # Use the categorical codes
    test_size=0.2,
    random_state=42,
    stratify=train_adata.obs[celltype_col]
)
# Initialize the XGBoost classifier
xgb_classifier = XGBClassifier(
    n_estimators=100,
    max_depth=6,
    learning_rate=0.1,
    random_state=42,
    use_label_encoder=False,
    eval_metric='mlogloss'
)
# Train the classifier
xgb_classifier.fit(X_train, y_train)
# Predict on the validation set
y_pred = xgb_classifier.predict(X_val)
# Evaluate the classifier
print("Classification Report:")
print(classification_report(y_val, y_pred))
print("Confusion Matrix:")
print(confusion_matrix(y_val, y_pred))
# Predict on the test_adata.obsm['X_scimilarity']
test_predictions_codes = xgb_classifier.predict(test_adata.obsm['X_scimilarity'])
# Convert predictions back to categorical labels using the category mapping
test_predictions_labels = celltype_categories[test_predictions_codes]
# Save both the numeric codes and categorical labels to the test_adata.obs
test_adata.obs['predicted_celltype_code'] = test_predictions_codes  # [0,1,2,3] format
test_adata.obs['pred_celltype'] = test_predictions_labels       # original categorical labels
# generate confusion matrix
confusion_matrix_result = confusion_matrix(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'])
# --- calculate ARI
ari_score = adjusted_rand_score(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'])
print(f"Adjusted Rand Index (ARI) score: {ari_score}")
test_adata.obs.to_csv(SCimilarity_result_dir+os.sep+'SCimilarity_predicted_celltypes.csv')


# # --- load tutorial data --- 
# data_path = SCimilarity_dir+os.sep+"GSE136831_subsample.h5ad"
# adams = sc.read(data_path) # the original data is log-normed
# adams = align_dataset(adams, ca.gene_order)
# adams = lognorm_counts(adams) # -> I do not understand why additional log-norm is important, but I kept this
# adams.obsm["X_scimilarity"] = ca.get_embeddings(adams.X)



# # --- let us try our data then ---
# Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'
# test_adata = utils.load_mtx(Chu_CD4T_dir+os.sep+'test_20_pancanT')
# test_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
# test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
# test_adata.layers['counts'] = test_adata.X.copy()  # save the original counts
# test_adata = align_dataset(test_adata, ca.gene_order)
# test_adata = lognorm_counts(test_adata)
# test_adata.obsm["X_scimilarity"] = ca.get_embeddings(test_adata.X)
# sc.pp.neighbors(test_adata, use_rep="X_scimilarity")
# sc.tl.umap(test_adata)
# sc.pl.umap(test_adata, color="celltype", legend_fontsize=5)
# plt.savefig(result_dir+os.sep+"tutorial/Chu_CD4T_test_umap.png", dpi=300, bbox_inches="tight")


# GSE179994_CD8_dir = data_dir+os.sep+'GSE179994_CD8T'
# train_adata = utils.load_mtx(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref')
# train_metadata = pd.read_csv(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref_metadata.csv', index_col=0)
# train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
# train_adata.layers['counts'] = train_adata.X.copy()  # save the original counts
# train_adata = align_dataset(train_adata, ca.gene_order)
# train_adata = lognorm_counts(train_adata)
# train_adata.obsm["X_scimilarity"] = ca.get_embeddings(train_adata.X)
# sc.pp.neighbors(train_adata, use_rep="X_scimilarity")
# sc.tl.umap(train_adata)
# sc.pl.umap(train_adata, color="celltype", legend_fontsize=5)
# plt.savefig(result_dir+os.sep+"tutorial/GSE179994_CD8T_train_umap.png", dpi=300, bbox_inches="tight")

# HNSC_CD8_dir = data_dir+os.sep+'HNSC_CD8T'
# test_adata = utils.load_mtx(HNSC_CD8_dir+os.sep+'HNSC_CD8T')
# test_metadata = pd.read_csv(HNSC_CD8_dir+os.sep+'HNSC_CD8T_metadata_curated.csv', index_col=0)
# test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
# test_adata.layers['counts'] = test_adata.X.copy()  # save the original counts
# test_adata = align_dataset(test_adata, ca.gene_order)
# test_adata = lognorm_counts(test_adata)
# test_adata.obsm["X_scimilarity"] = ca.get_embeddings(test_adata.X)
# sc.pp.neighbors(test_adata, use_rep="X_scimilarity")
# sc.tl.umap(test_adata)
# sc.pl.umap(test_adata, color="celltype", legend_fontsize=5)
# plt.savefig(result_dir+os.sep+"tutorial/HNSC_CD8T_test_umap.png", dpi=300, bbox_inches="tight")