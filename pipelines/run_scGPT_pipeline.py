# conda activate ~/envs/scGPT_bak/
import os, sys 
import argparse                                                                                                                                                                                                                                                                                                    
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=ResourceWarning)
import anndata
import scanpy as sc
import scib
import numpy as np
import pandas as pd

import scgpt as scg
import scGPT_utils
import data_utils

import matplotlib.pyplot as plt
plt.style.context('default')
# calculate ARI
from sklearn.metrics import adjusted_rand_score, accuracy_score, f1_score

data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data/"
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD/"
result_dir = project_dir+os.sep+'results'+os.sep+'scGPT'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# --- scGPT model directory ---
scGPT_dir = "/net/zootopia/disk1/wenjinma/LLM_models/scGPT"
# scGPT_model_dir = scGPT_dir+os.sep+'scGPT_CP'
scGPT_model_dir = scGPT_dir+os.sep+'scGPT_pancancer' # let use use pancancer dataset

# add argparse
parser = argparse.ArgumentParser(description='Run scGPT pipeline for single-slice data.')
parser.add_argument('--experiment', type=str,
                    help='Load experiment related dataset.')
args = parser.parse_args()
experiment = args.experiment

scGPT_result_dir = result_dir+os.sep+experiment
if not os.path.exists(scGPT_result_dir):
    os.makedirs(scGPT_result_dir)

if experiment == 'Chu_CD4_multibatch_validation':
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'
    train_adata = data_utils.load_mtx(Chu_CD4T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata = data_utils.load_mtx(Chu_CD4T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))  
    celltype_col = 'curated_celltype'

if experiment == 'Chu_CD8_multibatch_validation':
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    train_adata = data_utils.load_mtx(Chu_CD8T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata = data_utils.load_mtx(Chu_CD8T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))  
    test_adata.obs['celltype'] = test_adata.obs['curated_celltype'].values
    celltype_col = 'celltype'

if experiment == 'ProjecTILs_CD8_multibatch_to_GSE179994_CD8':
    ProjecTILs_CD8_dir = data_dir+os.sep+'ProjecTILs_CD8T'
    train_adata = data_utils.load_mtx(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T')
    train_metadata = pd.read_csv(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    GSE179994_CD8_dir = data_dir+os.sep+'GSE179994_CD8T'
    test_adata = data_utils.load_mtx(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref')
    test_metadata = pd.read_csv(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'

if experiment == 'Chu_CD8_multibatch_to_GSE179994_CD8':
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    train_adata = data_utils.load_mtx(Chu_CD8T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    GSE179994_CD8_dir = data_dir+os.sep+'GSE179994_CD8T'
    test_adata = data_utils.load_mtx(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref')
    test_metadata = pd.read_csv(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'

if experiment == 'ProjecTILs_CD8_multibatch_to_HNSC_CD8':
    ProjecTILs_CD8_dir = data_dir+os.sep+'ProjecTILs_CD8T'
    train_adata = data_utils.load_mtx(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T')
    train_metadata = pd.read_csv(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    HNSC_CD8_dir = data_dir+os.sep+'HNSC_CD8T'
    test_adata = data_utils.load_mtx(HNSC_CD8_dir+os.sep+'HNSC_CD8T')
    test_metadata = pd.read_csv(HNSC_CD8_dir+os.sep+'HNSC_CD8T_metadata_curated.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'

if experiment == 'Chu_CD8_multibatch_to_HNSC_CD8':
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    train_adata = data_utils.load_mtx(Chu_CD8T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    HNSC_CD8_dir = data_dir+os.sep+'HNSC_CD8T'
    test_adata = data_utils.load_mtx(HNSC_CD8_dir+os.sep+'HNSC_CD8T')
    test_metadata = pd.read_csv(HNSC_CD8_dir+os.sep+'HNSC_CD8T_metadata_curated.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'


if experiment == 'ProjecTILs_CD8_multibatch_to_Chu_CD8':
    ProjecTILs_CD8_dir = data_dir+os.sep+'ProjecTILs_CD8T'
    train_adata = data_utils.load_mtx(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T')
    train_metadata = pd.read_csv(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    test_adata = data_utils.load_mtx(Chu_CD8T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata.obs['celltype'] = test_adata.obs['curated_celltype'].values
    celltype_col = 'celltype'

if experiment == 'Chu_CD8_multibatch_to_ProjecTILs_CD8':
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    train_adata = data_utils.load_mtx(Chu_CD8T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x')) 
    ProjecTILs_CD8_dir = data_dir+os.sep+'ProjecTILs_CD8T'
    test_adata = data_utils.load_mtx(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T')
    test_metadata = pd.read_csv(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'

if experiment == 'Zheng_CD4_multibatch_to_Chu_CD4':
    Zheng_CD4T_dir = data_dir+os.sep+'Zheng_CD4T'
    train_adata = data_utils.load_mtx(Zheng_CD4T_dir+os.sep+'Zheng_CD4T')
    train_metadata = pd.read_csv(Zheng_CD4T_dir+os.sep+'Zheng_CD4T_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'
    test_adata = data_utils.load_mtx(Chu_CD4T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'

if experiment == 'Chu_CD4_multibatch_to_Zheng_CD4':
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'
    train_adata = data_utils.load_mtx(Chu_CD4T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    Zheng_CD4T_dir = data_dir+os.sep+'Zheng_CD4T'
    test_adata = data_utils.load_mtx(Zheng_CD4T_dir+os.sep+'Zheng_CD4T')
    test_metadata = pd.read_csv(Zheng_CD4T_dir+os.sep+'Zheng_CD4T_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'


# --- scGPT embeddings ---
ref_embed_adata = scg.tasks.embed_data(
    train_adata,
    scGPT_model_dir,
    gene_col='genes',
    batch_size=64,
)
test_embed_adata = scg.tasks.embed_data(
    test_adata,
    scGPT_model_dir,
    gene_col='genes',
    batch_size=64,
)
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
# Convert celltype to categorical and get the mapping
celltype_categorical = ref_embed_adata.obs[celltype_col].astype('category')
celltype_codes = celltype_categorical.cat.codes
celltype_categories = celltype_categorical.cat.categories

# Split the data into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(
    ref_embed_adata.obsm['X_scGPT'],
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
test_predictions_codes = xgb_classifier.predict(test_embed_adata.obsm['X_scGPT'])
# Convert predictions back to categorical labels using the category mapping
test_predictions_labels = celltype_categories[test_predictions_codes]
# Save both the numeric codes and categorical labels to the test_adata.obs
test_embed_adata.obs['predicted_celltype_code'] = test_predictions_codes  # [0,1,2,3] format
test_embed_adata.obs['pred_celltype'] = test_predictions_labels       # original categorical labels
# generate confusion matrix
confusion_matrix_result = confusion_matrix(test_embed_adata.obs[celltype_col], test_embed_adata.obs['pred_celltype'])
# --- calculate acc, macroF1, ARI
acc = accuracy_score(test_embed_adata.obs[celltype_col], test_embed_adata.obs['pred_celltype'])
F1 = f1_score(test_embed_adata.obs[celltype_col], test_embed_adata.obs['pred_celltype'], average='macro')
ARI = adjusted_rand_score(test_embed_adata.obs[celltype_col], test_embed_adata.obs['pred_celltype'])
# write to file
result_df = pd.DataFrame({
    'Acc': [acc],
    'macroF1': [F1], 
    'ARI': [ARI]
})
result_df.to_csv(
    os.path.join(scGPT_result_dir, 'scGPT_zeroshot_results_metrics.csv'),
    index=False,  # Don't write row indices
    quoting=0     # Don't quote strings (equivalent to quote=F in R)
)
test_embed_adata.obs.to_csv(scGPT_result_dir+os.sep+'scGPT_zeroshot_predicted_celltypes.csv')

# --- scGPT finetune ---
wandb_config = scGPT_utils.init_wandb_config(scGPT_model_dir)
scGPT_utils.run_scGPT_finetune(wandb_config, train_adata, test_adata, scGPT_result_dir,
                                train_celltype_col=celltype_col)
# load finetuned results
directories = os.listdir(scGPT_result_dir+os.sep+'save')
experiment_dir = scGPT_result_dir+os.sep+'save'+os.sep+directories[0]
# load pkl file
import pickle
if os.path.exists(experiment_dir+os.sep+'results.pkl'):
    with open(experiment_dir+os.sep+'results.pkl', "rb") as f:
        scGPT_finetune_results = pickle.load(f)
    acc = accuracy_score(scGPT_finetune_results['predictions'], scGPT_finetune_results['labels'])
    F1 = f1_score(scGPT_finetune_results['predictions'], scGPT_finetune_results['labels'], average='macro')
    ARI = adjusted_rand_score(scGPT_finetune_results['predictions'], scGPT_finetune_results['labels'])
    # write to file
    result_df = pd.DataFrame({
        'Acc': [acc],
        'macroF1': [F1], 
        'ARI': [ARI]
    })
    result_df.to_csv(
        os.path.join(scGPT_result_dir, 'scGPT_finetune_results_metrics.csv'),
        index=False,  # Don't write row indices
        quoting=0     # Don't quote strings (equivalent to quote=F in R)
    )
