# conda activate /net/mulan/home/wenjinma/miniconda3/envs/cellplm
import os, sys
import argparse
import numpy as np
import pandas as pd
import anndata
import scanpy as sc

import warnings
warnings.filterwarnings("ignore")

import hdf5plugin
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix
from CellPLM.utils import set_seed
from CellPLM.pipeline.cell_embedding import CellEmbeddingPipeline
from CellPLM.pipeline.cell_type_annotation import CellTypeAnnotationPipeline, CellTypeAnnotationDefaultPipelineConfig, CellTypeAnnotationDefaultModelConfig
# import rapids_singlecell as rsc

import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, f1_score, adjusted_rand_score

import data_utils

PRETRAIN_VERSION = '20230926_85M'
DEVICE = 'cuda:3'
set_seed(42)


# import os, sys 
# import argparse                                                                                                                                                                                                                                                                                                    
# from pathlib import Path
# import warnings
# warnings.filterwarnings("ignore", category=ResourceWarning)
# import anndata
# import scanpy as sc
# import scib
# import numpy as np
# import pandas as pd

# import scgpt as scg
# import scGPT_utils

# import matplotlib.pyplot as plt
# plt.style.context('default')
# from sklearn.cluster import KMeans
# # calculate ARI
# from sklearn.metrics import adjusted_rand_score

data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data/"
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD/"
result_dir = project_dir+os.sep+'results'+os.sep+'CellPLM'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# --- CellPLM model directory ---
CellPLM_dir = "/net/zootopia/disk1/wenjinma/LLM_models/CellPLM"

# add argparse
parser = argparse.ArgumentParser(description='Run scGPT pipeline for single-slice data.')
parser.add_argument('--experiment', type=str,
                    help='Load experiment related dataset.')
args = parser.parse_args()
experiment = args.experiment

CellPLM_result_dir = result_dir+os.sep+experiment
if not os.path.exists(CellPLM_result_dir):
    os.makedirs(CellPLM_result_dir)

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
    
# --- check embedding first ---
pipeline = CellEmbeddingPipeline(pretrain_prefix=PRETRAIN_VERSION, # Specify the pretrain checkpoint to load
                                 pretrain_directory=f'{CellPLM_dir}/ckpt')
pipeline.model

embedding = pipeline.predict(train_adata, # An AnnData object
                device=DEVICE) # Specify a gpu or cpu for model inference
train_adata.obsm['emb'] = embedding.cpu().numpy()
embedding = pipeline.predict(test_adata, # An AnnData object
                device=DEVICE) # Specify a gpu or cpu for model inference
test_adata.obsm['emb'] = embedding.cpu().numpy()
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
# Convert celltype to categorical and get the mapping
celltype_categorical = train_adata.obs[celltype_col].astype('category')
celltype_codes = celltype_categorical.cat.codes
celltype_categories = celltype_categorical.cat.categories
# Split the data into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(
    train_adata.obsm['emb'],
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
test_predictions_codes = xgb_classifier.predict(test_adata.obsm['emb'])
# Convert predictions back to categorical labels using the category mapping
test_predictions_labels = celltype_categories[test_predictions_codes]
# Save both the numeric codes and categorical labels to the test_adata.obs
test_adata.obs['predicted_celltype_code'] = test_predictions_codes  # [0,1,2,3] format
test_adata.obs['pred_celltype'] = test_predictions_labels       # original categorical labels
# generate confusion matrix
confusion_matrix_result = confusion_matrix(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'])
# --- calculate acc, macroF1, ARI
acc = accuracy_score(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'])
F1 = f1_score(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'], average='macro')
ARI = adjusted_rand_score(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'])
# write to file
result_df = pd.DataFrame({
    'Acc': [acc],
    'macroF1': [F1], 
    'ARI': [ARI]
})
result_df.to_csv(
    os.path.join(CellPLM_result_dir, 'CellPLM_zeroshot_results_metrics.csv'),
    index=False,  # Don't write row indices
    quoting=0     # Don't quote strings (equivalent to quote=F in R)
)
test_adata.obs.to_csv(CellPLM_result_dir+os.sep+'CellPLM_zeroshot_predicted_celltypes.csv')

# sc.pp.neighbors(data, use_rep='emb') # remove method='rapids' if rapids is not installed
# sc.tl.umap(data) # remove method='rapids' if rapids is not installed
# plt.rcParams['figure.figsize'] = (6, 6)
# sc.pl.umap(data, color='celltype', palette='Paired')
# plt.savefig(CellPLM_result_dir+os.sep+f'zeroshot_CellPLM_embedding.png', bbox_inches='tight')
# # sc.pl.umap(data, color='CancerType', palette='Paired')
# plt.savefig(CellPLM_result_dir+os.sep+f'zeroshot_CellPLM_embedding_cancertype.png', bbox_inches='tight')
# sc.pl.umap(data, color='TissueType', palette='Paired')
# plt.savefig(CellPLM_result_dir+os.sep+f'zeroshot_CellPLM_embedding_tissuetype.png', bbox_inches='tight')

# --- annotation
train_num = train_adata.shape[0]
data = ad.concat([train_adata, test_adata])
data.X = csr_matrix(data.X)
data.obs['celltype'] = data.obs[celltype_col]

data.obs['split'] = 'test'
tr = np.random.permutation(train_num) #torch.randperm(train_num).numpy()
data.obs['split'][tr[:int(train_num*0.9)]] = 'train'
data.obs['split'][tr[int(train_num*0.9):train_num]] = 'valid'

pipeline_config = CellTypeAnnotationDefaultPipelineConfig.copy()
pipeline_config['device'] = DEVICE

model_config = CellTypeAnnotationDefaultModelConfig.copy()
model_config['out_dim'] = data.obs['celltype'].nunique()
pipeline_config, model_config
pipeline = CellTypeAnnotationPipeline(pretrain_prefix=PRETRAIN_VERSION, # Specify the pretrain checkpoint to load
                                      overwrite_config=model_config,  # This is for overwriting part of the pretrain config
                                      pretrain_directory=f'{CellPLM_dir}/ckpt')
pipeline.fit(data, # An AnnData object
            pipeline_config, # The config dictionary we created previously, optional
            split_field = 'split', #  Specify a column in .obs that contains split information
            train_split = 'train',
            valid_split = 'valid',
            label_fields = ['celltype'],
            device=DEVICE) # Specify a column in .obs that contains cell type labels
pred_labels_codes = pipeline.predict(
                test_adata, # An AnnData object
                pipeline_config, # The config dictionary we created previously, optional
            )
pred_labels = pipeline.label_encoders['celltype'].inverse_transform(pred_labels_codes)
test_adata.obs['pred_celltype'] = pred_labels
# ari_score = adjusted_rand_score(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'])
# print(f"Finetune Adjusted Rand Index (ARI) score: {ari_score}")
acc = accuracy_score(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'])
F1 = f1_score(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'], average='macro')
ARI = adjusted_rand_score(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'])
# write to file
result_df = pd.DataFrame({
    'Acc': [acc],
    'macroF1': [F1], 
    'ARI': [ARI]
})
result_df.to_csv(
    os.path.join(CellPLM_result_dir, 'CellPLM_finetune_results_metrics.csv'),
    index=False,  # Don't write row indices
    quoting=0     # Don't quote strings (equivalent to quote=F in R)
)
test_adata.obs.to_csv(CellPLM_result_dir+os.sep+'CellPLM_finetune_predicted_celltypes.csv')

# pipeline.score(data, # An AnnData object
#                 pipeline_config, # The config dictionary we created previously, optional
#                 split_field = 'split', # Specify a column in .obs to specify train and valid split, optional
#                 target_split = 'test', # Specify a target split to predict, optional
#                 label_fields = ['celltype'])  # Specify a column in .obs that contains cell type labels
