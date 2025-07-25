# conda activate /net/zootopia/disk1/wenjinma/envs/GeneFormer
import os, sys
os.environ["CUDA_VISIBLE_DEVICES"] = "3"
import argparse
import numpy as np
import pandas as pd
import anndata
import scanpy as sc

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import anndata as ad
from geneformer import EmbExtractor
from geneformer import TranscriptomeTokenizer

import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score
import data_utils
import mygene # used to convert gene

data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data/"
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD/"
result_dir = project_dir+os.sep+'results'+os.sep+'GeneFormer'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# --- GeneFormer model directory ---
GeneFormer_dir = '/net/zootopia/disk1/wenjinma/LLM_models/Geneformer'
# GeneFormer_model_dir = GeneFormer_dir+os.sep+'Geneformer-V2-104M'
GeneFormer_model_dir = GeneFormer_dir+os.sep+'Geneformer-V2-316M'

# add argparse
parser = argparse.ArgumentParser(description='Run GeneFormer pipeline for single-slice data.')
parser.add_argument('--experiment', type=str,
                    help='Load experiment related dataset.')
args = parser.parse_args()
experiment = args.experiment

GeneFormer_result_dir = result_dir+os.sep+experiment
if not os.path.exists(GeneFormer_result_dir):
    os.makedirs(GeneFormer_result_dir)

if experiment == 'Chu_CD4T_validation':
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'
    train_adata = data_utils.load_mtx(Chu_CD4T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata = data_utils.load_mtx(Chu_CD4T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))      
    celltype_col = 'curated_celltype'

if experiment == 'Chu_CD8T_validation':
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    train_adata = data_utils.load_mtx(Chu_CD8T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata = data_utils.load_mtx(Chu_CD8T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))  
    test_adata.obs['celltype'] = test_adata.obs['curated_celltype'].values
    celltype_col = 'celltype'

if experiment == 'Zheng_CD4_to_Chu_CD4':
    Zheng_CD4_dir = data_dir+os.sep+'Zheng_CD4T'
    train_adata = data_utils.load_mtx(Zheng_CD4_dir+os.sep+'Zheng_CD4T')
    train_metadata = pd.read_csv(Zheng_CD4_dir+os.sep+'Zheng_CD4T_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'
    test_adata = data_utils.load_mtx(Chu_CD4T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'

if experiment == 'GSE179994_CD8_to_HNSC_CD8':
    GSE179994_CD8_dir = data_dir+os.sep+'GSE179994_CD8T'
    train_adata = data_utils.load_mtx(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref')
    train_metadata = pd.read_csv(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    HNSC_CD8_dir = data_dir+os.sep+'HNSC_CD8T'
    test_adata = data_utils.load_mtx(HNSC_CD8_dir+os.sep+'HNSC_CD8T')
    test_metadata = pd.read_csv(HNSC_CD8_dir+os.sep+'HNSC_CD8T_metadata_curated.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'

if experiment == 'ProjecTILs_CD8_to_HNSC_CD8':
    ProjecTILs_CD8_dir = data_dir+os.sep+'ProjecTILs_CD8T'
    train_adata = data_utils.load_mtx(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T')
    train_metadata = pd.read_csv(ProjecTILs_CD8_dir+os.sep+'ProjecTIL_CD8T_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    HNSC_CD8_dir = data_dir+os.sep+'HNSC_CD8T'
    test_adata = data_utils.load_mtx(HNSC_CD8_dir+os.sep+'HNSC_CD8T')
    test_metadata = pd.read_csv(HNSC_CD8_dir+os.sep+'HNSC_CD8T_metadata_curated.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    celltype_col = 'celltype'

if experiment == 'GSE179994_CD8_to_Chu_CD8T':
    GSE179994_CD8_dir = data_dir+os.sep+'GSE179994_CD8T'
    train_adata = data_utils.load_mtx(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref')
    train_metadata = pd.read_csv(GSE179994_CD8_dir+os.sep+'GSE179994_CD8T_ref_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'
    test_adata = data_utils.load_mtx(Chu_CD8T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD8T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata.obs['celltype'] = test_adata.obs['curated_celltype'].values
    celltype_col = 'celltype'  

if experiment == 'ProjecTILs_CD8_to_Chu_CD8T':
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
    
# --- check embedding first ---
mg = mygene.MyGeneInfo()
mapped_train_out = mg.querymany(train_adata.var_names.to_list(), scopes='symbol', fields='ensembl.gene', species='human')
mapped_test_out = mg.querymany(test_adata.var_names.to_list(), scopes='symbol', fields='ensembl.gene', species='human')

def map_genes(genelist, mapped_out):
    """
    Map genes in adata.var_names to ensembl IDs using the mapped_out from mygene.
    """
    gene_to_ensembl = {}
    for gene, result in zip(genelist, mapped_out):
        if 'ensembl' in result and 'notfound' not in result:
            if isinstance(result['ensembl'], list):
                # Take first Ensembl ID
                gene_to_ensembl[gene] = result['ensembl'][0]['gene']
            else:
                gene_to_ensembl[gene] = result['ensembl']['gene']
        else:
            gene_to_ensembl[gene] = np.nan
    return gene_to_ensembl

train_gene_to_ensembl = map_genes(train_adata.var_names.to_list(), mapped_train_out)
test_gene_to_ensembl = map_genes(test_adata.var_names.to_list(), mapped_test_out)

train_adata.var['ensembl_id'] = [train_gene_to_ensembl.get(gene, np.nan) for gene in train_adata.var_names.to_list()]
test_adata.var['ensembl_id'] = [test_gene_to_ensembl.get(gene, np.nan) for gene in test_adata.var_names.to_list()]
# --- remove NAN and duplicate ensemble IDs ---
train_adata_mapped = train_adata[:, train_adata.var['ensembl_id'].notna()].copy()
train_adata_mapped = train_adata_mapped[:, ~train_adata_mapped.var['ensembl_id'].duplicated(keep='first')].copy()
test_adata_mapped = test_adata[:, test_adata.var['ensembl_id'].notna()].copy()
test_adata_mapped = test_adata_mapped[:, ~test_adata_mapped.var['ensembl_id'].duplicated(keep='first')].copy()

# mapped genes
train_adata_mapped.obs['n_counts'] = train_adata_mapped.X.sum(axis=1).A1
test_adata_mapped.obs['n_counts'] = test_adata_mapped.X.sum(axis=1).A1
# create new directories for tokenization and embedding
GeneFormer_train_dir = GeneFormer_result_dir+os.sep+"train"
GeneFormer_test_dir = GeneFormer_result_dir+os.sep+"test"
if not os.path.exists(GeneFormer_train_dir):
    os.makedirs(GeneFormer_train_dir)
if not os.path.exists(GeneFormer_test_dir):
    os.makedirs(GeneFormer_test_dir)

train_adata_mapped.write_h5ad(GeneFormer_train_dir+os.sep+"train.h5ad")
test_adata_mapped.write_h5ad(GeneFormer_test_dir+os.sep+"test.h5ad")
tk = TranscriptomeTokenizer({'barcode': 'barcode'}, nproc=16)  # for V1 model, set model_version="V1"
tk.tokenize_data(GeneFormer_train_dir, 
                GeneFormer_train_dir, 
                "transformed", 
                file_format="h5ad")
tk.tokenize_data(GeneFormer_test_dir, 
                GeneFormer_test_dir, 
                "transformed", 
                file_format="h5ad")
embex = EmbExtractor(model_type="Pretrained",
                    emb_mode='cls',
                    model_version='V2',
                    max_ncells=None,
                    emb_label=['barcode'],
                    nproc=16)
train_embeddings = embex.extract_embs(f'{GeneFormer_model_dir}',
                 f'{GeneFormer_train_dir}/transformed.dataset', GeneFormer_train_dir, 'embeddings')
train_embeddings.index = train_embeddings['barcode'].values
train_embeddings.drop(columns=['barcode'], inplace=True, errors='ignore')
train_embeddings = train_embeddings.loc[train_adata_mapped.obs_names, :]
train_adata_mapped.obsm['X_GeneFormer'] = train_embeddings.values
train_adata_mapped.write_h5ad(GeneFormer_result_dir+os.sep+"GeneFormer_train_embedding.h5ad")

test_embeddings = embex.extract_embs(f'{GeneFormer_model_dir}',
                 f'{GeneFormer_test_dir}/transformed.dataset', GeneFormer_test_dir, 'embeddings')
test_embeddings.index = test_embeddings['barcode'].values
test_embeddings.drop(columns=['barcode'], inplace=True, errors='ignore')
test_embeddings = test_embeddings.loc[test_adata_mapped.obs_names, :]
# --- save embeddings to adata
test_adata_mapped.obsm['X_GeneFormer'] = test_embeddings.values
test_adata_mapped.write_h5ad(GeneFormer_result_dir+os.sep+"GeneFormer_test_embedding.h5ad")
# --- visualize
sc.pp.neighbors(train_adata_mapped, use_rep='X_GeneFormer')
sc.tl.umap(train_adata_mapped)
# sc.pl.umap(train_adata_mapped, color=[celltype_col], palette='Paired')
sc.pl.umap(train_adata_mapped, color=[celltype_col, 'TissueType', 'CancerType'], palette='Paired')
plt.savefig(GeneFormer_result_dir+os.sep+f'GeneFormer_train_embedding.png', bbox_inches='tight')
plt.close()

sc.pp.neighbors(test_adata_mapped, use_rep='X_GeneFormer')
sc.tl.umap(test_adata_mapped)
sc.pl.umap(test_adata_mapped, color=[celltype_col], palette='Paired')
plt.savefig(GeneFormer_result_dir+os.sep+f'GeneFormer_test_embedding.png', bbox_inches='tight')
plt.close()

from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
# Convert celltype to categorical and get the mapping
celltype_categorical = train_adata_mapped.obs[celltype_col].astype('category')
celltype_codes = celltype_categorical.cat.codes
celltype_categories = celltype_categorical.cat.categories
# Split the data into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(
    train_adata_mapped.obsm['X_GeneFormer'],
    celltype_codes,  # Use the categorical codes
    test_size=0.2,
    random_state=42,
    stratify=train_adata_mapped.obs[celltype_col]
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
test_predictions_codes = xgb_classifier.predict(test_adata_mapped.obsm['X_GeneFormer'])
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
test_adata.obs.to_csv(GeneFormer_result_dir+os.sep+'GeneFormer_zeroshot_predicted_celltypes.csv')


# --- finetune: GeneFormer finetune is not very friendly: https://huggingface.co/ctheodoris/Geneformer/blob/main/examples/cell_classification.ipynb
# train_num = train_adata.shape[0]
# data = ad.concat([train_adata, test_adata])
# data.X = csr_matrix(data.X)
# data.obs['celltype'] = data.obs[celltype_col]

# data.obs['split'] = 'test'
# tr = np.random.permutation(train_num) #torch.randperm(train_num).numpy()
# data.obs['split'][tr[:int(train_num*0.9)]] = 'train'
# data.obs['split'][tr[int(train_num*0.9):train_num]] = 'valid'

# pipeline_config = CellTypeAnnotationDefaultPipelineConfig.copy()
# pipeline_config['device'] = DEVICE

# model_config = CellTypeAnnotationDefaultModelConfig.copy()
# model_config['out_dim'] = data.obs['celltype'].nunique()
# pipeline_config, model_config
# pipeline = CellTypeAnnotationPipeline(pretrain_prefix=PRETRAIN_VERSION, # Specify the pretrain checkpoint to load
#                                       overwrite_config=model_config,  # This is for overwriting part of the pretrain config
#                                       pretrain_directory=f'{CellPLM_dir}/ckpt')
# pipeline.fit(data, # An AnnData object
#             pipeline_config, # The config dictionary we created previously, optional
#             split_field = 'split', #  Specify a column in .obs that contains split information
#             train_split = 'train',
#             valid_split = 'valid',
#             label_fields = ['celltype'],
#             device=DEVICE) # Specify a column in .obs that contains cell type labels
# pred_labels_codes = pipeline.predict(
#                 test_adata, # An AnnData object
#                 pipeline_config, # The config dictionary we created previously, optional
#             )
# pred_labels = pipeline.label_encoders['celltype'].inverse_transform(pred_labels_codes)
# test_adata.obs['pred_celltype'] = pred_labels
# ari_score = adjusted_rand_score(test_adata.obs[celltype_col], test_adata.obs['pred_celltype'])
# print(f"Finetune Adjusted Rand Index (ARI) score: {ari_score}")
# test_adata.obs.to_csv(CellPLM_result_dir+os.sep+'CellPLM_finetune_predicted_celltypes.csv')

# # pipeline.score(data, # An AnnData object
# #                 pipeline_config, # The config dictionary we created previously, optional
# #                 split_field = 'split', # Specify a column in .obs to specify train and valid split, optional
# #                 target_split = 'test', # Specify a target split to predict, optional
# #                 label_fields = ['celltype'])  # Specify a column in .obs that contains cell type labels
