# conda activate /net/zootopia/disk1/wenjinma/envs/CellFM
import os, sys
import argparse
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from tqdm import tqdm,trange
import mindspore as ms
from mindspore.train import Model, CheckpointConfig, ModelCheckpoint, LossMonitor
from scipy.sparse import csr_matrix as csm

import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score

import data_utils



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
result_dir = project_dir+os.sep+'results'+os.sep+'CellFM'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# --- CellFM model directory ---
CellFM_dir = "/net/zootopia/disk1/wenjinma/LLM_models/CellFM"
sys.path.append(f'{CellFM_dir}/CellFM')

ms.set_context(
    device_target='CPU', 
    mode=ms.GRAPH_MODE,
    device_id=0,
)
ms.set_seed(0)

from config import Config
from annotation_model import *
from metrics import annote_metric
from data_process import Prepare
import CellFM_annotation_utils
from utils import Wrapper

# --- process gene
gene_info=pd.read_csv(f'{CellFM_dir}/CellFM/csv/gene_info.csv',header=0,index_col=0)
geneset=gene_info.index
genemap={j:i+1 for i,j in enumerate(gene_info.index)}
hgcn=pd.read_csv(f'{CellFM_dir}/CellFM/csv/updated_hgcn.tsv',index_col=1,header=0,sep='\t')
hgcn=hgcn[hgcn['Status']=='Approved']
map_dict={}
alias=hgcn['Alias symbols']
prev=hgcn['Previous symbols']
for i in hgcn.index:
    if alias.loc[i] is not np.nan:
        for j in alias.loc[i].split(', '):
            if j not in hgcn.index:
                map_dict[j]=i
for i in hgcn.index:
    if prev.loc[i] is not np.nan:
        for j in prev.loc[i].split(', '):
            if j not in hgcn.index:
                map_dict[j]=i
egsn=pd.read_csv(f'{CellFM_dir}/CellFM/csv/updated_hgcn.tsv',index_col=None,header=0,sep='\t')
egsn=egsn.dropna(subset=['Ensembl gene ID'])
egsn=egsn.set_index('Ensembl gene ID')


# add argparse
parser = argparse.ArgumentParser(description='Run scGPT pipeline for single-slice data.')
parser.add_argument('--experiment', type=str,
                    help='Load experiment related dataset.')
args = parser.parse_args()
experiment = args.experiment

CellFM_result_dir = result_dir+os.sep+experiment
if not os.path.exists(CellFM_result_dir):
    os.makedirs(CellFM_result_dir)

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
    celltyoe_col = 'celltype'
    
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
    

train_adata.obs['cell_type'] = train_adata.obs[celltype_col].values
test_adata.obs['cell_type'] = test_adata.obs[celltype_col].values
train_adata.obs['train'] = 0
test_adata.obs['train']  = 2
adata = anndata.concat([train_adata, test_adata], join='outer')
print('origin shape:',adata.shape,len(adata.obs['cell_type'].unique()))

data=adata.X.astype(np.float32)
T=adata.X.sum(1)
data=csm(np.round(data/np.maximum(1,T/1e5,dtype=np.float32)))
data.eliminate_zeros()
adata.X=data

trainset = CellFM_annotation_utils.SCrna(adata, CellFM_dir, mode='train')
testset = CellFM_annotation_utils.SCrna(adata, CellFM_dir, mode='test')

cfg=Config()
cfg.num_cls=trainset.cls
cfg.enc_nlayers=2

prep=Prepare(
    cfg.nonz_len,pad=1,mask_ratio=0,random=False
)
train_loader=CellFM_annotation_utils.build_dataset(
    trainset,
    prep,
    128,
    drop=True,
    shuffle=True,
)
test_loader=CellFM_annotation_utils.build_dataset(
    testset,
    prep,
    1,
    drop=False,
    shuffle=False,
)
para=ms.load_checkpoint(f"{CellFM_dir}/CellFM_80M_weight.ckpt")
backbone=Backbone(len(trainset.geneset),cfg)
ms.load_param_into_net(backbone, para)
model=Net(backbone,cfg)
CellFM_annotation_utils.freeze_module(model.extractor)

optimizer=nn.Adam(model.trainable_params(),1e-4,weight_decay=1e-5)
update_cell=nn.DynamicLossScaleUpdateCell(1,2,1000)
wrapper=Wrapper(model,optimizer)
trainer=Model(
    wrapper,
    eval_network=model,
    amp_level='O0',
    metrics={
        'accuracy':annote_metric(trainset.cls,key='accuracy'),
    },
    eval_indexes=[0,1,2]
)
loss_cb = LossMonitor(20)

ckpt_config = CheckpointConfig(
    save_checkpoint_steps=len(train_loader),
    keep_checkpoint_max=1,
    integrated_save=False,
    async_save=False
)
ckpt_cb = ModelCheckpoint(
    prefix=f'{experiment}_zeroshot', 
    directory=f"{CellFM_result_dir}", 
    config=ckpt_config
)
cbs=[loss_cb,ckpt_cb]

trainer.train(30,train_loader,callbacks=cbs)
ms.load_param_into_net(model, ms.load_checkpoint(ckpt_cb.latest_ckpt_file_name))
trainer.eval(test_loader)
