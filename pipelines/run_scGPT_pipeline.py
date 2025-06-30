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

import matplotlib.pyplot as plt
plt.style.context('default')
from sklearn.cluster import KMeans
# calculate ARI
from sklearn.metrics import adjusted_rand_score

data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data/"
project_dir = "/net/mulan/home/wenjinma/projects/scPanKD/"
result_dir = project_dir+os.sep+'results'+os.sep+'scGPT'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# --- scGPT model directory ---
scGPT_dir = "/net/zootopia/disk1/wenjinma/LLM_models/scGPT"
# scGPT_model_dir = scGPT_dir+os.sep+'scGPT_CP'
scGPT_model_dir = scGPT_dir+os.sep+'scGPT_pancancer'

# add argparse
parser = argparse.ArgumentParser(description='Run scGPT pipeline for single-slice data.')
parser.add_argument('--experiment', type=str,
                    help='Load experiment related dataset.')
args = parser.parse_args()
experiment = args.experiment

scGPT_result_dir = result_dir+os.sep+experiment
if not os.path.exists(scGPT_result_dir):
    os.makedirs(scGPT_result_dir)

if experiment == 'Chu_CD4T_validation':
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'
    train_adata = scGPT_utils.load_mtx(Chu_CD4T_dir+os.sep+'train_80_pancanT')
    train_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'train_80_pancanT_metadata.csv', index_col=0)
    train_adata.obs = train_adata.obs.merge(train_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))
    test_adata = scGPT_utils.load_mtx(Chu_CD4T_dir+os.sep+'test_20_pancanT')
    test_metadata = pd.read_csv(Chu_CD4T_dir+os.sep+'test_20_pancanT_metadata.csv', index_col=0)
    test_adata.obs = test_adata.obs.merge(test_metadata, left_on="barcode", right_index=True, how='left', suffixes=('', '_x'))  
    wandb_config = scGPT_utils.init_wandb_config(scGPT_model_dir)
    scGPT_utils.run_scGPT_finetune(wandb_config, train_adata, test_adata, scGPT_result_dir,
                                   train_celltype_col='curated_celltype', device='cuda:2')
