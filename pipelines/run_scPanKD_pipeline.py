import os, glob, sys
import argparse
import torch

sys.path.append('/net/mulan/home/wenjinma/projects/scPanKD/scPanKD_code/scPanKD')
# import my packages
import train
import predict
from utils import _utils
_utils._set_seed(1993) # set seed for reproducibility

project_dir = "/net/mulan/home/wenjinma/projects/scPanKD"
data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data"
result_dir = project_dir+os.sep+'results/scPanKD'
parser = argparse.ArgumentParser()
parser.add_argument('--experiment', type=str, help='Experiments to run')
args = parser.parse_args()
print(args)

output_dir = result_dir+os.sep+args.experiment
os.makedirs(output_dir, exist_ok=True)
# --- train configure
model_config = {
    'fs':'F-test', # 'fs_strategy': 'normal',# F-test is better than using Seurat HVG
    'num_features': 1000, 'optimizer': 'adamW',
    'learning_rate': 0.001, 'batch_size': 32, 'max_epochs': 100, 
    "distillation_epochs": 30,
    'MLP_DIMS': [256, 64, 16], 'teacher_MLP_DIMS': [64, 16], 'student_MLP_DIMS': [64, 16],
    'Celltype_COLUMN': "celltype", 'PredCelltype_COLUMN': "pred_celltype",
    "entropy_quantile": 0.3}
# --- write config to file
config_file = output_dir+os.sep+'model_config.txt'
with open(config_file, 'w') as f:
    for key, value in model_config.items():
        f.write(f"{key}: {value}\n")
model_config['device'] = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')

# --- ProjecTILs_CD8T to HNSC_CD8T -> each
if args.experiment == 'ProjecTILs_each_CD8_to_HNSC_CD8':
    ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'trimmed'
    studies = ['EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 
               'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464']
    for study in studies:
        each_result_dir = output_dir+os.sep+study
        os.makedirs(each_result_dir, exist_ok=True)
        # --- train configure
        train_input_dict = {'batch_dir': ProjecTIL_dir,
                            'input': [f'ProjecTIL_{study}_CD8T_trimmed'],
                            'batch_info': [0],
                            'metadata': ProjecTIL_dir+os.sep+f'ProjecTIL_{study}_CD8T_trimmed_metadata.csv',
                            'output_dir': each_result_dir,
                            'prefix': 'train_'}
        train_args = argparse.Namespace(**train_input_dict)
        model = train.train_batch_MLP(train_args, model_config)
        # --- test configure
        test_input_dict = {'trained_model': each_result_dir+os.sep+'train_MLP_model',
                'input': data_dir+os.sep+'HNSC_CD8T'+os.sep+'HNSC_CD8T',
                'output_dir': each_result_dir,
                'oneround': False, 'prefix': 'scPanKD_predicted_'}
        predict_args = argparse.Namespace(**test_input_dict)
        predict.predict(predict_args, model, model_config)

if args.experiment == 'ProjecTILs_each_CD8_cancertrimmed_to_HNSC_CD8':
    ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'Tissue_trimmed'
    studies = ['HNSCC', 'Melanoma', 'SCC', 'Lung', 'Endometrial', 'Renal', 'Breast']
    for study in studies:
        each_result_dir = output_dir+os.sep+study
        os.makedirs(each_result_dir, exist_ok=True)
        # --- train configure
        train_input_dict = {'batch_dir': ProjecTIL_dir,
                            'input': [f'ProjecTIL_{study}_CD8T_trimmed'],
                            'batch_info': [0],
                            'metadata': ProjecTIL_dir+os.sep+f'ProjecTIL_{study}_CD8T_trimmed_metadata.csv',
                            'output_dir': each_result_dir,
                            'prefix': 'train_'}
        train_args = argparse.Namespace(**train_input_dict)
        model = train.train_batch_MLP(train_args, model_config)
        # --- test configure
        test_input_dict = {'trained_model': each_result_dir+os.sep+'train_MLP_model',
                'input': data_dir+os.sep+'HNSC_CD8T'+os.sep+'HNSC_CD8T',
                'output_dir': each_result_dir,
                'oneround': False, 'prefix': 'scPanKD_predicted_'}
        predict_args = argparse.Namespace(**test_input_dict)
        predict.predict(predict_args, model, model_config)

# --- ProjecTILs_CD8T to HNSC_CD8T -> batch effect counted
if args.experiment == 'ProjecTILs_CD8_multibatch_trimmed_to_HNSC_CD8':
    ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'trimmed'
    studies = ['EGAS00001004809', 'GSE123814', 'GSE139555', 'GSE159251', 
               'GSE176021', 'GSE179994', 'GSE180268', 'PRJNA705464']
    studies_input = []
    for study in studies:
        studies_input.append(f'ProjecTIL_{study}_CD8T_trimmed')
    batch_info = list(range(len(studies_input)))
    # --- train configure
    train_input_dict = {'batch_dir': ProjecTIL_dir,
                        'input': studies_input,
                        'batch_info': batch_info,
                        'metadata': os.path.dirname(ProjecTIL_dir)+os.sep+'ProjecTIL_CD8T_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'HNSC_CD8T'+os.sep+'HNSC_CD8T',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'ProjecTILs_CD8_multibatch_cancertrimmed_to_HNSC_CD8':
    ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'Tissue_trimmed'
    studies = ['HNSCC', 'Melanoma', 'SCC', 'Lung', 'Endometrial', 'Renal', 'Breast']
    studies_input = []
    for study in studies:
        studies_input.append(f'ProjecTIL_{study}_CD8T_trimmed')
    batch_info = list(range(len(studies_input)))
    # --- train configure
    train_input_dict = {'batch_dir': ProjecTIL_dir,
                        'input': studies_input,
                        'batch_info': batch_info,
                        'metadata': os.path.dirname(ProjecTIL_dir)+os.sep+'ProjecTIL_CD8T_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'HNSC_CD8T'+os.sep+'HNSC_CD8T',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    student = predict.predict(predict_args, model, model_config)
    # --- analyze data
    import anndata
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import matplotlib.pyplot as plt
    # train_adata = anndata.read_h5ad(output_dir+os.sep+'train_adata.h5ad')
    # x_train = _utils._extract_adata(train_adata)
    # _, y_pred_features = student.predict(torch.tensor(x_train, dtype=torch.float32, device=model_config['device']))
    # train_adata.obsm['embedding'] = y_pred_features.detach().cpu().numpy()
    # train_adata.write_h5ad(output_dir+os.sep+predict_args.prefix+'reference_features.h5ad')
    train_adata = anndata.read_h5ad(output_dir+os.sep+predict_args.prefix+'reference_features.h5ad')
    test_adata = anndata.read_h5ad(output_dir+os.sep+predict_args.prefix+'query_features.h5ad')
    test_pred_res = pd.read_csv(output_dir+os.sep+'scPanKD_predicted_celltypes-best.csv', index_col=0)
    test_adata.obs = test_pred_res.loc[test_adata.obs_names, :]
    embedding_df = np.concatenate((train_adata.obsm['embedding'], test_adata.X), axis=0)
    embedding_adata = anndata.AnnData(X=embedding_df, obs=train_adata.obs.append(test_adata.obs))
    # --- perform UMAP
    sc.pp.neighbors(embedding_adata, use_rep='X')
    sc.tl.umap(embedding_adata)
    embedding_adata.write_h5ad(output_dir+os.sep+'joint_embedding_adata.h5ad')
    # --- visualize
    train_embedding_adata = embedding_adata[train_adata.obs_names, :]
    test_embedding_adata = embedding_adata[test_adata.obs_names, :]
    with plt.rc_context({"figure.figsize": (5, 3)}):
        sc.pl.umap(train_embedding_adata, color=['Tissue', 'Cohort'], wspace=0.5)
        plt.savefig(output_dir+os.sep+'train_umap_Tissue_Cancertype.pdf', bbox_inches='tight')
    with plt.rc_context({"figure.figsize": (6, 3)}):
        sc.pl.umap(test_embedding_adata, color=['pred_celltype'])
        plt.savefig(output_dir+os.sep+'test_umap_Predcelltype.pdf', bbox_inches='tight')
    with plt.rc_context({"figure.figsize": (6, 3)}):
        sc.pl.umap(train_embedding_adata, color=['celltype'])
        plt.savefig(output_dir+os.sep+'train_umap_Predcelltype.pdf', bbox_inches='tight')
    

# --- run 10 other experiments:
if args.experiment == 'Chu_CD4_multibatch_validation': # Exp1
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'+os.sep+'CancerType'
    batch_files = glob.glob(Chu_CD4T_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': Chu_CD4T_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': Chu_CD4T_dir+os.sep+'train_80_pancanT_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'Chu_pancancer_CD4T'+os.sep+'test_20_pancanT',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'Chu_CD8_multibatch_validation': # Exp2
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'+os.sep+'CancerType'
    batch_files = glob.glob(Chu_CD8T_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': Chu_CD8T_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': Chu_CD8T_dir+os.sep+'train_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'Chu_pancancer_CD8T'+os.sep+'test_20_pancanT',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'ProjecTILs_CD8_multibatch_to_GSE179994_CD8': # Exp3
    ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'Tissue'
    batch_files = glob.glob(ProjecTIL_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': ProjecTIL_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': ProjecTIL_dir+os.sep+'train_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'GSE179994_CD8T'+os.sep+'GSE179994_CD8T_ref',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'Chu_CD8_multibatch_to_GSE179994_CD8':
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'+os.sep+'CancerType'
    batch_files = glob.glob(Chu_CD8T_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': Chu_CD8T_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': Chu_CD8T_dir+os.sep+'train_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'GSE179994_CD8T'+os.sep+'GSE179994_CD8T_ref',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'ProjecTILs_CD8_multibatch_to_HNSC_CD8':
    ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'Tissue'
    batch_files = glob.glob(ProjecTIL_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': ProjecTIL_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': ProjecTIL_dir+os.sep+'train_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'HNSC_CD8T'+os.sep+'HNSC_CD8T',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'Chu_CD8_multibatch_to_HNSC_CD8':
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'+os.sep+'CancerType'
    batch_files = glob.glob(Chu_CD8T_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': Chu_CD8T_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': Chu_CD8T_dir+os.sep+'train_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'HNSC_CD8T'+os.sep+'HNSC_CD8T',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'ProjecTILs_CD8_multibatch_to_Chu_CD8':
    ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'Tissue'
    batch_files = glob.glob(ProjecTIL_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': ProjecTIL_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': ProjecTIL_dir+os.sep+'train_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'Chu_pancancer_CD8T'+os.sep+'test_20_pancanT',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'Chu_CD8_multibatch_to_ProjecTILs_CD8':
    Chu_CD8T_dir = data_dir+os.sep+'Chu_pancancer_CD8T'+os.sep+'CancerType'
    batch_files = glob.glob(Chu_CD8T_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': Chu_CD8T_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': Chu_CD8T_dir+os.sep+'train_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'ProjecTIL_CD8T',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'Zheng_CD4_multibatch_to_Chu_CD4':
    Zheng_CD4T_dir = data_dir+os.sep+'Zheng_CD4T'+os.sep+'cancerType'
    batch_files = glob.glob(Zheng_CD4T_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': Zheng_CD4T_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': Zheng_CD4T_dir+os.sep+'train_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'Chu_pancancer_CD4T'+os.sep+'test_20_pancanT',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

if args.experiment == 'Chu_CD4_multibatch_to_Zheng_CD4':
    Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'+os.sep+'CancerType'
    batch_files = glob.glob(Chu_CD4T_dir+os.sep+'*.mtx.gz')
    input_batches = []
    for batch_file in batch_files:
        input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
    # --- train configure
    train_input_dict = {'batch_dir': Chu_CD4T_dir,
                        'input': input_batches,
                        'batch_info': list(range(len(input_batches))),
                        'metadata': Chu_CD4T_dir+os.sep+'train_80_pancanT_batch_metadata.csv',
                        'output_dir': output_dir,
                        'prefix': 'train_'}
    train_args = argparse.Namespace(**train_input_dict)
    model = train.train_batch_MLP(train_args, model_config)
    # --- test configure
    test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
            'input': data_dir+os.sep+'Zheng_CD4T'+os.sep+'Zheng_CD4T',
            'output_dir': output_dir,
            'oneround': False, 'prefix': 'scPanKD_predicted_'}
    predict_args = argparse.Namespace(**test_input_dict)
    predict.predict(predict_args, model, model_config)

# if args.experiment == 'GSE179994_CD8T_multibatch_to_HNSC_CD8': # Exp4
#     GSE179994_dir = data_dir+os.sep+'GSE179994_CD8T'+os.sep+'patient'
#     batch_files = glob.glob(GSE179994_dir+os.sep+'*.mtx.gz')
#     input_batches = []
#     for batch_file in batch_files:
#         input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
#     # --- train configure
#     train_input_dict = {'batch_dir': GSE179994_dir,
#                         'input': input_batches,
#                         'batch_info': list(range(len(input_batches))),
#                         'metadata': GSE179994_dir+os.sep+'train_batch_metadata.csv',
#                         'output_dir': output_dir,
#                         'prefix': 'train_'}
#     train_args = argparse.Namespace(**train_input_dict)
#     model = train.train_batch_MLP(train_args, model_config)
#     # --- test configure
#     test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
#             'input': data_dir+os.sep+'HNSC_CD8T'+os.sep+'HNSC_CD8T',
#             'output_dir': output_dir,
#             'oneround': False, 'prefix': 'scPanKD_predicted_'}
#     predict_args = argparse.Namespace(**test_input_dict)
#     predict.predict(predict_args, model, model_config)
# if args.experiment == 'GSE179994_CD8T_multibatch_to_Chu_CD8': # Exp6
#     GSE179994_dir = data_dir+os.sep+'GSE179994_CD8T'+os.sep+'patient'
#     batch_files = glob.glob(GSE179994_dir+os.sep+'*.mtx.gz')
#     input_batches = []
#     for batch_file in batch_files:
#         input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
#     # --- train configure
#     train_input_dict = {'batch_dir': GSE179994_dir,
#                         'input': input_batches,
#                         'batch_info': list(range(len(input_batches))),
#                         'metadata': GSE179994_dir+os.sep+'train_batch_metadata.csv',
#                         'output_dir': output_dir,
#                         'prefix': 'train_'}
#     train_args = argparse.Namespace(**train_input_dict)
#     model = train.train_batch_MLP(train_args, model_config)
#     # --- test configure
#     test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
#             'input': data_dir+os.sep+'Chu_pancancer_CD8T'+os.sep+'test_20_pancanT',
#             'output_dir': output_dir,
#             'oneround': False, 'prefix': 'scPanKD_predicted_'}
#     predict_args = argparse.Namespace(**test_input_dict)
#     predict.predict(predict_args, model, model_config)
# if args.experiment == 'ProjecTILs_CD8_multibatch_to_HNSC_CD8': # Exp5
#     ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'Cohort'
#     batch_files = glob.glob(ProjecTIL_dir+os.sep+'*.mtx.gz')
#     input_batches = []
#     for batch_file in batch_files:
#         input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
#     # --- train configure
#     train_input_dict = {'batch_dir': ProjecTIL_dir,
#                         'input': input_batches,
#                         'batch_info': list(range(len(input_batches))),
#                         'metadata': ProjecTIL_dir+os.sep+'train_batch_metadata.csv',
#                         'output_dir': output_dir,
#                         'prefix': 'train_'}
#     train_args = argparse.Namespace(**train_input_dict)
#     model = train.train_batch_MLP(train_args, model_config)
#     # --- test configure
#     test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
#             'input': data_dir+os.sep+'HNSC_CD8T'+os.sep+'HNSC_CD8T',
#             'output_dir': output_dir,
#             'oneround': False, 'prefix': 'scPanKD_predicted_'}
#     predict_args = argparse.Namespace(**test_input_dict)
#     predict.predict(predict_args, model, model_config)
# if args.experiment == 'ProjecTILs_CD8_multibatch_to_Chu_CD8': # Exp7
#     ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'Cohort'
#     batch_files = glob.glob(ProjecTIL_dir+os.sep+'*.mtx.gz')
#     input_batches = []
#     for batch_file in batch_files:
#         input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
#     # --- train configure
#     train_input_dict = {'batch_dir': ProjecTIL_dir,
#                         'input': input_batches,
#                         'batch_info': list(range(len(input_batches))),
#                         'metadata': ProjecTIL_dir+os.sep+'train_batch_metadata.csv',
#                         'output_dir': output_dir,
#                         'prefix': 'train_'}
#     train_args = argparse.Namespace(**train_input_dict)
#     model = train.train_batch_MLP(train_args, model_config)
#     # --- test configure
#     test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
#             'input': data_dir+os.sep+'Chu_pancancer_CD8T'+os.sep+'test_20_pancanT',
#             'output_dir': output_dir,
#             'oneround': False, 'prefix': 'scPanKD_predicted_'}
#     predict_args = argparse.Namespace(**test_input_dict)
#     predict.predict(predict_args, model, model_config)

# if args.experiment == 'Chu_GSE179994_CD8T_multibatch_to_HNSC_CD8':
#     pass
