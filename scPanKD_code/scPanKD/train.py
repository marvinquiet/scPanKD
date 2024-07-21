'''
Functions related to train models
'''
import os, sys

import anndata
import numpy as np
from sklearn.preprocessing import OneHotEncoder

## import my package
#from Cellcano.utils import _utils
from utils import _utils

## get the logger
#import logging
#logger = logging.getLogger(__name__)

def load_train_adata(args):
    ''' Load training data
    '''
    if args.anndata is not None:
        print("Load pre-processed anndata h5ad file..")
        train_adata = anndata.read_h5ad(args.anndata)
    else:
        if args.input is None or args.metadata is None:
            sys.exit("Please make sure that both gene score matrix and metadata are provided!")

        ## load input data
        print("Loading data... \n This may take a while depending on your data size..")
        if '.csv' in args.input:
            train_adata = _utils._csv_data_loader(args.input)
        else:
            train_adata = _utils._COOmtx_data_loader(args.input)
        metadata = _utils._metadata_loader(args.metadata)
        if _utils.Celltype_COLUMN not in metadata.columns:
            sys.exit("Column '%s' is not found in metadata. Please make sure to include cell type information in the metadata." % _utils.Celltype_COLUMN)

        common_cells = set(train_adata.obs_names).intersection(set(metadata.index))
        print("%d common cells found between input data and metadata." % len(common_cells))

        if len(common_cells) == 0:
            sys.exit("No common cells are found between input data and metadata, please check your data!")

        if len(common_cells) < 100:
            print("There are too few cells. Cellcano might not be accurate.")

        train_adata = train_adata[list(common_cells)]
        train_adata.obs = train_adata.obs.merge(metadata,
                left_on="barcode", right_index=True, how='left')

        ## preprocess data and select features
        train_adata = _utils._process_adata(train_adata, process_type='train')
        print("Data shape after processing: %d cells X %d genes" % (train_adata.shape[0], train_adata.shape[1]))
        train_adata = _utils._select_feature(train_adata,
                fs_method=args.fs, num_features=args.num_features)
        train_adata = _utils._scale_data(train_adata) ## center-scale
        _utils._visualize_data(train_adata, args.output_dir, prefix=args.prefix)
        _utils._save_adata(train_adata, args.output_dir, prefix=args.prefix)
    return train_adata

def load_train_batch_adata(args):
    ''' Load batch of train anndata
    '''
    if args.input is None or args.metadata is None:
        sys.exit("Please make sure that both gene score matrix and metadata are provided!")

    ## load list of input
    adata_list = []
    for batch in args.input:
        adata = _utils._COOmtx_data_loader(args.batch_dir+os.sep+batch)
        adata_list.append(adata)
    assert len(adata_list) == len(args.batch_info), "Input data must have the same length as batch labels"
    # concate adata
    batch_labels = []
    for i in range(len(adata_list)):
        num_cells = adata_list[i].shape[0]
        batch_labels.extend([args.batch_info[i]]*num_cells)
    train_adata = anndata.concat(adata_list, axis=0)
    train_adata.obs['batch'] = batch_labels
    metadata = _utils._metadata_loader(args.metadata)
    if _utils.Celltype_COLUMN not in metadata.columns:
        sys.exit("Column '%s' is not found in metadata. Please make sure to include cell type information in the metadata." % _utils.Celltype_COLUMN)
    # check common cells between input and metadata
    common_cells = set(train_adata.obs_names).intersection(set(metadata.index))
    print("%d common cells found between input data and metadata." % len(common_cells))
    if len(common_cells) == 0:
        sys.exit("No common cells are found between input data and metadata, please check your data!")
    if len(common_cells) < 100:
        print("There are too few cells. Cellcano might not be accurate.")
    train_adata = train_adata[list(common_cells)]
    train_adata.obs = train_adata.obs.merge(metadata,
            left_on="barcode", right_index=True, how='left', suffixes=('', '_x')) ## keep left batch information
    ## preprocess data and select features
    train_adata = _utils._process_adata(train_adata, process_type='train')
    print("Data shape after processing: %d cells X %d genes" % (train_adata.shape[0], train_adata.shape[1]))
    train_adata = _utils._select_feature(train_adata,
            fs_method=args.fs, num_features=args.num_features)
    train_adata = _utils._scale_data(train_adata) ## center-scale
    _utils._visualize_data(train_adata, args.output_dir, prefix=args.prefix)
    _utils._save_adata(train_adata, args.output_dir, prefix=args.prefix)
    return train_adata

def train_MLP(args):
    '''Train MLP model
        1. Load train data and metadata
        2. Feature selection
        3. Log-Norm and scale data
        4. Train model and save
    ---
    Input:
        - args: user's input arguments
        - train_adata: train anndata object

    '''
    MLP_DIMS = _utils.MLP_DIMS #if args.mlp_ns is None else args.mlp_ns

    train_adata = load_train_adata(args)
    ## get x_train, y_train
    x_train = _utils._extract_adata(train_adata)
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(train_adata.obs[[_utils.Celltype_COLUMN]]).toarray()
    print("Cell type categories: ", enc.categories_[0])

    mlp = _utils._init_MLP(x_train, y_train, dims=MLP_DIMS,
            seed=_utils.RANDOM_SEED)
    mlp.compile()
    mlp.fit(x_train, y_train)
    model_save_dir = args.output_dir+os.sep+args.prefix+'MLP_model'
    mlp.model.save(model_save_dir)

    ## save feature information along with mean and standard deviation
    train_adata.var.loc[:, ['mean', 'std']].to_csv(model_save_dir+os.sep+"features.txt", sep='\t')
    ## save enc information
    with open(model_save_dir+os.sep+"onehot_encoder.txt", 'w') as f:
        for idx, cat in enumerate(enc.categories_[0]):
            f.write('%d:%s\n' % (idx, cat))

def train_bath_MLP(args):
    '''Train MLP model based on concatenated batch vector
        1. Load train data and metadata from each batch
        2. Concat the anndata, perform the feature selectio, log-norm, scale
        3. Train MLP+batch model and save
    ---
    Input:
        - args: user's input argument
    '''
    train_adata = load_train_batch_adata(args)
    batch_info = np.array(train_adata.obs['batch']).astype('int')
    batch_mat = np.zeros((batch_info.size, batch_info.max()+1))
    batch_mat[np.arange(batch_info.size), batch_info] = 1
    n_batch = batch_mat.shape[1]
    # get MLP layers
    MLP_DIMS = _utils.MLP_DIMS #if args.mlp_ns is None else args.mlp_ns
    ## get x_train, y_train
    x_train = _utils._extract_adata(train_adata)
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(train_adata.obs[[_utils.Celltype_COLUMN]]).toarray()
    print("Cell type categories: ", enc.categories_[0])
    # initialize batch MLP model
    mlp = _utils._init_batch_MLP(x_train, batch_mat, y_train, dims=MLP_DIMS,
            seed=_utils.RANDOM_SEED)
    # fix node
    mlp.compile()
    mlp.fit(x_train, batch_mat, y_train)
    model_save_dir = args.output_dir+os.sep+args.prefix+'MLP_model'
    # save backbone model
    mlp_backbone = _utils._extract_backbone_mlp(mlp.model)
    mlp_backbone.save(model_save_dir)

    ## save feature information along with mean and standard deviation
    train_adata.var.loc[:, ['mean', 'std']].to_csv(model_save_dir+os.sep+"features.txt", sep='\t')
    ## save enc information
    with open(model_save_dir+os.sep+"onehot_encoder.txt", 'w') as f:
        for idx, cat in enumerate(enc.categories_[0]):
            f.write('%d:%s\n' % (idx, cat))

def train_KD(args):
    '''Train one step KD model
        1. Load train data and metadata
        2. Feature selection
        3. Log-Norm and scale data
        4. Train model and save
    ---
    Input:
        - args: user's input arguments
        - train_adata: train anndata object
    '''
    teacher_MLP_DIMS = _utils.Teacher_DIMS #if args.teacher_ns is None else args.teacher_ns
    student_MLP_DIMS = _utils.Student_DIMS #if args.student_ns is None else args.student_ns

    train_adata = load_train_adata(args)
    ## get x_train, y_train
    x_train = _utils._extract_adata(train_adata)
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(train_adata.obs[[_utils.Celltype_COLUMN]]).toarray()
    print("Cell type categories: ", enc.categories_)

    ## train a KD model
    teacher = _utils._init_MLP(x_train, y_train, dims=teacher_MLP_DIMS,
            seed=_utils.RANDOM_SEED)
    teacher.compile()
    teacher.fit(x_train, y_train, batch_size=_utils.BATCH_SIZE)
    student = _utils._init_MLP(x_train, y_train, dims=student_MLP_DIMS,
            seed=_utils.RANDOM_SEED)
    distiller = _utils._run_distiller(x_train, y_train,
            student_model=student.model,
            teacher_model=teacher.model)
    distiller.student.save(args.output_dir+os.sep+args.prefix+'KD_model')

#if __name__ == '__main__':
#    # --- train procedure
#    #input_dict = {'cmd_choice': 'train_batch',
#    #              'batch_dir': '/net/mulan/home/wenjinma/minicluster/projects/test_Cellcano_XuefengWang/CD8T_ref_data',
#    #              'input': ['HT2.2', 'HT3.3', 'SCT1.2', 'SCT1.3', 'RT14.1',
#    #                        'LT21.1', 'LT21.2', 'LT21.4',
#    #                        'LT26.2', 'LT26.3', 'LT35.2', 'LT39', 'LT40', 'LT58.2'],
#    #              'batch_info': [1, 1, 2, 2, 3,
#    #                             4, 4, 4,
#    #                             5, 5, 5, 5, 5, 5],
#    #              'metadata': '/net/mulan/home/wenjinma/minicluster/projects/test_Cellcano_XuefengWang/CD8T_ref_metadata.csv',
#    #              'model': 'MLP_batch',
#    #              'output_dir': '/net/mulan/home/wenjinma/minicluster/projects/test_Cellcano_XuefengWang/CD8T_ref_batch_result',
#    #              'prefix': 'train_',
#    #              'fs': 'F-test', 'num_features': 3000}
#
#    input_dict = {'cmd_choice': 'train_batch',
#                  'batch_dir': '/projects/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/data/ProjecTILs_CD8T',
#                  'input': ['ProjecTIL_GSE180268_CD8T_trimmed',
#                            'ProjecTIL_GSE159251_CD8T_trimmed',
#                            'ProjecTIL_GSE123814_CD8T_trimmed',
#                            'ProjecTIL_GSE176021_CD8T_trimmed',
#                            'ProjecTIL_GSE139555_CD8T_trimmed',
#                            'ProjecTIL_PRJNA705464_CD8T_trimmed',
#                            'ProjecTIL_GSE179994_CD8T_trimmed',
#                            'ProjecTIL_EGAS00001004809_CD8T_trimmed'],
#                  'batch_info': [1, 2, 3, 4, 5, 6, 7, 8],
#                  'metadata': '/projects/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/data/ProjecTILs_CD8T/ProjecTIL_CD8T_metadata.csv',
#                  'model': 'MLP_batch',
#                  'output_dir': '/projects/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/results/Cellcano_batch/ProjecTIL_CD8T_ref_batch_results',
#                  'prefix': 'train_',
#                  'fs': 'F-test', 'num_features': 3000}
#    import argparse
#    parser = argparse.ArgumentParser(description='Process existing arguments')
#    for key, value in input_dict.items():
#        parser.add_argument(f'--{key}', type=type(value), default=value, help=f'{key} description')
#    args = parser.parse_args()
#    if args.cmd_choice == 'train_batch':
#        train_bath_MLP(args)
#
#    # --- predict procedure
#    import predict
#    input_dict = {'trained_model': '/projects/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/results/Cellcano_batch/ProjecTIL_CD8T_ref_batch_results/train_MLP_model',
#            'input': '/projects/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/data/HNSC_CD8T/HNSC_CD8T',
#            'oneround': False,
#            'output_dir': '/projects/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/results/Cellcano_batch/ProjecTIL_CD8T_ref_batch_results',
#            'prefix': 'HNSC_CD8T'}
#    predict_parser = argparse.ArgumentParser(description='Process existing arguments')
#    for key, value in input_dict.items():
#        predict_parser.add_argument(f'--{key}', type=type(value), default=value, help=f'{key} description')
#    predict_args = predict_parser.parse_args()
#    predict.predict(predict_args)
#
