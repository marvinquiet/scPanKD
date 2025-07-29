'''
Functions related to train models
'''
import os, sys
import anndata
import numpy as np
from sklearn.preprocessing import OneHotEncoder

import torch
from tqdm.auto import tqdm
## import my package
from models.cancer_MLP import cancer_MLP
from utils import _utils

def load_train_batch_adata(args, model_config):
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
    if model_config['Celltype_COLUMN'] not in metadata.columns:
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
    train_adata = _utils._select_feature(train_adata, model_config)
    # TODO: test whether center scale for each batch would help.. -> not sure yet
    #train_adata = _utils._scale_data(train_adata) ## center-scale
    batch_adata_list = []
    for batch in set(train_adata.obs['batch']):
        batch_adata = train_adata[train_adata.obs['batch'] == batch]
        batch_adata = _utils._scale_data(batch_adata)
        batch_adata_list.append(batch_adata)
    train_adata = anndata.concat(batch_adata_list, axis=0)
    _utils._visualize_data(train_adata, args.output_dir, prefix=args.prefix, reduction="UMAP")
    _utils._save_adata(train_adata, args.output_dir, prefix=args.prefix)
    return train_adata

def train_batch_MLP(args, model_config):
    '''Train MLP model based on concatenated batch vector
        1. Load train data and metadata from each batch
        2. Concat the anndata, perform the feature selectio, log-norm, scale
        3. Train MLP+batch model and save
    ---
    Input:
        - args: user's input argument
    '''
    train_adata = load_train_batch_adata(args, model_config)
    cancer_info = np.array(train_adata.obs['batch']).astype('int')
    ## get x_train, y_train
    x_train = _utils._extract_adata(train_adata)
    # generate a list of gene expression matrix for each cancer type with x_train
    reference_gene_expr_by_cancer = []
    for cancer in set(cancer_info):
        cancer_idx = np.where(cancer_info == cancer)[0]
        reference_gene_expr_by_cancer.append(torch.tensor(x_train[cancer_idx, :]).float().to(model_config['device']))
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(train_adata.obs[[model_config['Celltype_COLUMN']]]).toarray()
    print("Cell type categories: ", enc.categories_[0])
    # initialize cancer MLP model
    cancer_mlp = cancer_MLP(dims=model_config['MLP_DIMS'], 
            input_dim=x_train.shape[1],
            n_batch=len(set(cancer_info)),
            n_classes=y_train.shape[1]).to(model_config['device'])
    if model_config['optimizer'] == 'adam':
        optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, cancer_mlp.parameters()), lr=model_config['learning_rate'])
    if model_config['optimizer'] == 'adamW':
        optimizer = torch.optim.AdamW(filter(lambda p: p.requires_grad, cancer_mlp.parameters()), lr=model_config['learning_rate'])
    criterion = torch.nn.CrossEntropyLoss()
    # triplet_loss_fn = torch.nn.TripletMarginLoss(margin=1.0, p=2)

    dataset = torch.utils.data.TensorDataset(torch.tensor(x_train, dtype=torch.float32, device=model_config['device']),
                                             torch.tensor(cancer_info, dtype=torch.long, device=model_config['device']),
                                             torch.tensor(y_train, dtype=torch.float32, device=model_config['device'])) 
    dataloader = torch.utils.data.DataLoader(dataset, batch_size=model_config['batch_size'], shuffle=True)
    cancer_mlp.train()
    for epoch in tqdm(range(model_config['max_epochs'])):
        total_loss, total_loss_cls, total_loss_adv = 0, 0, 0
        correct, total = 0, 0
        for batch_x, batch_c, batch_y in dataloader:
            optimizer.zero_grad()
            outputs, _ = cancer_mlp(batch_x)
            loss_cls = criterion(outputs, torch.argmax(batch_y, dim=1))
            # add batch adversarial loss
            cancer_logits = cancer_mlp(batch_x, reverse=True)
            loss_adv = criterion(cancer_logits, batch_c)
            loss = loss_cls + 10 * loss_adv 
            loss.backward()
            optimizer.step()
            # record training loss and accuracy
            total_loss += loss.item()
            total_loss_cls += loss_cls.item()
            total_loss_adv += loss_adv.item()
            # calculate accuracy
            _, predicted = torch.max(outputs.data, 1)
            correct += (predicted == torch.argmax(batch_y, dim=1)).sum().item()
            total += batch_y.size(0)
        if (epoch + 1) % 10 == 0:
            avg_loss = total_loss / len(dataloader)
            avg_cls_loss = total_loss_cls / len(dataloader)
            avg_adv_loss = total_loss_adv / len(dataloader)
            accuracy = 100 * correct / total
            print(f"Epoch [{epoch + 1}/{model_config['max_epochs']}], Total Loss: {avg_loss:.4f}, Classification Loss: {avg_cls_loss:.4f}, Adversatial Loss: {avg_adv_loss:.4f}, Accuracy: {accuracy:.2f}%")

    model_save_dir = args.output_dir+os.sep+args.prefix+'MLP_model'
    if not os.path.exists(model_save_dir):
        os.makedirs(model_save_dir)
    torch.save(cancer_mlp.state_dict(), model_save_dir+os.sep+'model.pth')
    
    
    # visualize predicted features for training data
    _, features, _ = cancer_mlp.predict(torch.tensor(x_train, dtype=torch.float32, device=model_config['device']),
                                        reference_gene_expr_by_cancer,
                                        vote_strategy='majority')
    feature_adata = anndata.AnnData(X=features.cpu().detach().numpy(), obs=train_adata.obs)
    _utils._visualize_embedding(feature_adata, args.output_dir, color_columns=['celltype', 'batch'], prefix=args.prefix+'embedding_', reduction="UMAP")

    ## save feature information along with mean and standard deviation
    train_adata.var.to_csv(model_save_dir+os.sep+"features.txt", sep='\t')
    ## save enc information
    with open(model_save_dir+os.sep+"onehot_encoder.txt", 'w') as f:
        for idx, cat in enumerate(enc.categories_[0]):
            f.write('%d:%s\n' % (idx, cat))
    return cancer_mlp