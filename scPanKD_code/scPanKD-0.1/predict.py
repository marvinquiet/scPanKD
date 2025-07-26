'''
Functions related to predict cell types
'''
import os, sys 

import torch
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

## import my package
from utils import _utils
from models.MLP import MLP

import anndata

def predict(args, model, model_config):
    feature_file = args.trained_model+os.sep+'features.txt'
    encoder_file = args.trained_model+os.sep+'onehot_encoder.txt'
    if not os.path.exists(feature_file) or not os.path.exists(encoder_file):
        sys.exit("Feature file or encoder mapping does not exist! Please check your trained model was trained successfully.")

    features = pd.read_csv(feature_file, sep='\t', header=0, index_col=0)
    encoders = {}
    with open(encoder_file) as f:
        for line in f:
            line_info = line.strip().split(':')
            encoders[int(line_info[0])] = line_info[1]

    ## load input data
    print("Loading data... \n This may take a while depending on your data size..")
    if '.csv' in args.input:
        test_adata = _utils._csv_data_loader(args.input)
    else:
        test_adata = _utils._COOmtx_data_loader(args.input)
    ## process test adata
    test_adata = _utils._process_adata(test_adata, process_type='test')

    ## fill in the data with the same order of features
    feature_idx = []
    NA_idx = []
    for f_idx, feature in enumerate(features.index):
        find_flag = False
        for test_idx, gene in enumerate(test_adata.var_names):
            if gene == feature:
                feature_idx.append(test_idx)
                find_flag = True
                break
        if not find_flag:
            feature_idx.append(-1)
            NA_idx.append(f_idx)
    print("%d genes from reference data are found in target.\n" % (len(features)-len(NA_idx)))

    if len(NA_idx) > 0.1 * len(features):
        print("Warnings: too few genes found in target and this will result in inaccurate prediction.")
    if -1 in feature_idx:
        print("Warnings: since some feature does not exist in target dataset. We will fill in 0s for those columns.")
        ## first replace those unique genes with index 
        curated_feature_idx = np.array(feature_idx)
        curated_feature_idx[NA_idx] = 0
        test_adata = test_adata[:, curated_feature_idx].copy()
        test_adata.var_names.values[NA_idx] = ["GenesNotFound-"+str(i) for i, NA_item in enumerate(NA_idx)]  ## change gene names
        test_adata_X = test_adata.X
        test_adata_X[:, NA_idx] = 0
        test_adata.X = test_adata_X
    else:
        test_adata = test_adata[:, feature_idx]
    print("Data shape after processing: %d cells X %d genes"  % (test_adata.shape[0], test_adata.shape[1]))

    test_adata = _utils._scale_data(test_adata)
    test_data_mat = _utils._extract_adata(test_adata)
    test_data_tensor = torch.tensor(test_data_mat, dtype=torch.float32, device=model_config['device'])

    y_pred, y_pred_features = model.predict(test_data_tensor)
    y_pred = y_pred.detach().cpu().numpy()
    model.eval()
    y_pred_batch = model.batch_predictor(y_pred_features).detach().cpu().numpy()
    y_pred_batch = np.argmax(y_pred_batch, axis=1)
    test_adata.obs['pred_batch'] = y_pred_batch
    pred_celltypes = _utils._prob_to_label(y_pred, encoders)
    test_adata.obs[model_config['PredCelltype_COLUMN']] = pred_celltypes

    # --- visualize predicted features\
    metadata = pd.read_csv(args.input+'_metadata.csv', header=0, index_col=0)
    if 'curated_celltype' in metadata.columns:
        test_adata.obs['true_label'] = metadata['curated_celltype'].values
    else:
        test_adata.obs['true_label'] = metadata['celltype'].values

    y_pred_features = y_pred_features.detach().cpu().numpy()
    y_pred_features_adata = anndata.AnnData(X=y_pred_features, obs=test_adata.obs)
    _utils._visualize_embedding(y_pred_features_adata, args.output_dir, color_columns=[model_config['PredCelltype_COLUMN'], 'true_label'], prefix=args.prefix+'embedding_pred_', reduction="UMAP")
    
    ## if less than 1000 in test data
    if test_adata.shape[0] < 1000:
        print("Your input cell is less than 1000. For performance, we will not perform two-round strategy on your data.")
        test_adata.obs[['pred_celltype']].to_csv(args.output_dir+os.sep+args.prefix+'celltypes.csv')
    ## when cell number is large enough
    else:
        firstround_COLUMN = 'firstround_' + model_config['PredCelltype_COLUMN']
        test_adata.obs[firstround_COLUMN] = pred_celltypes
        entropy = [-np.nansum(y_pred[i]*np.log(y_pred[i])) for i in range(y_pred.shape[0])]
        test_adata.obs['entropy'] = entropy
        test_adata = _utils._select_confident_cells(
                test_adata, celltype_col=firstround_COLUMN, entropy=model_config['entropy_quantile'])

        low_entropy_cells = test_adata.obs_names[np.where(test_adata.obs['entropy_status'] == 'low')].tolist()
        high_entropy_cells = test_adata.obs_names[np.where(test_adata.obs['entropy_status'] == 'high')].tolist()
        test_ref_adata = test_adata[low_entropy_cells]
        test_tgt_adata = test_adata[high_entropy_cells]

        ## enlarge reference dataset for second round
        sampled_ref_adata = _utils._oversample_cells(test_ref_adata, 
                celltype_col=firstround_COLUMN)
        x_tgt_train = _utils._extract_adata(sampled_ref_adata)
        x_tgt_train_batch = sampled_ref_adata.obs['pred_batch'].values
        y_tgt_train = _utils._label_to_onehot(sampled_ref_adata.obs[firstround_COLUMN].tolist(),
                encoders=encoders)
        
        # add training data information
        tgt_dataset = torch.utils.data.TensorDataset(torch.tensor(x_tgt_train, dtype=torch.float32, device=model_config['device']),
                                                     torch.tensor(x_tgt_train_batch, dtype=torch.long, device=model_config['device']),
                                                     torch.tensor(y_tgt_train, dtype=torch.float32, device=model_config['device'])) 
        tgt_dataloader = torch.utils.data.DataLoader(tgt_dataset, batch_size=model_config['batch_size'], shuffle=True)

        # --- retrain model as teacher model
        if model_config['optimizer'] == 'adam':
            teacher_optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=3e-4)
        if model_config['optimizer'] == 'adamW':
            teacher_optimizer = torch.optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=3e-4)
        criterion = torch.nn.CrossEntropyLoss()
        model.train()
        for _ in tqdm(range(model_config['max_epochs'])):
            for batch_x, batch_b, batch_y in tgt_dataloader:
                teacher_optimizer.zero_grad()
                outputs, _ = model(batch_x, batch_b)
                loss_cls = criterion(outputs, torch.argmax(batch_y, dim=1))
                # add batch adversarial loss
                batch_logits = model(batch_x, batch_b, reverse=True)
                loss_adv = criterion(batch_logits, batch_b)
                loss = loss_cls + 10 * loss_adv 
                loss.backward()
                teacher_optimizer.step()

        student = MLP(dims=model_config['student_MLP_DIMS'], input_dim=x_tgt_train.shape[1], n_classes=y_tgt_train.shape[1])
        student = student.to(model_config['device'])
        student_optimizer = torch.optim.AdamW(student.parameters(), lr=3e-4)
        model.eval()
        student.train()
        for epoch in tqdm(range(model_config['distillation_epochs'])):
            total_loss = 0
            for batch_x, batch_b, batch_y in tgt_dataloader:
                # Forward pass
                teacher_outputs, _ = model(batch_x, batch_b)
                teacher_outputs = teacher_outputs.detach()
                student_outputs = student(batch_x)
                loss = _utils._distillation_loss(student_outputs, teacher_outputs, batch_y, alpha=0.1, T=3)
                # Backward pass and optimization
                student_optimizer.zero_grad()
                loss.backward()
                student_optimizer.step()
                total_loss += loss.item()
            # Print loss every 10 epochs
            if (epoch + 1) % 10 == 0:
                avg_loss = total_loss / len(tgt_dataloader)
                print(f"Epoch [{epoch + 1}/{model_config['distillation_epochs']}], Loss: {avg_loss:.4f}")

        x_tgt_test = _utils._extract_adata(test_tgt_adata)
        y_pred_tgt = student.predict(torch.tensor(x_tgt_test, dtype=torch.float32, device=model_config['device']))
        y_pred_tgt = y_pred_tgt.detach().cpu().numpy()

        pred_celltypes = _utils._prob_to_label(y_pred_tgt, encoders)
        test_adata.obs.loc[high_entropy_cells, model_config['PredCelltype_COLUMN']] = pred_celltypes
        ## select certain columns and store to the file
        test_adata.obs[['pred_celltype', 'firstround_pred_celltype', 'entropy']].to_csv(args.output_dir+os.sep+args.prefix+'celltypes.csv')


        # # teacher/student model on original celltype label
        # teacher = MLP(dims=model_config['teacher_MLP_DIMS'], input_dim=x_tgt_train.shape[1], n_classes=y_tgt_train.shape[1])
        # teacher = teacher.to(model_config['device'])
        # criterion = torch.nn.CrossEntropyLoss()
        # teacher_optimizer = torch.optim.AdamW(teacher.parameters(), lr=model_config['learning_rate'])
        # teacher.train()
        # for _ in tqdm(range(model_config['max_epochs'])):
        #     for batch_x, batch_y in tgt_dataloader:
        #         teacher_optimizer.zero_grad()
        #         outputs = teacher(batch_x)
        #         loss = criterion(outputs, batch_y)
        #         loss.backward()
        #         teacher_optimizer.step()
        # student = MLP(dims=model_config['student_MLP_DIMS'], input_dim=x_tgt_train.shape[1], n_classes=y_tgt_train.shape[1])
        # student = student.to(model_config['device'])
        # # Initialize and compile distiller
        # trained_student = _utils._run_distiller(
        #     tgt_dataloader,
        #     student_model=student,
        #     teacher_model=teacher,
        #     epochs=model_config['distillation_epochs'])

        # x_tgt_test = _utils._extract_adata(test_tgt_adata)
        # y_pred_tgt = trained_student.predict(torch.tensor(x_tgt_test, dtype=torch.float32, device=model_config['device']))
        # y_pred_tgt = y_pred_tgt.detach().cpu().numpy()

        # pred_celltypes = _utils._prob_to_label(y_pred_tgt, encoders)
        # test_adata.obs.loc[high_entropy_cells, model_config['PredCelltype_COLUMN']] = pred_celltypes
        # ## select certain columns and store to the file
        # test_adata.obs[['pred_celltype', 'firstround_pred_celltype', 'entropy']].to_csv(args.output_dir+os.sep+args.prefix+'celltypes.csv')