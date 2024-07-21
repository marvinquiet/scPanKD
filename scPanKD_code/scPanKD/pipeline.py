import os, glob
import argparse
from train import *
from predict import *

if __name__ == '__main__':
    project_dir = '/panfs/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang'
    data_dir = project_dir+os.sep+'data'
    result_dir = project_dir+os.sep+'results/scPanKD'
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment', type=str, help='Experiments to run')
    args = parser.parse_args()

    if args.experiment == 'Chu_CD4T_multibatch_validation':
        output_dir = result_dir+os.sep+args.experiment
        os.makedirs(output_dir, exist_ok=True)
        Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'+os.sep+'CancerType'
        batch_files = glob.glob(Chu_CD4T_dir+os.sep+'*.mtx.gz')
        input_batches = []
        for batch_file in batch_files:
            input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
        # --- train configure
        train_input_dict = {'cmd_choice': 'train_batch',
            'batch_dir': Chu_CD4T_dir,
            'input': input_batches,
            'batch_info': list(range(len(input_batches))),
            'metadata': Chu_CD4T_dir+os.sep+'train_80_pancanT_batch_metadata.csv',
            'model': 'MLP_batch',
            'output_dir': output_dir,
            'prefix': 'train_', 'fs': 'F-test', 'num_features': 3000}
        train_parser = argparse.ArgumentParser(description='Process existing arguments')
        for key, value in train_input_dict.items():
            train_parser.add_argument(f'--{key}', type=type(value), default=value, help=f'{key} description')
        train_args = train_parser.parse_args()
        if train_args.cmd_choice == 'train_batch':
            train_bath_MLP(train_args)
        if train_args.cmd_choice == 'train':
            train_MLP(train_args)
        # --- test configure
        test_input_dict = {'trained_model': output_dir+os.sep+'train_MLP_model',
             'input': data_dir+os.sep+'Chu_pancancer_CD4T'+os.sep+'test_20_pancanT',
             'output_dir': output_dir,
             'oneround': False, 'prefix': 'Cellcano_predicted_'}
        predict_parser = argparse.ArgumentParser(description='Process existing arguments')
        for key, value in test_input_dict.items():
            predict_parser.add_argument(f'--{key}', type=type(value), default=value, help=f'{key} description')
        predict_args = predict_parser.parse_args()
        predict(predict_args)
    if args.experiment == 'Chu_CD8T_multibatch_validation':
        pass

    if args.experiment == 'Chu_GSE179994_CD8T_multibatch_to_HNSC_CD8':
        pass

    if args.experiment == 'ProjecTILs_CD8_multibatch_to_HNSC_CD8':
        pass

    if args.experiment == 'ProjecTILs_CD8_multibatch_to_Chu_CD8':
        pass


