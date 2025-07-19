import os, glob
import argparse
import torch

# import my packages
import train
import predict
from utils import _utils

_utils._set_seed(1993)

if __name__ == '__main__':
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
    model_config = {'fs': 'F-test', 'num_features': 3000, 'optimizer': 'adamW',
                    'learning_rate': 0.001, 'batch_size': 32, 'max_epochs': 100, "distillation_epochs": 30,
                    'MLP_DIMS': [256, 64, 16], 'teacher_MLP_DIMS': [64, 16], 'student_MLP_DIMS': [64, 16],
                    'Celltype_COLUMN': "celltype", 'PredCelltype_COLUMN': "pred_celltype",
                    "entropy_quantile": 0.4}
    model_config['device'] = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
    if args.experiment == 'Chu_CD4T_multibatch_validation': # Exp1
        Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'+os.sep+'CancerType'
        #Chu_CD4T_dir = data_dir+os.sep+'Chu_pancancer_CD4T'+os.sep+'batch'
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
    if args.experiment == 'Chu_CD8T_multibatch_validation': # Exp2
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
    if args.experiment == 'Zheng_CD4_multibatch_to_Chu_CD4': # Exp3
        Zheng_dir = data_dir+os.sep+'Zheng_CD4T'+os.sep+'patient'
        batch_files = glob.glob(Zheng_dir+os.sep+'*.mtx.gz')
        input_batches = []
        for batch_file in batch_files:
            input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
        # --- train configure
        train_input_dict = {'batch_dir': Zheng_dir,
                            'input': input_batches,
                            'batch_info': list(range(len(input_batches))),
                            'metadata': Zheng_dir+os.sep+'train_batch_metadata.csv',
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
    if args.experiment == 'GSE179994_CD8T_multibatch_to_HNSC_CD8': # Exp4
        GSE179994_dir = data_dir+os.sep+'GSE179994_CD8T'+os.sep+'patient'
        batch_files = glob.glob(GSE179994_dir+os.sep+'*.mtx.gz')
        input_batches = []
        for batch_file in batch_files:
            input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
        # --- train configure
        train_input_dict = {'batch_dir': GSE179994_dir,
                            'input': input_batches,
                            'batch_info': list(range(len(input_batches))),
                            'metadata': GSE179994_dir+os.sep+'train_batch_metadata.csv',
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
    if args.experiment == 'GSE179994_CD8T_multibatch_to_Chu_CD8': # Exp6
        GSE179994_dir = data_dir+os.sep+'GSE179994_CD8T'+os.sep+'patient'
        batch_files = glob.glob(GSE179994_dir+os.sep+'*.mtx.gz')
        input_batches = []
        for batch_file in batch_files:
            input_batches.append(os.path.basename(batch_file).replace('.mtx.gz', ''))
        # --- train configure
        train_input_dict = {'batch_dir': GSE179994_dir,
                            'input': input_batches,
                            'batch_info': list(range(len(input_batches))),
                            'metadata': GSE179994_dir+os.sep+'train_batch_metadata.csv',
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
    if args.experiment == 'ProjecTILs_CD8_multibatch_to_HNSC_CD8': # Exp5
        ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'Cohort'
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
    if args.experiment == 'ProjecTILs_CD8_multibatch_to_Chu_CD8': # Exp7
        ProjecTIL_dir = data_dir+os.sep+'ProjecTILs_CD8T'+os.sep+'Cohort'
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

    if args.experiment == 'Chu_GSE179994_CD8T_multibatch_to_HNSC_CD8':
        pass


