import os

import fire
import pandas as pd
import celltypist
from celltypist import models


def do_train(
    mtx_path: str,
    label_path: str,
    gene_path: str,
    out_dir: str,
    is_gene_by_cell: bool = True,
):
    label_df = pd.read_csv(label_path)
    mini_batch_flag = True if len(label_df) >= 1000 else False
    new_model = celltypist.train(
        mtx_path, 
        labels=label_path, 
        genes=gene_path, 
        transpose_input=is_gene_by_cell,
        n_jobs=4,
        use_SGD=True,
        mini_batch=mini_batch_flag,
        )
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_path = os.path.join(out_dir, 'default_model.pkl')
    print(f'Saving model to {out_path}')
    new_model.write(out_path)


def extract_label(
    metadata_path: str,
    col: str = 'celltype',
):
    out_path = metadata_path.replace('.csv', '_label.csv')
    print(f'Extracting label from {metadata_path} to {out_path}')
    df = pd.read_csv(metadata_path)
    print(f'{df.shape=}')
    out_df = df[[col]]
    out_df.to_csv(out_path, sep='\t', index=False, header=False)


if __name__ == '__main__':
    fire.Fire()