import os

import fire
import celltypist
from celltypist import models


def do_predict(
    mtx_path: str,
    gene_path: str,
    cell_path: str,
    out_dir: str,
    is_gene_by_cell: bool = True,
    model_name: str = 'Immune_All_Low.pkl',
):
    """
    # mtx_file: need to take a look inside to know gene by cell or cell by gene
    # gene_file: genes.tsv
    # cell_file: barcodes.tsv

    In [2]: mtx = scipy.io.mmread('datasets/Chu_CD4T_validation/test_20_pancanT.mtx.gz')
    Out[3]: (12191, 18914)  -> gene by cell matrix
    """
    model = models.Model.load(model=model_name)
    predictions = celltypist.annotate(
        mtx_path, 
        model=model, 
        transpose_input=is_gene_by_cell, 
        gene_file=gene_path, 
        cell_file=cell_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    predictions.to_table(folder=out_dir, prefix='')


if __name__ == '__main__':
    fire.Fire(do_predict)