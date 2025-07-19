import anndata
import pandas as pd

# --- load data ---
def load_mtx(mtx_prefix):
    print(f"Loading data from {mtx_prefix}")
    adata = anndata.read_mtx(mtx_prefix+'.mtx.gz').T
    genes = pd.read_csv(mtx_prefix+'_genes.tsv', header=None, sep='\t')
    adata.var["genes"] = genes[0].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None
    cells = pd.read_csv(mtx_prefix+'_barcodes.tsv', header=None, sep='\t')
    adata.obs["barcode"] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None
    adata = adata[:, adata.var_names.notnull()]
    adata.var_names=[i.upper() for i in list(adata.var_names)]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None
    return adata

