# convert genes from Gene Symbols to Ensembl IDs
.libPaths("/net/mulan/home/wenjinma/Rlib")
# --- obtain python through command line
python_path = system("which python", intern = TRUE)
# library(biomaRt)
library("org.Hs.eg.db")
# library(reticulate)
Sys.setenv(RETICULATE_PYTHON = python_path)
library(anndata)

args = commandArgs(trailingOnly = TRUE)
data_dir = args[1]
train_adata = read_h5ad(file.path(data_dir, "train_adata.h5ad"))
train_mapped_genes = ensembldb::select(org.Hs.eg.db, keys= train_adata$var_names, 
                    keytype = "SYMBOL", 
                    columns = c("SYMBOL","ENSEMBL"))
# deduplicate and remove na genes
match_idx = match(train_adata$var_names, train_mapped_genes$SYMBOL)
names(match_idx) = 1:length(train_adata$var_names)
match_idx = match_idx[!is.na(match_idx)] # remove NA entries
dedup_match_idx = match_idx[!duplicated(match_idx)]
ensembl_ids = train_mapped_genes[unname(dedup_match_idx), "ENSEMBL"]
train_adata$var['ensembl_id'] = ensembl_ids
nonna_ensembl_idx = which(!is.na(ensembl_ids))
write_h5ad(train_adata[, nonna_ensembl_idx], file.path(data_dir, "train_adata_mapped.h5ad"))

# repeat for test data
test_adata = read_h5ad(file.path(data_dir, "test_adata.h5ad"))
test_mapped_genes = ensembldb::select(org.Hs.eg.db, keys= test_adata$var_names, 
                    keytype = "SYMBOL", 
                    columns = c("SYMBOL","ENSEMBL"))
match_idx = match(test_adata$var_names, test_mapped_genes$SYMBOL)
names(match_idx) = 1:length(test_adata$var_names)
match_idx = match_idx[!is.na(match_idx)] # remove NA entries
dedup_match_idx = match_idx[!duplicated(match_idx)]
ensembl_ids = test_mapped_genes[unname(dedup_match_idx), "ENSEMBL"]
test_adata$var['ensembl_id'] = ensembl_ids
nonna_ensembl_idx = which(!is.na(ensembl_ids))
write_h5ad(test_adata[, nonna_ensembl_idx], file.path(data_dir, "test_adata_mapped.h5ad"))
