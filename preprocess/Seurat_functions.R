predict_Seurat = function(ref_obj, tgt_obj, ref_celltype_ind='curated_celltype') {
    k.filter = ifelse(ncol(tgt_obj) > 200, 200, NA)
    ## if from the same modality, use pcaproject; otherwise, use cca
    transfer.anchors <- FindTransferAnchors(reference = ref_obj, query = tgt_obj, dims = 1:30,
                                                reference.reduction = "pca", k.filter=k.filter)
    celltype.predictions = TransferData(anchorset = transfer.anchors, refdata = ref_obj@meta.data[, ref_celltype_ind],
                                             dims=1:30)
    tgt_obj = AddMetaData(tgt_obj, metadata=celltype.predictions)
    return(tgt_obj)
}

load_Seurat_data = function(experiment, seed=1234) {
    data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data"
    require(Seurat)
    require(Matrix)
    # if (experiment == 'Chu_CD4T_validation') {
    #     # NOTE: For Chu CD4T cell, the original data is given by log-normalized
    #     Chu_data_dir = file.path(data_dir, 'Chu.PancanerT_raw')
    #     load(file.path(Chu_data_dir, 'CD4.keep.RData'))
    #     pancancerT_obj = readRDS(file.path(Chu_data_dir, 'CD4'))
    #     CD4_keep_cells = rownames(CD4.keep)
    #     common_cells = intersect(colnames(pancancerT_obj), CD4_keep_cells) # 94571 cells
    #     pancancerT_obj = pancancerT_obj[, common_cells]
    #     pancancerT_obj$curated_celltype = CD4.keep[common_cells, 'Annotation']
    #     counts = GetAssayData(pancancerT_obj, slot="counts", assay="RNA")  # seems to be log-normalized
    #     genes_percent = rowMeans(counts > 0)
    #     genes_selected = genes_percent > 0.01 # remaining 12,191 genes
    #     pancancerT_obj = pancancerT_obj[genes_selected]
    #     # split by 80 as train, 20 as target
    #     set.seed(seed)
    #     train_idx = sample(1:ncol(pancancerT_obj), round(0.8*ncol(pancancerT_obj)))
    #     test_idx = setdiff(1:ncol(pancancerT_obj), train_idx)
    #     train_obj = pancancerT_obj[, train_idx]
    #     train_obj = FindVariableFeatures(train_obj)
    #     train_obj = ScaleData(train_obj)
    #     train_obj = RunPCA(train_obj)
    #     test_obj = pancancerT_obj[, test_idx] # already log-normalized
    # }
    # if (experiment == 'Chu_CD8T_validation') {
    #     # NOTE: For Chu CD8T cell, the original data is given by log-normalized
    #     Chu_data_dir = file.path(data_dir, 'Chu.PancanerT_raw')
    #     load(file.path(Chu_data_dir, 'CD8.keep.RData'))
    #     pancancerT_obj = readRDS(file.path(Chu_data_dir, 'CD8'))
    #     CD8_keep_cells = rownames(CD8.keep)
    #     common_cells = intersect(colnames(pancancerT_obj), CD8_keep_cells) # 41809 cells
    #     pancancerT_obj = pancancerT_obj[, common_cells]
    #     pancancerT_obj$curated_celltype = CD8.keep[common_cells, 'Annotation']
    #     counts = GetAssayData(pancancerT_obj, slot="counts", assay="RNA")  # seems to be log-normalized
    #     genes_percent = rowMeans(counts > 0)
    #     genes_selected = genes_percent > 0.01 # remaining 12,175 genes
    #     pancancerT_obj = pancancerT_obj[genes_selected]
    #     # split by 80 as train, 20 as target
    #     set.seed(seed)
    #     train_idx = sample(1:ncol(pancancerT_obj), round(0.8*ncol(pancancerT_obj)))
    #     test_idx = setdiff(1:ncol(pancancerT_obj), train_idx)
    #     train_obj = pancancerT_obj[, train_idx]
    #     train_obj = FindVariableFeatures(train_obj)
    #     train_obj = ScaleData(train_obj)
    #     train_obj = RunPCA(train_obj)
    #     test_obj = pancancerT_obj[, test_idx] # already log-normalized
    # }
    # if (experiment == 'Zheng_CD4_to_Chu_CD4') {
    #     Zheng_data_dir = file.path(data_dir, 'Zheng_CD4T')
    #     Zheng_obj = readRDS(file.path(Zheng_data_dir, 'Zheng.PancanerT_raw', 'CD4.rds'))
    #     selected_idx = Zheng_obj$Annotation %in% c('CD4.EOMES+', 'CD4.Naivelike', 'Th17', 'Tfh', 'Treg')
    #     Zheng_obj = Zheng_obj[,selected_idx] ## selected cell types, 10911 cells remained
    #     Zheng_obj$curated_celltype = Zheng_obj$Annotation
    #     integrated_counts = GetAssayData(Zheng_obj, slot="data", assay="integrated")   
    #     selected_genes = rownames(integrated_counts)
    #     counts = GetAssayData(Zheng_obj, slot="counts", assay="RNA")
    #     counts = counts[selected_genes,] ## 800 genes x 10911 cells
    #     train_obj = CreateSeuratObject(counts, meta.data=Zheng_obj@meta.data)
    #     train_obj = NormalizeData(train_obj)
    #     train_obj = FindVariableFeatures(train_obj)
    #     train_obj = ScaleData(train_obj)
    #     train_obj = RunPCA(train_obj)
    #     Chu_data_dir = file.path(data_dir, 'Chu.PancanerT_raw')
    #     load(file.path(Chu_data_dir, 'CD4.keep.RData'))
    #     pancancerT_obj = readRDS(file.path(Chu_data_dir, 'CD4'))
    #     CD4_keep_cells = rownames(CD4.keep)
    #     common_cells = intersect(colnames(pancancerT_obj), CD4_keep_cells) # 94571 cells
    #     pancancerT_obj = pancancerT_obj[, common_cells]
    #     pancancerT_obj$curated_celltype = CD4.keep[common_cells, 'Annotation']
    #     counts = GetAssayData(pancancerT_obj, slot="counts", assay="RNA")  # seems to be log-normalized
    #     genes_percent = rowMeans(counts > 0)
    #     genes_selected = genes_percent > 0.01 # remaining 12,191 genes
    #     pancancerT_obj = pancancerT_obj[genes_selected]
    #     # split by 80 as train, 20 as target
    #     set.seed(seed)
    #     train_idx = sample(1:ncol(pancancerT_obj), round(0.8*ncol(pancancerT_obj)))
    #     test_idx = setdiff(1:ncol(pancancerT_obj), train_idx)
    #     test_obj = pancancerT_obj[, test_idx] # already log-normalized
    #     #Chu_data_dir = file.path(data_dir, 'Chu.PancanerT_raw')
    #     #load(file.path(Chu_data_dir, 'CD4.keep.RData'))
    #     #pancancerT_obj = readRDS(file.path(Chu_data_dir, 'CD4'))
    #     #CD4_keep_cells = rownames(CD4.keep)
    #     #common_cells = intersect(colnames(pancancerT_obj), CD4_keep_cells) # 94571 cells
    #     #pancancerT_obj = pancancerT_obj[, common_cells]
    #     #pancancerT_obj$curated_celltype = CD4.keep[common_cells, 'Annotation']
    #     #idx = pancancerT_obj$curated_celltype == 'D4.Naivelike'
    #     #pancancerT_obj$curated_celltype[idx] = 'CD4.Naivelike'
    #     #selected_idx = pancancerT_obj$curated_celltype %in% c('CD4.CTL', 'CD4.Naivelike', 'Tfh', 'Th17', 'Treg') # 54079 cells
    #     #pancancerT_obj = pancancerT_obj[, selected_idx]
    #     #counts = GetAssayData(pancancerT_obj, slot="counts", assay="RNA")  # seems to be log-normalized
    #     #genes_percent = rowMeans(counts > 0)
    #     #genes_selected = genes_percent > 0.01 # remaining 12,191 genes
    #     #pancancerT_obj = pancancerT_obj[genes_selected]
    #     #test_obj = pancancerT_obj
    # }
    # if (experiment == 'GSE179994_CD8_to_HNSC_CD8') {
    #     GSE179994_data_dir = file.path(data_dir, 'GSE179994_CD8T')
    #     load(file.path(GSE179994_data_dir, 'raw', 'GSE179994.CD8T.ref.rawcount.RData'))
    #     rownames(GSE179994.CD8T.ref.label) = GSE179994.CD8T.ref.label$CellID
    #     colnames(GSE179994.CD8T.ref.label) = c('CellID', 'curated_celltype', 'Sample')
    #     GSE179994.CD8T.ref.label$batch = unlist(sapply(strsplit(GSE179994.CD8T.ref.label$CellID, split='_'), '[', 1))
    #     train_obj = CreateSeuratObject(counts=GSE179994.CD8T.ref.rawcount.mat, meta.data=GSE179994.CD8T.ref.label) # 27421 genes x 1694 cells
    #     train_obj = NormalizeData(train_obj)
    #     train_obj = FindVariableFeatures(train_obj)
    #     train_obj = ScaleData(train_obj)
    #     train_obj = RunPCA(train_obj)
    #     HNSC_data_dir = file.path(data_dir, 'HNSC_CD8T')
    #     test_counts = readMM(file.path(HNSC_data_dir, 'HNSC_CD8T.mtx.gz'))
    #     test_barcodes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_barcodes.tsv'), what=character())
    #     test_genes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_genes.tsv'), what=character())
    #     rownames(test_counts) = test_genes
    #     colnames(test_counts) = test_barcodes
    #     test_metadata = read.csv(file.path(HNSC_data_dir, 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
    #     test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata)
    #     test_obj = NormalizeData(test_obj)
    # }
    # if (experiment == 'ProjecTILs_CD8_to_HNSC_CD8') {
    #     ProjecTIL_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
    #     train_obj = readRDS(file.path(ProjecTIL_data_dir, 'raw', 'CD8T_human_ref_v1.rds'))
    #     DefaultAssay(train_obj) = 'RNA'
    #     train_obj$curated_celltype = train_obj$functional.cluster
    #     train_obj = NormalizeData(train_obj)
    #     train_obj = FindVariableFeatures(train_obj)
    #     train_obj = ScaleData(train_obj)
    #     train_obj = RunPCA(train_obj)
    #     HNSC_data_dir = file.path(data_dir, 'HNSC_CD8T')
    #     test_counts = readMM(file.path(HNSC_data_dir, 'HNSC_CD8T.mtx.gz'))
    #     test_barcodes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_barcodes.tsv'), what=character())
    #     test_genes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_genes.tsv'), what=character())
    #     rownames(test_counts) = test_genes
    #     colnames(test_counts) = test_barcodes
    #     test_metadata = read.csv(file.path(HNSC_data_dir, 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
    #     test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata)
    #     test_obj = NormalizeData(test_obj)
    # }
    # if (experiment == 'GSE179994_CD8_to_Chu_CD8T') {
    #     GSE179994_data_dir = file.path(data_dir, 'GSE179994_CD8T')
    #     load(file.path(GSE179994_data_dir, 'raw', 'GSE179994.CD8T.ref.rawcount.RData'))
    #     rownames(GSE179994.CD8T.ref.label) = GSE179994.CD8T.ref.label$CellID
    #     colnames(GSE179994.CD8T.ref.label) = c('CellID', 'curated_celltype', 'Sample')
    #     GSE179994.CD8T.ref.label$batch = unlist(sapply(strsplit(GSE179994.CD8T.ref.label$CellID, split='_'), '[', 1))
    #     train_obj = CreateSeuratObject(counts=GSE179994.CD8T.ref.rawcount.mat, meta.data=GSE179994.CD8T.ref.label) # 27421 genes x 1694 cells
    #     train_obj = NormalizeData(train_obj)
    #     train_obj = FindVariableFeatures(train_obj)
    #     train_obj = ScaleData(train_obj)
    #     train_obj = RunPCA(train_obj)
    #     Chu_data_dir = file.path(data_dir, 'Chu.PancanerT_raw')
    #     load(file.path(Chu_data_dir, 'CD8.keep.RData'))
    #     pancancerT_obj = readRDS(file.path(Chu_data_dir, 'CD8'))
    #     CD8_keep_cells = rownames(CD8.keep)
    #     common_cells = intersect(colnames(pancancerT_obj), CD8_keep_cells) # 41809 cells
    #     pancancerT_obj = pancancerT_obj[, common_cells]
    #     pancancerT_obj$curated_celltype = CD8.keep[common_cells, 'Annotation']
    #     counts = GetAssayData(pancancerT_obj, slot="counts", assay="RNA")  # seems to be log-normalized
    #     genes_percent = rowMeans(counts > 0)
    #     genes_selected = genes_percent > 0.01 # remaining 12,175 genes
    #     pancancerT_obj = pancancerT_obj[genes_selected]
    #     # split by 80 as train, 20 as target
    #     set.seed(seed)
    #     train_idx = sample(1:ncol(pancancerT_obj), round(0.8*ncol(pancancerT_obj)))
    #     test_idx = setdiff(1:ncol(pancancerT_obj), train_idx)
    #     test_obj = pancancerT_obj[, test_idx] # already log-normalized
    # }
    # if (experiment == 'ProjecTILs_CD8_to_Chu_CD8T') {
    #     ProjecTIL_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
    #     train_obj = readRDS(file.path(ProjecTIL_data_dir, 'raw', 'CD8T_human_ref_v1.rds'))
    #     DefaultAssay(train_obj) = 'RNA'
    #     train_obj$curated_celltype = train_obj$functional.cluster
    #     train_obj = NormalizeData(train_obj)
    #     train_obj = FindVariableFeatures(train_obj)
    #     train_obj = ScaleData(train_obj)
    #     train_obj = RunPCA(train_obj)
    #     Chu_data_dir = file.path(data_dir, 'Chu.PancanerT_raw')
    #     load(file.path(Chu_data_dir, 'CD8.keep.RData'))
    #     pancancerT_obj = readRDS(file.path(Chu_data_dir, 'CD8'))
    #     CD8_keep_cells = rownames(CD8.keep)
    #     common_cells = intersect(colnames(pancancerT_obj), CD8_keep_cells) # 41809 cells
    #     pancancerT_obj = pancancerT_obj[, common_cells]
    #     pancancerT_obj$curated_celltype = CD8.keep[common_cells, 'Annotation']
    #     counts = GetAssayData(pancancerT_obj, slot="counts", assay="RNA")  # seems to be log-normalized
    #     genes_percent = rowMeans(counts > 0)
    #     genes_selected = genes_percent > 0.01 # remaining 12,175 genes
    #     pancancerT_obj = pancancerT_obj[genes_selected]
    #     # split by 80 as train, 20 as target
    #     set.seed(seed)
    #     train_idx = sample(1:ncol(pancancerT_obj), round(0.8*ncol(pancancerT_obj)))
    #     test_idx = setdiff(1:ncol(pancancerT_obj), train_idx)
    #     test_obj = pancancerT_obj[, test_idx] # already log-normalized
    # }

    if (experiment == 'ProjecTILs_CD8_multibatch_trimmed_to_HNSC_CD8') {
        # --- load train dataset
        ProjecTIL_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
        train_obj = readRDS(file.path(ProjecTIL_data_dir, 'ProjecTILs_CD8T_trimmed_Cohort_integrated.RDS'))
        # --- load test dataset
        HNSC_data_dir = file.path(data_dir, 'HNSC_CD8T')
        test_counts = readMM(file.path(HNSC_data_dir, 'HNSC_CD8T.mtx.gz'))
        test_barcodes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_barcodes.tsv'), what=character())
        test_genes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_genes.tsv'), what=character())
        rownames(test_counts) = test_genes
        colnames(test_counts) = test_barcodes
        test_metadata = read.csv(file.path(HNSC_data_dir, 'HNSC_CD8T_metadata_curated.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata)
    }

    # === multiple batches
    if (experiment == 'Chu_CD4_multibatch_validation') {
        # NOTE: For Chu CD4T cell, the original data is given by log-normalized
        Chu_CD4_data_dir = file.path(data_dir, 'Chu_pancancer_CD4T')
        train_obj = readRDS(file.path(Chu_CD4_data_dir, 'Chu_CD4T_train_80_integrated.RDS'))
        test_counts = readMM(file.path(Chu_CD4_data_dir, 'test_20_pancanT.mtx.gz'))
        test_genes = scan(file.path(Chu_CD4_data_dir, 'test_20_pancanT_genes.tsv'), what=character())
        test_barcodes = scan(file.path(Chu_CD4_data_dir, 'test_20_pancanT_barcodes.tsv'), what=character())
        rownames(test_counts) = test_genes
        colnames(test_counts) = test_barcodes
        test_metadata = read.csv(file.path(Chu_CD4_data_dir, 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata[colnames(test_counts), ])
    }
    if (experiment == 'Chu_CD8_multibatch_validation') {
        # NOTE: For Chu CD8T cell, the original data is given by log-normalized
        Chu_CD8_data_dir = file.path(data_dir, 'Chu_pancancer_CD8T')
        train_obj = readRDS(file.path(Chu_CD8_data_dir, 'Chu_CD8T_train_80_integrated.RDS'))
        test_counts = readMM(file.path(Chu_CD8_data_dir, 'test_20_pancanT.mtx.gz'))
        test_genes = scan(file.path(Chu_CD8_data_dir, 'test_20_pancanT_genes.tsv'), what=character())
        test_barcodes = scan(file.path(Chu_CD8_data_dir, 'test_20_pancanT_barcodes.tsv'), what=character())
        rownames(test_counts) = test_genes
        colnames(test_counts) = test_barcodes
        test_metadata = read.csv(file.path(Chu_CD8_data_dir, 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata[colnames(test_counts), ])
    }

    if (experiment == 'ProjecTILs_CD8_multibatch_to_GSE179994_CD8') {
        # --- load train dataset
        ProjecTIL_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
        train_obj = readRDS(file.path(ProjecTIL_data_dir, 'ProjecTILs_CD8T_Tissue_integrated.RDS'))
        # --- load test dataset
        GSE179994_data_dir = file.path(data_dir, 'GSE179994_CD8T')
        GSE179994_counts = readMM(file.path(GSE179994_data_dir, 'GSE179994_CD8T_ref.mtx.gz'))
        GSE179994_genes = scan(file.path(GSE179994_data_dir, 'GSE179994_CD8T_ref_genes.tsv'), what=character())
        GSE179994_barcodes = scan(file.path(GSE179994_data_dir, 'GSE179994_CD8T_ref_barcodes.tsv'), what=character())
        rownames(GSE179994_counts) = GSE179994_genes
        colnames(GSE179994_counts) = GSE179994_barcodes
        GSE179994_metadata = read.csv(file.path(GSE179994_data_dir, 'GSE179994_CD8T_ref_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=GSE179994_counts, meta.data=GSE179994_metadata[colnames(GSE179994_counts), ])
    }

    if (experiment == 'Chu_CD8_multibatch_to_GSE179994_CD8') {
        # --- load train dataset
        Chu_CD8_data_dir = file.path(data_dir, 'Chu_pancancer_CD8T')
        train_obj = readRDS(file.path(Chu_CD8_data_dir, 'Chu_CD8T_train_80_integrated.RDS'))
        # --- load test dataset
        GSE179994_data_dir = file.path(data_dir, 'GSE179994_CD8T')
        GSE179994_counts = readMM(file.path(GSE179994_data_dir, 'GSE179994_CD8T_ref.mtx.gz'))
        GSE179994_genes = scan(file.path(GSE179994_data_dir, 'GSE179994_CD8T_ref_genes.tsv'), what=character())
        GSE179994_barcodes = scan(file.path(GSE179994_data_dir, 'GSE179994_CD8T_ref_barcodes.tsv'), what=character())
        rownames(GSE179994_counts) = GSE179994_genes
        colnames(GSE179994_counts) = GSE179994_barcodes
        GSE179994_metadata = read.csv(file.path(GSE179994_data_dir, 'GSE179994_CD8T_ref_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=GSE179994_counts, meta.data=GSE179994_metadata[colnames(GSE179994_counts), ])
    }
    
    if (experiment == "ProjecTILs_CD8_multibatch_to_HNSC_CD8") {
        # --- load train dataset
        ProjecTIL_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
        train_obj = readRDS(file.path(ProjecTIL_data_dir, 'ProjecTILs_CD8T_Tissue_integrated.RDS'))
        # --- load test dataset
        HNSC_data_dir = file.path(data_dir, 'HNSC_CD8T')
        test_counts = readMM(file.path(HNSC_data_dir, 'HNSC_CD8T.mtx.gz'))
        test_barcodes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_barcodes.tsv'), what=character())
        test_genes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_genes.tsv'), what=character())
        rownames(test_counts) = test_genes
        colnames(test_counts) = test_barcodes
        test_metadata = read.csv(file.path(HNSC_data_dir, 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata)
        test_obj = NormalizeData(test_obj)
    }

    if (experiment == 'Chu_CD8_multibatch_to_HNSC_CD8') {
        # --- load train dataset
        Chu_data_dir = file.path(data_dir, 'Chu_pancancer_CD8T')
        train_obj = readRDS(file.path(Chu_data_dir, 'Chu_CD8T_train_80_integrated.RDS'))
        # --- load test dataset
        HNSC_data_dir = file.path(data_dir, 'HNSC_CD8T')
        test_counts = readMM(file.path(HNSC_data_dir, 'HNSC_CD8T.mtx.gz'))
        test_barcodes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_barcodes.tsv'), what=character())
        test_genes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_genes.tsv'), what=character())
        rownames(test_counts) = test_genes
        colnames(test_counts) = test_barcodes
        test_metadata = read.csv(file.path(HNSC_data_dir, 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata)
    }

    if (experiment == 'ProjecTILs_CD8_multibatch_to_Chu_CD8') {
        # --- load train dataset
        ProjecTIL_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
        train_obj = readRDS(file.path(ProjecTIL_data_dir, 'ProjecTILs_CD8T_Tissue_integrated.RDS'))
        # --- load test dataset
        test_counts = readMM(file.path(Chu_CD8_data_dir, 'test_20_pancanT.mtx.gz'))
        test_genes = scan(file.path(Chu_CD8_data_dir, 'test_20_pancanT_genes.tsv'), what=character())
        test_barcodes = scan(file.path(Chu_CD8_data_dir, 'test_20_pancanT_barcodes.tsv'), what=character())
        rownames(test_counts) = test_genes
        colnames(test_counts) = test_barcodes
        test_metadata = read.csv(file.path(Chu_CD8_data_dir, 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata[colnames(test_counts), ])
    }

    if (experiment == 'Chu_CD8_multibatch_to_ProjecTILs_CD8') {
        # --- load train dataset
        Chu_data_dir = file.path(data_dir, 'Chu_pancancer_CD8T')
        train_obj = readRDS(file.path(Chu_data_dir, 'Chu_CD8T_train_80_integrated.RDS'))
        # --- load test dataset
        ProjecTIL_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
        test_counts = readMM(file.path(ProjecTIL_data_dir, 'ProjecTILs_CD8T.mtx.gz'))
        test_barcodes = scan(file.path(ProjecTIL_data_dir, 'ProjecTILs_CD8T_barcodes.tsv'), what=character())
        test_genes = scan(file.path(ProjecTIL_data_dir, 'ProjecTILs_CD8T_genes.tsv'), what=character())
        rownames(test_counts) = test_genes
        colnames(test_counts) = test_barcodes
        test_metadata = read.csv(file.path(ProjecTIL_data_dir, 'ProjecTILs_CD8T_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata[colnames(test_counts), ])
    }

    if (experiment == 'Zheng_CD4_multibatch_to_Chu_CD4') {
        # --- load train dataset
        Zheng_data_dir = file.path(data_dir, 'Zheng_CD4T')
        train_obj = readRDS(file.path(Zheng_data_dir, 'Zheng_CD4T_cancerType_integrated.RDS'))
        # --- load test dataset
        Chu_data_dir = file.path(data_dir, 'Chu_pancancer_CD4T')
        test_counts = readMM(file.path(Chu_data_dir, 'test_20_pancanT.mtx.gz'))
        test_genes = scan(file.path(Chu_data_dir, 'test_20_pancanT_genes.tsv'), what=character())
        test_barcodes = scan(file.path(Chu_data_dir, 'test_20_pancanT_barcodes.tsv'), what=character())
        rownames(test_counts) = test_genes
        colnames(test_counts) = test_barcodes
        test_metadata = read.csv(file.path(Chu_data_dir, 'test_20_pancanT_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata[colnames(test_counts), ])
    }

    if (experiment == 'Chu_CD4_multibatch_to_Zheng_CD4') {
        # --- load train dataset
        Chu_data_dir = file.path(data_dir, 'Chu_pancancer_CD4T')
        train_obj = readRDS(file.path(Chu_data_dir, 'Chu_CD4T_train_80_integrated.RDS'))
        # --- load test dataset
        Zheng_data_dir = file.path(data_dir, 'Zheng_CD4T')
        test_counts = readMM(file.path(Zheng_data_dir, 'Zheng_CD4T.mtx.gz'))
        test_barcodes = scan(file.path(Zheng_data_dir, 'Zheng_CD4T_barcodes.tsv'), what=character())
        test_genes = scan(file.path(Zheng_data_dir, 'Zheng_CD4T_genes.tsv'), what=character())
        rownames(test_counts) = test_genes
        colnames(test_counts) = test_barcodes
        test_metadata = read.csv(file.path(Zheng_data_dir, 'Zheng_CD4T_metadata.csv'), header=T, row.names=1)
        test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata[colnames(test_counts), ])
    }

    # if (experiment == 'Chu_GSE179994_CD8T_multibatch_to_HNSC_CD8') {
    #     # NOTE: For Chu CD8T cell, the original data is given by log-normalized
    #     Chu_data_dir = file.path(data_dir, 'Chu.PancanerT_raw')
    #     load(file.path(Chu_data_dir, 'CD8.keep.RData'))
    #     pancancerT_obj = readRDS(file.path(Chu_data_dir, 'CD8'))
    #     CD8_keep_cells = rownames(CD8.keep)
    #     common_cells = intersect(colnames(pancancerT_obj), CD8_keep_cells) # 41809 cells
    #     pancancerT_obj = pancancerT_obj[, common_cells]
    #     pancancerT_obj$curated_celltype = CD8.keep[common_cells, 'Annotation']
    #     counts = GetAssayData(pancancerT_obj, slot="counts", assay="RNA")  # seems to be log-normalized
    #     genes_percent = rowMeans(counts > 0)
    #     genes_selected = genes_percent > 0.01 # remaining 12,175 genes
    #     pancancerT_obj = pancancerT_obj[genes_selected]
    #     GSE179994_data_dir = file.path(data_dir, 'GSE179994_CD8T')
    #     load(file.path(GSE179994_data_dir, 'raw', 'GSE179994.CD8T.ref.rawcount.RData'))
    #     rownames(GSE179994.CD8T.ref.label) = GSE179994.CD8T.ref.label$CellID
    #     colnames(GSE179994.CD8T.ref.label) = c('CellID', 'curated_celltype', 'Sample')
    #     GSE179994.CD8T.ref.label$batch = unlist(sapply(strsplit(GSE179994.CD8T.ref.label$CellID, split='_'), '[', 1))
    #     GSE179994_obj = CreateSeuratObject(counts=GSE179994.CD8T.ref.rawcount.mat, meta.data=GSE179994.CD8T.ref.label) # 27421 genes x 1694 cells
    #     GSE179994_obj = NormalizeData(GSE179994_obj)
    #     # integrate reference dataset
    #     pancancerT_obj$dataset = 'Chu'
    #     GSE179994_obj$dataset = 'GSE179994'
    #     train_obj_list = list('Chu'=pancancerT_obj, 'GSE179994'=GSE179994_obj)
    #     train_obj_list = lapply(X = train_obj_list, FUN = function(x) {
    #                         x=FindVariableFeatures(x)
    #                      })
    #     features = SelectIntegrationFeatures(object.list=train_obj_list)
    #     anchors = FindIntegrationAnchors(object.list = train_obj_list, anchor.features = features)
    #     train_obj_integrated = IntegrateData(anchorset = anchors)
    #     # integrate trainining object using tissue as batch
    #     train_obj_integrated = ScaleData(train_obj_integrated)
    #     train_obj_integrated = RunPCA(train_obj_integrated)
    #     train_obj = train_obj_integrated
    #     HNSC_data_dir = file.path(data_dir, 'HNSC_CD8T')
    #     test_counts = readMM(file.path(HNSC_data_dir, 'HNSC_CD8T.mtx.gz'))
    #     test_barcodes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_barcodes.tsv'), what=character())
    #     test_genes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_genes.tsv'), what=character())
    #     rownames(test_counts) = test_genes
    #     colnames(test_counts) = test_barcodes
    #     test_metadata = read.csv(file.path(HNSC_data_dir, 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
    #     test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata)
    #     test_obj = NormalizeData(test_obj)
    # }
    # if (experiment == 'ProjecTILs_CD8_multibatch_to_HNSC_CD8') {
    #     ProjecTIL_data_dir = file.path(data_dir, 'ProjecTILs_CD8T')
    #     train_obj = readRDS(file.path(ProjecTIL_data_dir, 'raw', 'CD8T_human_ref_v1.rds'))
    #     DefaultAssay(train_obj) = 'RNA'
    #     train_obj$curated_celltype = train_obj$functional.cluster
    #     # perform integration using Seurat
    #     counts = train_obj@assays$RNA@counts
    #     genes_percent = rowMeans(counts > 0)
    #     genes_selected = genes_percent > 0.01 # remaining 10,116 genes
    #     train_obj = train_obj[genes_selected]
    #     train_obj$Cohort = as.factor(train_obj$Cohort) 
    #     train_obj_list = SplitObject(train_obj, split.by = "Cohort")
    #     train_obj_list = lapply(X = train_obj_list, FUN = function(x) {
    #                         x=FindVariableFeatures(x)
    #                      })
    #     features = SelectIntegrationFeatures(object.list=train_obj_list)
    #     anchors = FindIntegrationAnchors(object.list = train_obj_list, anchor.features = features)
    #     train_obj_integrated = IntegrateData(anchorset = anchors)
    #     # integrate trainining object using tissue as batch
    #     train_obj_integrated = ScaleData(train_obj_integrated)
    #     train_obj_integrated = RunPCA(train_obj_integrated)
    #     train_obj = train_obj_integrated
    #     HNSC_data_dir = file.path(data_dir, 'HNSC_CD8T')
    #     test_counts = readMM(file.path(HNSC_data_dir, 'HNSC_CD8T.mtx.gz'))
    #     test_barcodes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_barcodes.tsv'), what=character())
    #     test_genes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_genes.tsv'), what=character())
    #     rownames(test_counts) = test_genes
    #     colnames(test_counts) = test_barcodes
    #     test_metadata = read.csv(file.path(HNSC_data_dir, 'HNSC_CD8T_metadata.csv'), header=T, row.names=1)
    #     test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata)
    #     test_obj = NormalizeData(test_obj)
    # }
    test_obj = NormalizeData(test_obj) # normalize test data
    return(list(train_obj=train_obj, test_obj=test_obj))
}

## --- Seurat for ProjecTIL each to HNSC data
load_Seurat_data_for_ProjecTIL_to_HNSC = function(ref_study) {
    data_dir = "/net/zootopia/disk1/wenjinma/data/scPanKD_data"
    require(Seurat)
    require(Matrix)
    # --- load train data
    ProjecTIL_data_dir = file.path(data_dir, 'ProjecTILs_CD8T/trimmed')
    train_counts = readMM(file.path(ProjecTIL_data_dir, paste('ProjecTIL', ref_study, 'CD8T_trimmed.mtx.gz', sep='_')))
    train_genes = scan(file.path(ProjecTIL_data_dir, paste('ProjecTIL', ref_study, 'CD8T_trimmed_genes.tsv', sep='_')), what=character())
    train_barcodes = scan(file.path(ProjecTIL_data_dir, paste('ProjecTIL', ref_study, 'CD8T_trimmed_barcodes.tsv', sep='_')), what=character())
    rownames(train_counts) = train_genes
    colnames(train_counts) = train_barcodes
    train_metadata = read.csv(file.path(ProjecTIL_data_dir, paste('ProjecTIL', ref_study, 'CD8T_trimmed_metadata.csv', sep='_')), header=T, row.names=1)
    # --- load test data
    # HNSC_data_dir = file.path(data_dir, 'HNSC_CD8T')
    # test_counts = readMM(file.path(HNSC_data_dir, 'HNSC_CD8T.mtx.gz'))
    # test_barcodes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_barcodes.tsv'), what=character())
    # test_genes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_genes.tsv'), what=character())
    # rownames(test_counts) = test_genes
    # colnames(test_counts) = test_barcodes
    # test_metadata = read.csv(file.path(HNSC_data_dir, 'HNSC_CD8T_metadata_curated.csv'), header=T, row.names=1)
    # # --- find common cell types, no need to do it anymore
    # common_celltypes = unique(intersect(train_metadata$curated_celltype, test_metadata$celltype))
    # train_idx = which(train_metadata$curated_celltype %in% common_celltypes)
    # train_counts = train_counts[, train_idx]
    # train_metadata = train_metadata[train_idx, ]
    # test_idx = which(test_metadata$celltype %in% common_celltypes)
    # test_counts = test_counts[, test_idx]
    # test_metadata = test_metadata[test_idx, ]
    train_obj = CreateSeuratObject(counts=train_counts, meta.data=train_metadata)
    train_obj = NormalizeData(train_obj)
    train_obj = FindVariableFeatures(train_obj)
    train_obj = ScaleData(train_obj)
    train_obj = RunPCA(train_obj)
    HNSC_data_dir = file.path(data_dir, 'HNSC_CD8T')
    test_counts = readMM(file.path(HNSC_data_dir, 'HNSC_CD8T.mtx.gz'))
    test_barcodes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_barcodes.tsv'), what=character())
    test_genes = scan(file.path(HNSC_data_dir, 'HNSC_CD8T_genes.tsv'), what=character())
    rownames(test_counts) = test_genes
    colnames(test_counts) = test_barcodes
    test_metadata = read.csv(file.path(HNSC_data_dir, 'HNSC_CD8T_metadata_curated.csv'), header=T, row.names=1)
    test_obj = CreateSeuratObject(counts=test_counts, meta.data=test_metadata)
    test_obj = NormalizeData(test_obj)
    return(list(train_obj=train_obj, test_obj=test_obj))
}