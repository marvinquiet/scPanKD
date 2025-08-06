# Jun 26

# conda activate ~/envs/spatialATAC
#for study in GSE159251
# for study in GSE176021 GSE123814 GSE179994 EGAS00001004809 PRJNA705464 GSE180268 GSE139555
# do
#     # step 1
#     python train.py extract_label \
#       --metadata_path ./datasets_common_celltype/ProjecTIL_${study}_CD8T_trimmed_metadata.csv \
#       --col curated_celltype
#     # step 2
#     python train.py do_train \
#       --mtx_path ./datasets_common_celltype/ProjecTIL_${study}_CD8T_trimmed.mtx.gz \
#       --label_path ./datasets_common_celltype/ProjecTIL_${study}_CD8T_trimmed_metadata_label.csv \
#       --gene_path ./datasets_common_celltype/ProjecTIL_${study}_CD8T_trimmed_genes.tsv \
#       --out_dir ./exps_common_celltype_Jun26/${study}
#     # step 3
#     python predict.py do_predict \
#       --mtx_path datasets_common_celltype/HNSC_CD8T/HNSC_CD8T.mtx.gz \
#       --gene_path datasets_common_celltype/HNSC_CD8T/HNSC_CD8T_genes.tsv \
#       --cell_path datasets_common_celltype/HNSC_CD8T/HNSC_CD8T_barcodes.tsv \
#       --model_name ./exps_common_celltype_Jun26/${study}/default_model.pkl  \
#       --out_dir ./exps_common_celltype_Jun26/${study}
# done

#for study in GSE159251
data_dir='/net/zootopia/disk1/wenjinma/data/scPanKD_data/'
result_dir='/net/mulan/home/wenjinma/projects/scPanKD/results/CellTypist/ProjecTILs_each_CD8_cancertrimmed_to_HNSC_CD8'
for study in HNSCC Melanoma SCC Lung Endometrial Renal Breast
do
    # step 1
    python train.py extract_label \
      --metadata_path ${data_dir}/ProjecTILs_CD8T/Tissue_trimmed/ProjecTIL_${study}_CD8T_trimmed_metadata.csv \
      --col celltype
    # step 2
    python train.py do_train \
      --mtx_path ${data_dir}/ProjecTILs_CD8T/Tissue_trimmed/ProjecTIL_${study}_CD8T_trimmed.mtx.gz \
      --label_path ${data_dir}/ProjecTILs_CD8T/Tissue_trimmed/ProjecTIL_${study}_CD8T_trimmed_metadata_label.csv \
      --gene_path ${data_dir}/ProjecTILs_CD8T/Tissue_trimmed/ProjecTIL_${study}_CD8T_trimmed_genes.tsv \
      --out_dir ${result_dir}/${study}
    # step 3
    python predict.py do_predict \
      --mtx_path ${data_dir}/HNSC_CD8T/HNSC_CD8T.mtx.gz \
      --gene_path ${data_dir}/HNSC_CD8T/HNSC_CD8T_genes.tsv \
      --cell_path ${data_dir}/HNSC_CD8T/HNSC_CD8T_barcodes.tsv \
      --model_name ${result_dir}/${study}/default_model.pkl  \
      --out_dir ${result_dir}/${study}
done