#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=Cellcano_pipeline
#SBATCH --nodes=1
#SBATCH --partition=largemem
#SBATCH --output=Cellcano_pipeline.out

source ~/.bashrc
conda activate Cellcano

cd /panfs/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/pipelines
data_dir="/projects/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/data"
result_dir="/projects/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/results/Cellcano"

# === ProjecTIL each study -> HNSC, common celltypes
exprs="ProjecTILs_each_CD8_to_HNSC_CD8"
studies=('EGAS00001004809' 'GSE123814' 'GSE139555' 'GSE159251' 'GSE176021' 'GSE179994' 'GSE180268' 'PRJNA705464')
for study in "${studies[@]}"
do 
    echo $study
    each_result_dir=$result_dir/$exprs/$study
    mkdir -p $each_result_dir
    Cellcano train -i $data_dir/ProjecTILs_CD8T/ProjecTIL_${study}_CD8T_trimmed \
        -m $data_dir/ProjecTILs_CD8T/ProjecTIL_${study}_CD8T_trimmed_metadata.csv \
        -o $each_result_dir \
        --prefix Cellcano_trained_
    Cellcano predict -i $data_dir/HNSC_CD8T/HNSC_CD8T \
        --trained_model $each_result_dir/Cellcano_trained_MLP_model \
        -o $each_result_dir \
        --prefix Cellcano_predicted_
done


# === Other studies -> inconsistent cell types
# need to edit column curated_celltype to celltype for running Cellcano...
#exprs="Chu_CD4T_validation"
#Cellcano train -i $data_dir/Chu_pancancer_CD4T/train_80_pancanT \
#    -m $data_dir/Chu_pancancer_CD4T/train_80_pancanT_metadata.csv \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_trained_
#Cellcano predict -i $data_dir/Chu_pancancer_CD4T/test_20_pancanT \
#    --trained_model $result_dir/$exprs/Cellcano_trained_MLP_model \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_predicted_
#
#exprs="Zheng_CD4_to_Chu_CD4"
#Cellcano train -i $data_dir/Zheng_CD4T/Zheng_CD4T \
#    -m $data_dir/Zheng_CD4T/Zheng_CD4T_metadata.csv \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_trained_
#Cellcano predict -i $data_dir/Chu_pancancer_CD4T/test_20_pancanT \
#    --trained_model $result_dir/$exprs/Cellcano_trained_MLP_model \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_predicted_

#exprs="Chu_CD8T_validation"
#Cellcano train -i $data_dir/Chu_pancancer_CD8T/train_80_pancanT \
#    -m $data_dir/Chu_pancancer_CD8T/train_80_pancanT_metadata.csv \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_trained_
#Cellcano predict -i $data_dir/Chu_pancancer_CD8T/test_20_pancanT \
#    --trained_model $result_dir/$exprs/Cellcano_trained_MLP_model \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_predicted_
#
#exprs="GSE179994_CD8_to_HNSC_CD8"
#Cellcano train -i $data_dir/GSE179994_CD8T/GSE179994_CD8T_ref \
#    -m $data_dir/GSE179994_CD8T/GSE179994_CD8T_ref_metadata.csv \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_trained_
#Cellcano predict -i $data_dir/HNSC_CD8T/HNSC_CD8T \
#    --trained_model $result_dir/$exprs/Cellcano_trained_MLP_model \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_predicted_
#
#exprs="ProjecTILs_CD8_to_HNSC_CD8"
#Cellcano train -i $data_dir/ProjecTILs_CD8T/ProjecTIL_CD8T \
#    -m $data_dir/ProjecTILs_CD8T/ProjecTIL_CD8T_metadata.csv \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_trained_
#Cellcano predict -i $data_dir/HNSC_CD8T/HNSC_CD8T \
#    --trained_model $result_dir/$exprs/Cellcano_trained_MLP_model \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_predicted_
#
#exprs="GSE179994_CD8_to_Chu_CD8T"
#mkdir $result_dir/$exprs
#Cellcano predict -i $data_dir/Chu_pancancer_CD8T/test_20_pancanT \
#    --trained_model $result_dir/GSE179994_CD8_to_HNSC_CD8/Cellcano_trained_MLP_model \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_predicted_
#
#exprs="ProjecTILs_CD8_to_Chu_CD8T"
#mkdir $result_dir/$exprs
#Cellcano predict -i $data_dir/Chu_pancancer_CD8T/test_20_pancanT \
#    --trained_model $result_dir/ProjecTILs_CD8_to_HNSC_CD8/Cellcano_trained_MLP_model \
#    -o $result_dir/$exprs \
#    --prefix Cellcano_predicted_

# === multibatch
#Rscript run_Seurat_pipeline.R --experiment Chu_CD4T_multibatch_validation
#Rscript run_Seurat_pipeline.R --experiment Chu_CD8T_multibatch_validation
#Rscript run_Seurat_pipeline.R --experiment Chu_GSE179994_CD8T_multibatch_to_HNSC_CD8
#Rscript run_Seurat_pipeline.R --experiment ProjecTILs_CD8_multibatch_to_HNSC_CD8
#Rscript run_Seurat_pipeline.R --experiment ProjecTILs_CD8_multibatch_to_Chu_CD8T

