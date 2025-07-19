#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=SCimilarity_pipeline
#SBATCH --nodes=1
#SBATCH --partition=mulan-gpu
#SBATCH --output=SCimilarity_pipeline.out

source ~/.bashrc
conda activate /net/zootopia/disk1/wenjinma/envs/SCimilarity/

cd /net/mulan/home/wenjinma/projects/scPanKD/pipelines
data_dir="/net/zootopia/disk1/wenjinma/data/scPanKD_data"
result_dir="/net/mulan/home/wenjinma/projects/scPanKD/results/SCimilarity"

# # === ProjecTIL each study -> HNSC, common celltypes
# exprs="ProjecTILs_each_CD8_to_HNSC_CD8"
# studies=('EGAS00001004809' 'GSE123814' 'GSE139555' 'GSE159251' 'GSE176021' 'GSE179994' 'GSE180268' 'PRJNA705464')
# for study in "${studies[@]}"
# do 
#     echo $study
#     each_result_dir=$result_dir/$exprs/$study
#     mkdir -p $each_result_dir
#     Cellcano train -i $data_dir/ProjecTILs_CD8T/ProjecTIL_${study}_CD8T_trimmed \
#         -m $data_dir/ProjecTILs_CD8T/ProjecTIL_${study}_CD8T_trimmed_metadata.csv \
#         -o $each_result_dir \
#         --prefix Cellcano_trained_
#     Cellcano predict -i $data_dir/HNSC_CD8T/HNSC_CD8T \
#         --trained_model $each_result_dir/Cellcano_trained_MLP_model \
#         -o $each_result_dir \
#         --prefix Cellcano_predicted_
# done

python run_SCimilarity_pipeline.py --experiment Chu_CD4T_validation
python run_SCimilarity_pipeline.py --experiment Chu_CD8T_validation
python run_SCimilarity_pipeline.py --experiment Zheng_CD4_to_Chu_CD4
python run_SCimilarity_pipeline.py --experiment GSE179994_CD8_to_HNSC_CD8
python run_SCimilarity_pipeline.py --experiment ProjecTILs_CD8_to_HNSC_CD8
python run_SCimilarity_pipeline.py --experiment GSE179994_CD8_to_Chu_CD8T
python run_SCimilarity_pipeline.py --experiment ProjecTILs_CD8_to_Chu_CD8T


# === multibatch
#Rscript run_Seurat_pipeline.R --experiment Chu_CD4T_multibatch_validation
#Rscript run_Seurat_pipeline.R --experiment Chu_CD8T_multibatch_validation
#Rscript run_Seurat_pipeline.R --experiment Chu_GSE179994_CD8T_multibatch_to_HNSC_CD8
#Rscript run_Seurat_pipeline.R --experiment ProjecTILs_CD8_multibatch_to_HNSC_CD8
#Rscript run_Seurat_pipeline.R --experiment ProjecTILs_CD8_multibatch_to_Chu_CD8T

