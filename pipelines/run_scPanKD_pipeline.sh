#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=scPanKD_pipeline
#SBATCH --nodes=1
#SBATCH --partition=mulan-gpu
#SBATCH --output=scPanKD_pipeline.out

source ~/.bashrc
conda activate ~/envs/spatialATAC

cd /net/mulan/home/wenjinma/projects/scPanKD/pipelines

# === ProjecTIL each study -> HNSC, common celltypes
python run_scPanKD_pipeline.py --experiment ProjecTILs_each_CD8_to_HNSC_CD8 # celltype
python run_scPanKD_pipeline.py --experiment ProjecTILs_CD8_multibatch_trimmed_to_HNSC_CD8 # celltype

# === Other studies -> inconsistent cell types
python run_scPanKD_pipeline.py --experiment Chu_CD4_multibatch_validation
python run_scPanKD_pipeline.py --experiment Chu_CD8_multibatch_validation
python run_scPanKD_pipeline.py --experiment ProjecTILs_CD8_multibatch_to_GSE179994_CD8
python run_scPanKD_pipeline.py --experiment Chu_CD8_multibatch_to_GSE179994_CD8
python run_scPanKD_pipeline.py --experiment ProjecTILs_CD8_multibatch_to_HNSC_CD8
python run_scPanKD_pipeline.py --experiment Chu_CD8_multibatch_to_HNSC_CD8
python run_scPanKD_pipeline.py --experiment ProjecTILs_CD8_multibatch_to_Chu_CD8
python run_scPanKD_pipeline.py --experiment Chu_CD8_multibatch_to_ProjecTILs_CD8
python run_scPanKD_pipeline.py --experiment Zheng_CD4_multibatch_to_Chu_CD4
python run_scPanKD_pipeline.py --experiment Chu_CD4_multibatch_to_Zheng_CD4