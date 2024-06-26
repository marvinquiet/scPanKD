#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=ProjecTIL_pipeline
#SBATCH --nodes=1
#SBATCH --partition=day-long-cpu
#SBATCH --output=ProjecTIL_pipeline.out

source ~/.bashrc
conda activate Seurat

cd /panfs/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/pipelines

# === ProjecTIL each study -> HNSC, common celltypes
Rscript run_ProjecTIL_pipeline.R --experiment ProjecTILs_each_CD8_to_HNSC_CD8

# === Other studies -> inconsistent cell types
#Rscript run_ProjecTIL_pipeline.R --experiment Chu_CD4T_validation
#Rscript run_ProjecTIL_pipeline.R --experiment Chu_CD8T_validation
#Rscript run_ProjecTIL_pipeline.R --experiment Zheng_CD4_to_Chu_CD4
#Rscript run_ProjecTIL_pipeline.R --experiment GSE179994_CD8_to_HNSC_CD8
#Rscript run_ProjecTIL_pipeline.R --experiment ProjecTILs_CD8_to_HNSC_CD8
#Rscript run_ProjecTIL_pipeline.R --experiment GSE179994_CD8_to_Chu_CD8T
#Rscript run_ProjecTIL_pipeline.R --experiment ProjecTILs_CD8_to_Chu_CD8T

# === multibatch
#Rscript run_ProjecTIL_pipeline.R --experiment Chu_CD4T_multibatch_validation
#Rscript run_ProjecTIL_pipeline.R --experiment Chu_CD8T_multibatch_validation
#Rscript run_ProjecTIL_pipeline.R --experiment Chu_GSE179994_CD8T_multibatch_to_HNSC_CD8
#Rscript run_ProjecTIL_pipeline.R --experiment ProjecTILs_CD8_multibatch_to_HNSC_CD8
#Rscript run_ProjecTIL_pipeline.R --experiment ProjecTILs_CD8_multibatch_to_Chu_CD8T


