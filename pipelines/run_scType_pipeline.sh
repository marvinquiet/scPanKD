#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=scType_pipeline
#SBATCH --nodes=1
#SBATCH --partition=day-long-cpu
#SBATCH -o scType_pipeline.log

source ~/.bashrc
conda activate Seurat

cd /panfs/compbio/users/wma36/test_Cellcano_XuefengWang/test_Cellcano_XuefengWang/pipelines

#Rscript run_scType_pipeline.R --experiment Chu_CD4T_validation
Rscript run_scType_pipeline.R --experiment Chu_CD8T_validation
Rscript run_scType_pipeline.R --experiment Zheng_CD4_to_Chu_CD4
Rscript run_scType_pipeline.R --experiment GSE179994_CD8_to_HNSC_CD8

# No need to check multibatch
#Rscript run_scType_pipeline.R --experiment Chu_CD4T_multibatch_validation
#Rscript run_scType_pipeline.R --experiment Chu_CD8T_multibatch_validation
