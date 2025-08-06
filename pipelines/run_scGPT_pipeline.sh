#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --job-name=scGPT_pipeline
#SBATCH --nodes=1
#SBATCH --partition=mulan-gpu
#SBATCH --output=scGPT_pipeline.out

source ~/.bashrc
conda activate ~/envs/scGPT_bak/

cd /net/mulan/home/wenjinma/projects/scPanKD/pipelines
data_dir="/net/zootopia/disk1/wenjinma/data/scPanKD_data"
result_dir="/net/mulan/home/wenjinma/projects/scPanKD/results/scGPT"


python run_scGPT_pipeline.py --experiment Chu_CD4_multibatch_validation
python run_scGPT_pipeline.py --experiment Chu_CD8_multibatch_validation
python run_scGPT_pipeline.py --experiment ProjecTILs_CD8_multibatch_to_GSE179994_CD8
python run_scGPT_pipeline.py --experiment Chu_CD8_multibatch_to_GSE179994_CD8
python run_scGPT_pipeline.py --experiment ProjecTILs_CD8_multibatch_to_HNSC_CD8
python run_scGPT_pipeline.py --experiment Chu_CD8_multibatch_to_HNSC_CD8
python run_scGPT_pipeline.py --experiment ProjecTILs_CD8_multibatch_to_Chu_CD8
python run_scGPT_pipeline.py --experiment Chu_CD8_multibatch_to_ProjecTILs_CD8
python run_scGPT_pipeline.py --experiment Zheng_CD4_multibatch_to_Chu_CD4
python run_scGPT_pipeline.py --experiment Chu_CD4_multibatch_to_Zheng_CD4