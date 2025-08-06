##!/bin/bash
##SBATCH --time=24:00:00
##SBATCH --job-name=draw_Sankey
##SBATCH --nodes=1
##SBATCH --partition=day-long-cpu
##SBATCH -o draw_Sankey_pipeline.log

#source ~/.bashrc
#conda activate ~/envs/SpaDOT
#
#cd /net/mulan/home/wenjinma/projects/scPanKD/analysis

declare -a exprs=(
    #"Chu_CD4T_validation"
    #"Chu_CD8T_validation"
    # "Zheng_CD4_to_Chu_CD4"
    # "GSE179994_CD8_to_HNSC_CD8"
    # "ProjecTILs_CD8_to_HNSC_CD8"
    # "GSE179994_CD8_to_Chu_CD8T"
    # "ProjecTILs_CD8_to_Chu_CD8T"
    #"Chu_CD4T_multibatch_validation"
    #"Chu_CD8T_multibatch_validation"
    #"GSE179994_CD8T_multibatch_to_HNSC_CD8"
    #"ProjecTILs_CD8_multibatch_to_Chu_CD8"
    #"Zheng_CD4_multibatch_to_Chu_CD4"
    #"GSE179994_CD8T_multibatch_to_Chu_CD8"
    #"ProjecTILs_CD8_multibatch_to_HNSC_CD8"
    "Chu_CD4_multibatch_validation"
    "Chu_CD8_multibatch_validation"
    "ProjecTILs_CD8_multibatch_to_GSE179994_CD8"
    "Chu_CD8_multibatch_to_GSE179994_CD8"
    "ProjecTILs_CD8_multibatch_to_HNSC_CD8"
    "Chu_CD8_multibatch_to_HNSC_CD8"
    "ProjecTILs_CD8_multibatch_to_Chu_CD8"
    "Chu_CD8_multibatch_to_ProjecTILs_CD8"
    "Zheng_CD4_multibatch_to_Chu_CD4"
    "Chu_CD4_multibatch_to_Zheng_CD4"
)
declare -a methods=(
    # "Seurat"
    # "scType"
    "ProjecTIL"
    #"CellTypist"
    # "scPanKD-best"
    "scGPT"
)

for method in ${methods[@]}; do
    echo "Running Sankey analysis for $method"
    for pred_expr in ${exprs[@]}; do
        Rscript draw_Sankey.R --experiment $pred_expr --method $method
    done
done

