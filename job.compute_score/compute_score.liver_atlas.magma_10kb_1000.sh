#!/bin/bash
#SBATCH -n 1                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-01:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=16000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=1-14          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.out # Standard output
#SBATCH -e /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.err # Standard error



BATCH_NUM=$SLURM_ARRAY_TASK_ID
H5AD_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/liver_atlas/obj_raw.h5ad
GS_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/gs_file/magma_10kb_top1000_zscore.75_traits.batch/batch$BATCH_NUM.gs
OUT_FOLDER=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/score_file/score.liver_atlas.magma_10kb_top1000_zscore

scdrs compute-score \
    --h5ad_file $H5AD_FILE\
    --h5ad_species human\
    --gs_file $GS_FILE\
    --gs_species human\
    --flag_filter_data True\
    --flag_raw_count True\
    --n_ctrl 1000\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER
    
    
    
# BATCH_NUM=$SLURM_ARRAY_TASK_ID
# # BATCH_NUM=0
# H5AD_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/liver_atlas/obj_raw.h5ad
# GS_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/gs_file/magma_10kb_1000.74_traits.gs.batch/magma_10kb_1000.batch$BATCH_NUM.gs
# OUT_FOLDER=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/score_file/score.liver_atlas.magma_10kb_1000

# python /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/compute_score.py \
#     --h5ad_file $H5AD_FILE\
#     --h5ad_species human\
#     --gs_file $GS_FILE\
#     --gs_species human\
#     --flag_filter True\
#     --flag_raw_count True\
#     --n_ctrl 1000\
#     --flag_return_ctrl_raw_score False\
#     --flag_return_ctrl_norm_score True\
#     --out_folder $OUT_FOLDER