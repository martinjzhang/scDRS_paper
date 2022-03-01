#!/bin/bash
#SBATCH -n 1                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-01:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=64000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0         # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.out # Standard output
#SBATCH -e /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.err # Standard error

# I_LINE=$SLURM_ARRAY_TASK_ID
# read -a TEMP_LINE <<< $( head ./data_list.txt -n $I_LINE | tail -n -1)
# DNAME=${TEMP_LINE[0]}
# DPATH=${TEMP_LINE[1]}

# H5AD_FILE=${DPATH}
# OUT_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/gene_stats/${DNAME}.tsv

# python compute_data_stats.py \
#     --h5ad_file $H5AD_FILE\
#     --out_file $OUT_FILE
    
    
# Ayhan (log1p-ed)
I_LINE=10
read -a TEMP_LINE <<< $( head ./data_list.txt -n $I_LINE | tail -n -1)
DNAME=${TEMP_LINE[0]}
DPATH=${TEMP_LINE[1]}

H5AD_FILE=${DPATH}
OUT_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/gene_stats/${DNAME}.tsv

echo $H5AD_FILE
echo $OUT_FILE

python compute_data_stats.py \
    --h5ad_file $H5AD_FILE\
    --flag_filter False\
    --flag_raw_count False\
    --out_file $OUT_FILE