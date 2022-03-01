#!/bin/bash
#SBATCH -n 1                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-02:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=64000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.out # Standard output
#SBATCH -e /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.err # Standard error

# ZSCORE_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/gene_annotation/MAGMA-v108/MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt
# OUT_FOLDER=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/gs_file

# for weight in zscore uniform
# do
#     # Top nmax genes
#     for nmax in 100 500 1000 2000
#         do
#         scdrs munge-gs\
#             --zscore_file ${ZSCORE_FILE}\
#             --out-file ${OUT_FOLDER}/magma_10kb_top${nmax}_${weight}.all_traits.rv1.gs\
#             --weight ${weight}\
#             --n-max ${nmax}
#         done
    
#     # FDR<0.01, capping between 100 and 2000
#     scdrs munge-gs\
#         --zscore_file ${ZSCORE_FILE}\
#         --out-file ${OUT_FOLDER}/magma_10kb_fdr001_cap100n2000_${weight}.all_traits.rv1.gs\
#         --weight ${weight}\
#         --fdr 0.01\
#         --n-min 100\
#         --n-max 2000
    
#     # FWER<0.05, capping between 100 and 2000
#     scdrs munge-gs\
#         --zscore_file ${ZSCORE_FILE}\
#         --out-file ${OUT_FOLDER}/magma_10kb_fwer005_cap100n2000_${weight}.all_traits.rv1.gs\
#         --weight ${weight}\
#         --fwer 0.05\
#         --n-min 100\
#         --n-max 2000
# done

# ZSCORE_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/gene_annotation/MAGMA-v108/MAGMA_v108_GENE_10_ZSTAT_for_scDRS_subsampled_refmt.txt
# OUT_FOLDER=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/gs_file

# scdrs munge-gs\
#     --zscore_file ${ZSCORE_FILE}\
#     --out-file ${OUT_FOLDER}/magma_10kb_top1000_zscore.UKB_subsample.rv1.gs\
#     --weight zscore\
#     --n-max 1000
    
ZSCORE_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/gene_annotation/MAGMA-v108/MAGMA_v108_GENE_10_ZSTAT.UKB145K_26trait.txt
OUT_FOLDER=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/gs_file

scdrs munge-gs\
    --zscore_file ${ZSCORE_FILE}\
    --out-file ${OUT_FOLDER}/magma_10kb_top1000_zscore.UKB145K_26trait.rv1.gs\
    --weight zscore\
    --n-max 1000