Code and data for the paper Zhang*, Hou*, et al. [Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data](https://www.biorxiv.org/content/10.1101/2021.09.24.461597v2).

Also check out the scDRS [software](https://github.com/martinjzhang/scDRS) and [documentation](https://martinjzhang.github.io/scDRS/). 


### Versions
- [revision 1 (rv1)](https://www.biorxiv.org/content/10.1101/2021.09.24.461597v2): current version.
- [Initial submission](https://www.biorxiv.org/content/10.1101/2021.09.24.461597v1): see [readme_inital_sub.md](https://github.com/martinjzhang/scDRS_paper/blob/master/readme_initial_sub.md) file.

# Subset of code and data to reproduce main results of the paper

Codes are at `./job.reproduce`. To run the code, set `DATA_PATH` (if in the code) to your local folder of the main data `scDRS_data_release_XXXXXX` and set `SCORE_FILE_PATH` (if in the code) to your local folder of the scDRS score files `scDRS_data_release_XXXXXX.score_file_tmsfacs`.
- Revision 1 (rv1): Main data [scDRS_data_release_030122](https://figshare.com/articles/dataset/scDRS_data_release_030122/19312583) (3.8 GB). scDRS score files for TMS FACS + 74 diseases/traits [scDRS_data_release_030122.score_file_tmsfacs](https://figshare.com/articles/dataset/scDRS_data_release_030122_score_file_tmsfacs/19312607) (36.8 GB).
- Initial submission: Main data [scDRS_data_release_092121](https://figshare.com/articles/dataset/scDRS_data_release_092121/16664080) (3.6 GB). scDRS score files for TMS FACS + 74 diseases/traits [scDRS_data_release_092121.score_file_tmsfacs](https://figshare.com/articles/dataset/scDRS_data_release_092121_score_file_tmsfacs/16664077) (36.3 GB).


**Compute scDRS scores for TMS FACS + 74 diseases/traits**
- Revision 1 (rv1): Score files are in `scDRS_data_release_030122.score_file_tmsfacs`. Compute them yourself using `reproduce_compute_score.tms_facs_with_cov.magma_10kb_1000.rv1.sh`
- Initial submission: Score files were included in `scDRS_data_release_092121.score_file_tmsfacs`. Compute them yourself using `reproduce_compute_score.tms_facs_with_cov.magma_10kb_1000.sh`

**Cell type-level analysis (Fig. 3)**
- Revision 1 (rv1): `reproduce_celltype.rv1.ipynb`.
- Initial submission: `reproduce_celltype.ipynb`.

**T cell analysis (Fig. 4)**
- Revision 1 (rv1) (Fig. 4A-E): `reproduce_tcell.rv1.ipynb`.
- Initial submission (Fig. 4A-C): `reproduce_tcell.ipynb`.

**T cell gene prioritization (Fig. 4)**
- Revision 1 (rv1) (Fig. 4F): `reproduce_tcell_gene.rv1.ipynb`.
- Initial submission (Fig. 4D): `reproduce_tcell_gene.ipynb`.

**Neuron analysis (Fig. 5A,B)**
- Revision 1 (rv1): `reproduce_neuron.rv1.ipynb`.
- Initial submission: `reproduce_neuron.ipynb`.

**Hepatocyte analysis (Fig. 5C,D)**
- Revision 1 (rv1): `reproduce_hep.rv1.ipynb`.
- Initial submission: `reproduce_hep.ipynb`.



# Complete code

## Data curation: `job.curate_data`
Curate information for 74 diseases/traits: 
- Curate information for the 74 diseases: `job.curate_data/get_trait_list.ipynb`
- Curate information for the 74 diseases (rv1): `job.curate_data/get_trait_list.rv1.ipynb`

Curate gene set (.gs) files:
- .gs file for 74 diseases: `job.curate_data/curate_gs_file.ipynb`
- .gs file for 74 diseases (rv1): `job.curate_data/curate_gs_file.rv1.ipynb`
- .gs file for T cell signatures: `job.curate_data/curate_gs.tcell_signature.ipynb`
- .gs file for ploidy signatures: `job.curate_data/curate_ploidy_gs.ipynb`
- .gs file for zonation signatures: `job.curate_data/curate_zonation_gs.ipynb`
- .gs file for metabolic pathways: `job.curate_data/curate_gs.metabolic.ipynb`

Curate scRNA-seq data sets:
- TS FACS: `job.curate_data/curate_ts_data.ipynb`
- Cano-Gamez & Soskic et al.: `job.curate_data/curate_canogamez_tcell_data.ipynb`
- Nathan et al.: `job.curate_data/curate_nathan_tcell_data.ipynb`
- Aizarani et al.: `job.curate_data/curate_aizarani_liver_atlas_data.ipynb`
- Halpern & Shenhav et al.: `job.curate_data/curate_halpern_mouse_liver_data.ipynb`
- Richter & Deligiannis et al.: `job.curate_data/curate_richter_hepatocyte_data.ipynb`
- Compute gene-level statistics: `compute_data_stats.py` and `compute_data_stats.sh`
- Get meta information of data sets: `get_data_info.ipynb`
- Get meta information of data sets: `get_data_info.rv1.ipynb`


## Compute scDRS scores: `job.compute_score`
- TMS FACS + 74 diseases: `job.compute_score/compute_score.tms_facs_with_cov.magma_10kb_1000.sh`
- TMS FACS + T cell signatures (using scDRS scripts instead of CLI): `job.compute_score/compute_score.tms_facs_with_cov.tcell_sig.sh`
- TMS FACS + metabolic (using scDRS scripts instead of CLI): `job.compute_score/compute_score.tms_facs_with_cov.hep_metabolic.sh`
- TMS droplet + 74 diseases: `job.compute_score/compute_score.tms_droplet_with_cov.magma_10kb_1000.sh`
- TS FACS + 74 diseases: `job.compute_score/compute_score.ts_facs_with_cov.magma_10kb_1000.sh`


## Scehma (Fig. 1): `job.schema`
Make schematic figures.

## Simulation (Fig. 2): `job.simulation`
Data generation:
- Generate the TMS FACS 10K subsampled data and null gene sets: `job.simulation/generate_null_simulation_data.ipynb`
- Generate the TMS FACS 10K subsampled data and null gene sets (rv1): `job.simulation/generate_null_simulation_data.rv1.ipynb`
- Generate location-matched gene sets and make figures (rv1): `job.simulation/simulation.other_null.ipynb`
- Generate causal gene sets and perturbation configurations: `job.simulation/generate_causal_simulation_data.ipynb`
- Generate 20 reps of subsampled TMS FACS 10K data (rv1): `job.simulation/generate_subsampled_tms_data.rv1.ipynb`

Compute results: 
- Compute scDRS scores for null simulations: `job.simulation/compute_simu_score.sh`
- Compute scDRS scores for null simulations (rv1): `job.simulation/compute_simu_score.rv1.sh`
- Compute scDRS scores for null simulations using the `-adj-prop` option (rv1): `job.simulation/compute_simu_score.adj_prop.rv1.sh`
- Compute Seurat scores for null simulations: `job.simulation/compute_simu_score_scanpy.sh`
- Compute Seurat scores for null simulations (rv1): `job.simulation/compute_simu_score_scanpy.rv1.sh`
- Compute Vision scores for null simulations: `job.simulation/compute_simu_score_vision.sh`
- Compute Vision scores for null simulations (rv1): `job.simulation/compute_simu_score_vision.rv1.sh`
- Compute VAM scores for null simulations: `job.simulation/call_R_vam.sh`
- Compute VAM scores for null simulations (rv1): `job.simulation/call_R_vam.rv1.sh`
- Compute scores (scDRS/Seurat/Vision) for causal simulations (500 random causal cells): `job.simulation/compute_perturb_simu_score.sh`
- Compute scores (scDRS/Seurat/Vision) for causal simulations (B cells causal): `job.simulation/compute_perturb_simu_score_Bcell.sh`

Make figures:
- Make figures for null simulations: `job.simulation/make_figure.null_simulation.ipynb`
- Make figures for null simulations (rv1): `job.simulation/make_figure.null_simulation.rv1.ipynb`
- Make figures for causal simulations (500 random causal cells): `job.simulation/make_figure.causal_simulation.ipynb`
- Make figures for causal simulations (B cells causal): `job.simulation/make_figure.causal_simulation_Bcell.ipynb`
- Make figures for subsampled UKB (rv1): `job.simulation/make_figure.UKB_subsample.rv1.ipynb`


## Cell type-level results (Fig. 3): `job.celltype_association`
- Summary of the cell-type association results: `job.celltype_association/summary_ct.ipynb`
- Main analysis: `job.celltype_association/main_figure.ipynb`
- Main analysis (rv1): `job.celltype_association/main_figure.rv1.ipynb`
- Comparison of cell-type association for three atlas datasets: TMS FACS, TMS droplet, TS FACS: `job.celltype_association/atlas_compare.ipynb`
- Relationship between scDRS power and heritability, polygenicity: `job.celltype_association/optim_param.ipynb`
- Comparison of cell-type association to LDSC-SEG: `job.celltype_association/ldsc_compare.ipynb`
- Comparison of cell-type association to alternative methods (rv1): `job.celltype_association/methods_compare.rv1.ipynb`
- Effects of gene sets for scDRS power: `job.celltype_association/vary_geneset.ipynb`
- Evaluation of alternative versions of scDRS using control traits and cell types (rv1): `job.continuous_score/` (see the directory for more details)

## T cell example (Fig. 4): `job.case_tcell`
- Reprocess TMS T cells and assign effectorness gradients: `job.case_tcell/s1_reprocess_tms_tcell.ipynb`
- Main analysis: `job.case_tcell/s3_analysis_tcell.ipynb`
- Main analysis (rv1): `job.case_tcell/s3_analysis_tcell.rv1.ipynb`
- Replication in Cano-Gamez & Soskic et al. and Nathan et al. data: `job.case_tcell/s4_analysis_tcell.replication.ipynb`
- Replication in Cano-Gamez & Soskic et al. and Nathan et al. data (rv1): `job.case_tcell/s4_analysis_tcell.replication.rv1.ipynb`
- Cluster-level LDSC-SEG analysis: `job.case_tcell/s5_compare_ldsc_cluster_4res.ipynb`
- Cluster-level LDSC-SEG analysis (rv1): `job.case_tcell/s5_compare_ldsc_cluster_4res.rv1.ipynb`
- Disease gene prioritization: `job.case_tcell/s6_gene_prioritization.ipynb`
- Disease gene prioritization (rv1): `job.case_tcell/s6_gene_prioritization.rv1.ipynb`
- Gene set overlap and disease score correlation between traits (rv1): `job.case_tcell/s7_contrast_traits.rv1.ipynb`
- Automatic annotation using ProjecTILE in R (rv1): `job.case_tcell/s8_map_tms_tcell_ProjecTILs.rv1.ipynb`

## Neuron example (Fig. 5AB):  `job.ca1_pyramidal`
- Main analysis (Fig. 5AB): `job.ca1_pyramidal/main_figure.ipynb`
- Main analysis (Fig. 5AB) (rv1): `job.ca1_pyramidal/main_figure.rv1.ipynb`
- Analysis of neurons in TMS FACS dataset: `job.ca1_pyramidal/tms.ipynb` 
- Analysis of Zeisel et al. 2015 dataset: `job.ca1_pyramidal/zeisel.ipynb`
- Verification of the inferred spatial coordinates: `job.ca1_pyramidal/spatial_verify.ipynb`

## Hepatocyte example (Fig. 5CD): `job.case_hepatocyte`
- Reprocess TMS hepatocytes: `job.case_hepatocyte/s1_reprocess_tms_hep.ipynb`
- Main analysis: `job.case_hepatocyte/s3_analysis_hep.ipynb`
- Main analysis (rv1): `job.case_hepatocyte/s3_analysis_hep.rv1.ipynb`
