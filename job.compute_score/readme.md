### Constants 
`dir=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/score_file`

### Score files for 74 GWAS traits 
- TMS FACS (`cov=const+n_genes+sex_male+age`): `dir/score.tms_facs_with_cov.magma_10kb_top1000_zscore`
- TMS FACS (`cov=const+n_genes+sex_male+age`; bin-gs): `dir/score.tms_facs_with_cov.magma_10kb_1000`
- TMS droplet (`cov=const+n_genes+sex_male+age`): `dir/score.tms_droplet_with_cov.magma_10kb_top1000_zscore`
- TMS droplet (`cov=const+n_genes+sex_male+age`; bin-gs): `dir/score.tms_droplet_with_cov.magma_10kb_1000`
- TS FACS (`cov=const+n_genes+donor`): `dir/score.ts_facs_with_cov.magma_10kb_top1000_zscore`
- TS FACS (`cov=const+n_genes+donor`; bin-gs): `dir/score.ts_facs_with_cov.magma_10kb_1000`
- Cano-Gamez & Soskic et al. (`cov=const+n_genes`): `dir/score.canogamez_with_cov.magma_10kb_top1000_zscore`
- Nathan et al. (`cov=const+n_genes+batch`): `dir/score.nathan_ni_2021_with_cov.magma_10kb_top1000_zscore`

### Other score files
- TMS FACS + MSigDB T cell signatures (`cov=const+n_genes+sex_male+age`): `dir/score.tms_facs_with_cov.tcell_sig`
- TMS FACS + Metabolic signatures (`cov=const+n_genes+sex_male+age`): `dir/score.tms_facs_with_cov.hep_metabolic`
