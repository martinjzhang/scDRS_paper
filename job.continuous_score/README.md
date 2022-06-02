# Experiments to explore the utility of using continuous score

- `01-compile-geneset.ipynb`: compile the gene set files for various configurations.
- `submit.compute_subsample_score.sh` / `submit.compute_full_score.sh`: calculate scDRS score for each of the scDRS gene set version.
- `02-power-analysis.ipynb`: calculate the the t-stats comparing positive control cell types and negative control cell types for each of the comiled trait.
- `03-plot-power.ipynb`: produce the summary of the results.
- `03-plot-power.rv_final.ipynb`: produce the summary of the results (plotting styles updated).
- `04-config-robustness.ipynb`: assess the robustness of the results across different scDRS versions.
- `04-config-robustness.rv_final.ipynb`: assess the robustness of the results across different scDRS versions (plotting styles updated)