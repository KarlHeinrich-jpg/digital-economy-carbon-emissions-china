# Code

Two stacks are archived here.

## Python — table replication

`python/replicate_tables.py` reproduces the baseline and robustness estimates from the processed panel:

```bash
pip install pandas numpy linearmodels
python code/python/replicate_tables.py
```

It prints the baseline two-way fixed-effects estimates with 0/6 controls, the carbon-intensity outcome, the industrial-structure variants, lagged DEI, exclusion of the four municipalities, 1% winsorization, and descriptive statistics. The script reads `data/panel_AI5.csv` by a path relative to its own location, so it can be run from any directory.

## Stata — archived do-files

`stata/` contains the do-files that produced the paper's tables in Stata 17/18. They were archived from the working project and adapted: the working directory is now the repository's `data/` folder, and the panel file is `panel_AI5.csv`.

Before running, either:

```stata
cd "path/to/digital-economy-carbon-emissions-china/data"
```

or edit the `cd` line at the top of the do-file.

| File | Paper table / purpose |
|---|---|
| `01_baseline_panel_fe.do` | Baseline two-way FE with 0/6 controls (Table 4) |
| `02_mechanism_energy_efficiency.do` | Energy-efficiency channel (Table 7) |
| `03_spatial_sdm_sar_sem.do` | Spatial models under economic distance |
| `04_composite_weight_sdm.do` | Spatial models under the composite weight |
| `05_system_gmm.do` | System GMM (Table 12) |
| `06_global_moran.do` | Global Moran's I and Geary's c (Table 8) |

The spatial do-files use `xsmle` and `spatwmat`; install them with `ssc install xsmle` and `ssc install spatwmat` if needed. The do-files retain their original Chinese comments.

## Figures

The six figures are archived as PDF and PNG in `manuscript/figures/`. Figure 1 is a TikZ source (`fig01_framework.tex`, compile with `xelatex`). The maps and LISA panels were produced from the panel and weight files in `data/` with the authors' Python mapping workflow; the final vector files are included so the published figures can be reproduced exactly without the auxiliary geographic assets.
