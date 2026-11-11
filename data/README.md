# Data

Processed data used in the paper. Raw yearbook files are not redistributed; all series below are the processed versions.

## Main panel

`panel_AI5.csv` — balanced panel of 30 provincial-level regions, 2011–2022 (N = 360). Key columns:

| Column | Meaning |
|---|---|
| `sh` | Province (short Chinese name) |
| `年份` | Year |
| `DEI` | Digital economy index (13 indicators, combined entropy + CRITIC–Spearman weights) |
| `CEI` | Carbon emission pressure index (4 indicators, same combined weighting) |
| `CE` | Carbon emission intensity |
| `lnCE`, `lnDEI`, `lnCEI` | Natural logs |
| `lnpgdp` | Log per capita GDP |
| `fisc` | Fiscal expenditure / GDP |
| `dens` | Population density (people/km²) |
| `fin` | (Deposits + loans) / GDP |
| `urb` | Urbanization rate |
| `mkt` | Marketization index |
| `inds` | Industrial structure upgrading (tertiary / secondary) |
| `ee` | Energy efficiency (output per unit of energy) |

`panel_with_pca.csv` adds PCA-based DEI and CEI measures used in the measurement-robustness table.

## Spatial weights and Moran output

| File | Content |
|---|---|
| `W_economic_distance.csv` | Pure economic-distance weight matrix (row-standardized) |
| `W_geographic_contiguity.csv` | Geographic contiguity weight matrix (row-standardized) |
| `Wadj.txt`, `Wadj_bin.txt`, `Wdist.txt`, `Wecon.txt` | Weight matrices in Stata `spatwmat` text format |
| `FWadj.dta`, `FWdist.dta`, `FWecon.dta` | Weight matrices in Stata `.dta` format |
| `panel.txt` | Panel in space-delimited Stata input format for the Moran do-file |
| `local_moran.txt` | Local Moran output, 2011–2022 |
| `moran_width.txt` | Global Moran's I, DEI/CEI × three weights (paper Table 8) |
| `lisa_quadrants.csv` | Province-year LISA quadrant used in Figure 5 |

## Sources

China Statistical Yearbook; China Energy Statistical Yearbook; China Internet Network Information Center; Ministry of Industry and Information Technology; China Industrial Statistical Yearbook; Peking University Digital Financial Inclusion Index. Carbon emissions are calculated from 17 fossil fuels using the IPCC sectoral approach and are consistent with the CEADs accounting framework.

License: CC BY 4.0.
