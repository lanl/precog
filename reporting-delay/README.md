# Code for “A method to quantify infectious disease variant proportion bias arising from reporting delay”

**reporting-delay** is a subdirectory of the *Precog* project containing the code used to quantify evidence of bias in near-real-time SARS-CoV-2 variant proportions arising from sequence-reporting delay. For questions, issues, or clarifications, please contact Marina Mancuso at <mmancuso@lanl.gov>.

---
## Cloning the repository

The following commands download only the `reporting-delay` directory from the Precog repository:

```bash
git clone --filter=blob:none --sparse https://github.com/lanl/precog.git
cd precog
git sparse-checkout set reporting-delay
```

## Data

The data used in the analyses are available on [Zenodo](https://doi.org/10.5281/zenodo.20128960). Download the file as `reporting_delay_data.csv` and place it directly in `reporting-delay/`. Each row gives a sequence collection date, collection location (`Admin0`), Pango lineage (`pango`), reporting delay in days (`delay_days`), and sequence count (`counts`). The underlying SARS-CoV-2 sequences are publicly accessible through [GISAID](https://weekly.chinacdc.cn/en/article/doi/10.46234/ccdcw2021.255) as part of [EPI_SET_250929ge](https://doi.org/10.55876/gis8.250929ge).

The scripts expect the following directories to exist:

```text
reporting-delay/
├── reporting_delay_data.csv
├── results/
├── plots/
│   ├── main_body/
│   └── supplement/
└── *.R
```

Paths in the scripts are defined relative to the parent `precog/` directory, so run the scripts from that directory. Scripts are numbered in workflow order.

## Analysis workflow

### `0_calculate_bias_metrics.R`

**Input:** `reporting-delay/reporting_delay_data.csv`

For each location, collection date, and 7-, 14-, 21-, or 30-day reporting threshold, this script compares near-real-time variant proportions with validation proportions calculated after all sequences have been reported. Calculations use a seven-day collection window by default and include Cohen's $w$, Fei, the L1, L2, and L-infinity vector norms, near-real-time sample size (`n`), reporting rate (`r`), and number of co-circulating lineages (`K`). It also performs sensitivity analyses using three- and 14-day collection windows and broader, collapsed Omicron categories.

**Outputs:**

- `results/all_metric_results.csv`: metrics from the primary seven-day collection window.
- `results/all_metric_results.RData`: filtered primary results, including `p_nrt` and `p_val`, used as input to the null simulations.
- `results/all_metric_results_3_window.csv`: three-day collection-window sensitivity analysis for Brazil, Denmark, and the United States.
- `results/all_metric_results_14_window.csv`: 14-day collection-window sensitivity analysis for Brazil, Denmark, and the United States.
- `results/all_metric_results_collapse_omicron.csv`: seven-day collection-window sensitivity analysis using collapsed Omicron categories for Brazil, Denmark, and the United States.

### `1_simulations_under_null.R`

**Input:** `results/all_metric_results.RData`

For every eligible location-date-delay observation, this script draws multinomial samples of size `n` from the corresponding validation composition (`p_val`). It calculates null distributions, Monte Carlo p-values, and 95th percentiles for Cohen's $w$, Fei, and the three vector norms.

**Outputs:**

- `results/all_sim_results.RData`: primary results from 1,000 draws per observation.
- `results/sim_results_M_500.RData`: sensitivity results from 500 draws per observation.
- `results/sim_results_M_5000.RData`: sensitivity results from 5,000 draws per observation.

Each file contains the `metric_comb` object, which combines observed metrics with their simulated p-values and 95th percentiles and adds sample-size, reporting-rate, and p-value categories.

### `2_find_emerging_variants.R`

**Input:** `reporting-delay/reporting_delay_data.csv`

For each location, collection date, lineage, and 7-, 14-, 21-, or 30-day reporting threshold, this script calculates the validation proportion (`p_val`), near-real-time proportion (`p_nrt`), and their difference (`p_val - p_nrt`) using a seven-day collection window.

**Output:** `results/emerging_variant_ts.csv`

### `3_dominant_variant_results.R`

**Inputs:**

- `results/emerging_variant_ts.csv`
- `results/all_metric_results.csv`

This script identifies the dominant validation lineage for each location and date, resolves ties sequentially, calculates its daily rate of change, and flags windows surrounding changes in the dominant lineage. It then joins these results to the bias metrics.

**Outputs:**

- `results/dominant_variant_rate_of_change.csv`: dominant lineage, validation proportion, daily rate of change, and three turnover-window indicators by location and date.
- `results/dominant_variant_w_merged.csv`: dominant-variant results joined to the bias metrics for November 1, 2020 through December 31, 2022.

### `4_main_body_figures.R`

**Inputs:**

- `reporting-delay/reporting_delay_data.csv`
- `results/all_metric_results.csv`
- `results/all_sim_results.RData`
- `results/emerging_variant_ts.csv`
- `results/dominant_variant_w_merged.csv`

This script generates the manuscript's main-text figures and prints several supporting numerical summaries to the R console.

**Outputs in `plots/main_body/`:**

- `compare_delays.png`
- `n_r_cats_all_countries.png`
- `w_vs_reporting_rate.png`
- `Compare_w_Locations.png`
- `US_14_emerge.png`
- `Denmark_14_emerge.png`
- `dominant_variant.png`
- `variant_turnover.png`

### `5_supplement_figures.R`

**Inputs:**

- `reporting-delay/reporting_delay_data.csv`
- `results/all_metric_results.csv`
- `results/all_metric_results_collapse_omicron.csv`
- `results/all_sim_results.RData`
- `results/emerging_variant_ts.csv`

This script generates supplemental analyses and figures, including comparisons among bias metrics, null-distribution thresholds, emerging-variant time series, and the Omicron-category sensitivity analysis. It also prints metric correlations and full-versus-collapsed Omicron correlations to the R console.

**Outputs in `plots/supplement/`:**

- `n_r_cats.png`
- `w_vs_Norms_box.png`
- `L1_vs_reporting_rate.png`
- `q95.png`
- `US_7_emerge.png`
- `Denmark_7_emerge.png`
- `Brazil_7_emerge.png`
- `Brazil_14_emerge.png`
- `category_sensitivity.png`

## Running the analysis

From the parent `precog/` directory, run the scripts in order:

```bash
Rscript reporting-delay/0_calculate_bias_metrics.R
Rscript reporting-delay/1_simulations_under_null.R
Rscript reporting-delay/2_find_emerging_variants.R
Rscript reporting-delay/3_dominant_variant_results.R
Rscript reporting-delay/4_main_body_figures.R
Rscript reporting-delay/5_supplement_figures.R
```

---

## Citation
M. Mancuso, L.J. Beesley, D. Osthus, L.A. Castro. (202x). A method to quantify infectious disease variant proportion bias arising from reporting delay.  _In Review._ 

## Release

This software has been approved for open source release and has been assigned **O4726** 

## Copyright

© 2024. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

## License

This code repository is distributed under the BSD-3 License:

Copyright 2024. Triad National Security, LLC.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
