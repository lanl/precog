---
output:
  pdf_document: default
  html_document: default
---
Code for Experiment in Publication “A method to quantify infectious disease variant proportion bias arising from reporting delay”
========

**reporting-delay** is a sub repository of the *Precog* project specific to quantifying reporting delay bias of SARS-CoV-2 variant proportions. Specifically, this repository contains all code used to perform the experiments in the cited publication.  For questions, issues, or clarifications please reach out to Marina Mancuso: <mmancuso@lanl.gov>.

---
## Cloning the repository

The following commands will permit you to downlaod the Reporting Delay sub repository from the Precog project:

```bash

git clone --filter=blob:none --sparse https://github.com/lanl/precog.git
cd precog
git sparse-checkout set reporting-delay

```
---

## 🔬 Data

The data used in the experiments are located on [Zenodo](https://doi.org/10.5281/zenodo.20128960). The file should be named `reporting_delay_data.csv` and saved to this directory. This dataset contains dates of sequence collection, location of sequence collection, sequence pangolineage, number of days between sample collection and sequence submission, and number of sequences. The SARS-CoV-2 sequences are publically accessible on [GISAID](https://weekly.chinacdc.cn/en/article/doi/10.46234/ccdcw2021.255) and part of [EPI_SET_250929ge](https://doi.org/10.55876/gis8.250929ge).

---

## 🚀 Code

R scripts are numbered by workflow order. 

### 0_calculate_bias_metrics.R

This script calculates Cramer's V and vector norm metrics for quantifying evidence of SARS-CoV-2 variant composition bias. It compares similarity between validation (final) variant proportions with reported proportions at 7, 14, 21, and 30 days post-collection. Results are saved to `all_cramer_results.RData` and `all_cramer_results.csv` files.

### 1_simulations_under_null.R

This script determines statistical significance of bias quantification metrics. It generates 1,000 sample simulations under the null hypothesis (validation composition) and compares them to the observed metric values from `all_cramer_results.RData` in terms of p-values and 95% percentiles. Results are saved to `all_sim_results.RData`.

### 2_find_emerging_variants.R

This script finds validation and near-real-time variant proportions for each location, collection date, and pangolineage at 7, 14, 21, and 30 days post-collection. Results are saved to `emerging_variants.csv`.

### 3_plots_for_pub.R

This script generates visualizations in the cited publication.

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
