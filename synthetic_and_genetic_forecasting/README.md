# README for “Leveraging Synthetic and Genetic Data to Improve Epidemic Forecasting”

**Osthus et al.**

## Overview

This repository contains a series of R scripts in the `code/` directory used to run the full analysis pipeline.

Script filenames follow the convention `(#)_...`, and the scripts should be run in numeric order.

Below is a description of what each script does, along with its inputs and outputs.

---

## `(0)_make_train_and_test_data.R`

### Description
Formats the training and testing data.

### Inputs
- All files in `~/data_raw/mutantigen/`
- `~/data_raw/real_respiratory_data_complete.RDS`
- `~/data_raw/dfreal_var_attr.RDS`

### Outputs
- `~/data_clean/train_syn_tc.RDS`
- `~/data_clean/train_syn_vac.RDS`
- `~/data_clean/train_real.RDS`
- `~/data_clean/test_covid.RDS`

---

## `(1)_train_model.R`

### Description
Trains the transformer models for `M(r,.)`, `M(st,.)`, `M(sv,.)`, and `M(a,.)`. This script takes many hours to run.

### Inputs
- `~/data_clean/train_syn_tc.RDS`
- `~/data_clean/train_syn_vac.RDS`
- `~/data_clean/train_real.RDS`

### Outputs
- All files in `~/trained_models/`

---

## `(2)_forecast_model.R`

### Description
Uses the trained models from `(1)_train_model.R` to generate forecasts on the COVID test set.

### Inputs
- `~/data_clean/test_covid.RDS`
- All files in `~/trained_models/`

### Outputs
- `~/output/forecasts.csv`

---

## `(3)_make_results_and_figs.R`

### Description
Plots most of the figures in the paper.

### Inputs
- `~/data_clean/test_covid.RDS`
- `~/output/forecasts.csv`
- All config files in `~/output/`

### Outputs
- Most of the EPS files in `~/figs/`

---

## `(4)_perform_by_ts_trends.R`

### Description
Performs the analysis needed to generate Figure 16 in the paper.

### Inputs
- `~/data_clean/test_covid.RDS`
- `~/output/forecasts.csv`
- All files in `~/data_raw/ts_phases/`

### Outputs
- `~/figs/phase_of_outbreak.eps`

---

## `(5)_umap_analysis.R`

### Description
Performs UMAP analysis and builds a classifier used to generate Figure 13.

### Inputs
- `~/trained_models/cfg_all.rds`
- `~/data_clean/train_real.RDS`
- `~/data_clean/train_syn_tc.RDS`
- `~/data_clean/test_covid.RDS`

### Outputs
- `~/figs/prob_syn_tc_by_data_types.eps`
- `~/figs/prob_syn_tc_allstates.eps`

---

## `(6)_show_mutantigen_stochasticity.R`

### Description
Performs the analysis needed to generate Figure S1.

### Inputs
- Selected files in `~/data_raw/mutantigen/out_XX.timeseries`

### Outputs
- `~/figs/mutantigen_stochasticity.eps`

---

## `(7a)_train_models_with_training_subsets.R`

### Description
Fits all transformer models to subsets of the training data to support the analysis presented in Figures 11 and 12. This script takes many hours to run.

### Inputs
- `~/data_clean/train_syn_tc.RDS`
- `~/data_clean/train_syn_vac.RDS`
- `~/data_clean/train_real.RDS`

### Outputs
- All files in `~/trained_models_subsets/`

---

## `(7b)_forecast_model_with_training_subsets.R`

### Description
Performs forecasts using all models created in `(7a)_train_models_with_training_subsets.R` to support the analysis presented in Figures 11 and 12. This script takes over an hour to run.

### Inputs
- `~/data_clean/test_covid.RDS`
- All files in `~/trained_models_subsets/`

### Outputs
- All CSV files beginning with `~/output/forecasts_model_subsets_mod_`

---

## `(7c)_make_results_and_figs_training_subsets.R`

### Description
Performs the analysis needed to generate Figures 11 and 12.

### Inputs
- `~/data_clean/test_covid.RDS`
- All CSV files beginning with `~/output/forecasts_model_subsets_mod_`

### Outputs
- `~/figs/training_subsets_paired_differences.eps`
- `~/figs/training_subsets_metrics.eps`


Software Release

This software has been approved for open source release and has been assigned O4726.

Copyright

© 2024. Triad National Security, LLC. All rights reserved. This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

License

See individual data stream READMEs for dataset-specific licenses.

This code repository is distributed under the BSD-3 License: Copyright 2024. Triad National Security, LLC.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

