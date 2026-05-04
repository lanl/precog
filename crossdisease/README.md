
# Code for Publication "Transfer Learning using 66 Diseases for Disease Forecasting Applications"

### **Basic Information**

This repository contains code associated the manuscript cited below. Raw data used in this analysis has been provided in https://github.com/lanl/precog/tree/main/infectious_timeseries_repo and not duplicated here. This work was conducted as part of the Directed Research project entitled, "Merging Sequence Data with Epidemiological Forecasting to Predict Future Pandemics." This project began in Fiscal Year 2024 at Los Alamos National Lab under the Laboratory Directed Research and Development (LDRD) program. 

After downloading the raw data, users can replicate our analyses by running code in the following folders in order:

- features (get_lstm_embeddings.R, get_moa_embeddings.R, make_gbm_real.R, evaluate_forecastability.R)
- models (train_gbm.R, train_lstm.R)
- evaluations (evaluate_gbm.R, evaluate_lstm.R, evaluate_moa.R)
- uq (uq.R)
- summaries (get_summary_files.R, get_time_series_features.R, get_summary_files_moacovid.R, summarize_all.R, summarize_moa_covid.R, viz_raw_data.R)
- summaries_uq (get_summary_files_uq.R, summarize_all.T)



### **Citation**

**Lauren J Beesley, Alexander C Murph, Dave Osthus, and Lauren A Castro. "Transfer Learning using 66 Diseases for Disease Forecasting Applications." In preparation.**

### **How to Cite this Repository**

When citing this repository, please cite both the corresponding manuscript and the corresponding compiled data in https://github.com/lanl/precog/tree/main/infectious_timeseries_repo. 

**Lauren J Beesley, Alexander C Murph, Dave Osthus, and Lauren A Castro. "Transfer Learning using 66 Diseases for Disease Forecasting Applications." In preparation.**


### **Software Release**


This software has been approved for open source release and has been assigned O4726.

### **Copyright**


© 2024. Triad National Security, LLC. All rights reserved. This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

### **License**


See individual data stream READMEs for dataset-specific licenses. 


This code repository is distributed under the BSD-3 License:
Copyright 2024. Triad National Security, LLC.


Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1.	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
	
2.	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
	
3.	Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


### **Contact**

This repository was compiled by researchers at Los Alamos National Laboratory to facilitate the development of time series forecast models. Questions may be directed to **Lauren J Beesley** at lvandervort@lanl.gov. 
