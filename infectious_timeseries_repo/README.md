# Data Repository for "Tools for Analyzing Infectious Disease Sequence Data and Forecasting Public Health Outcomes"


### **Basic Information**


This repository contains a collection infectious disease time series data obtained from a variety of publicly-available sources compiled as part of the Directed Research project entitled, "Merging Sequence Data with Epidemiological Forecasting to Predict Future Pandemics." This project began in Fiscal Year 2024 at Los Alamos National Lab under the Laboratory Directed Research and Development (LDRD) program. The project leads are:

- Lauren Castro (PI) lcastro@lanl.gov;
- Dave Osthus (co-PI) dosthus@lanl.gov;
- Will Fischer (co-PI) wfischer@lanl.gov.

The organized and cleaned data are provided in **Organized_Lists.zip**. These lists are structured such that each element of the list contains an infectious disease time series for a particular location, disease, data source, reporting cadence (e.g., daily vs weekly) and endpoint. For more information about the structure of these lists, see **Summarize_Data.pptx**. 

Within each data source folder, we provide information about how we accessed the data and code we used to clean the data from each source. To facilitate modeling, observed weekly data were assigned to the nearest regular week or recast to follow an exactly weekly cycle (with zero-count weeks added as appropriate). This did not always align exactly with the actual reporting dates for the individual data streams. In general, see each data source’s **README** for information about data cleaning choices made and how the original reporting dates were recorded. Users wanting a more careful look at the reporting dates in the original data (for example, to harmonize these data with other data streams or to do detailed epidemiological studies) should dig into the data cleaning code or return to the original data source rather than use our processed data. The goal of this repository is to provide a database of infectious disease time series data to use for training and validating predictive models, and our data processing choices reflect this. 

### **How to Cite this Repository**


All data in this repository are publicly available, but individual data sources may have additional limitations regarding citation, acceptable use cases/licensing, etc. We have done our best to summarize the information provided for each data stream in the individual data stream **READMEs**. Users wanting to use these data for commercial purposes should consult the links provided for the original data to verify acceptable use cases. Users of this repository should also cite this GitHub repository. 

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
