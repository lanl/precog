# epiFFORMA: Data-Driven Ensemble Reweighting for Forecast Modeling
Code associated with epiFFORMA, a method for estimating time series forecast ensemble weights as a function of time series features. This repository provides code for implementing the method and for reproducing many of the figures in a corresponding manuscript. This documentation will be updated when that manuscript is publicly available. epiFFORMA builds on the FFORMA (Feature FORcasting Model Averaging; https://robjhyndman.com/papers/fforma.pdf) framework by catering decisions towards infectious disease forecasting applications.

---
## Summary

Accurate forecasting of infectious diseases drives modern public health interventions that reduce morbidity and mortality. However, accurate forecasting in real-time remains a challenge for the modeling community. Ensembling has emerged as a critical tool for accurate forecasting by leveraging multiple ``component" (individual) models into a single weighted average. Traditional ensembling strategies have relied on bespoke component models that weight the contributions of individual models according to extensive historical data for specific diseases. This is impractical for an emerging disease, since there would be very little -- if any -- data.   We propose an ensembling strategy, called epiFFORMA, that determines component weights for an ensemble model without historical data and is therefore disease-agnostic.
    
This strategy builds upon the FFORMA model from the M4 forecasting competition to harness epidemiological dynamics through synthetic data.  We demonstrate that epiFFORMA performs better than a naive, equal-weighting strategy when forecasting outbreaks of COVID-19, Diphtheria, ILI, Dengue, Measles, Mumps, Polio, Rubella, Smallpox, and Chikungunya.  We further show that epiFFORMA on average performs better than the individual component models in the ensemble.

## System requirements

The code is supported on all operating systems for which the requisite downloads (see below) are possible. The example code was tested on a MacBook Pro running macOS Sonoma 14.7.1, using R version 4.2.0.

## Installation

To downloading and install software and packages:
 - R (>= 4.1.2) follow instructions at https://www.r-project.org/

Installation should take less than 15 minutes on a normal desktop computer.

## Instructions for use

After R is installed, run **epifforma_run.R** to reproduce the epiFFORMA training and forecasting pipeline. No other scripts need to be run by the user; everything else is called by **epifforma_run.R**. Due to data storage limitations imposed by GitHub, we were unable to provide intermediate results (e.g., the evaluations of the real data) here. However, all of these results are reproducible with the provided data pipeline. 
  
## Release

These software have been approved for open source release and has been assigned **O4726**.

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
