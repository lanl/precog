# US HHS Data

### **Download Link**


https://github.com/reichlab/flusion/blob/main/data-raw/influenza-hhs/download-hhs.R

### **Citation**


@misc{ray2024flusionintegratingmultipledata,
      title={Flusion: Integrating multiple data sources for accurate influenza predictions}, 
      author={Evan L. Ray and Yijin Wang and Russell D. Wolfinger and Nicholas G. Reich},
      year={2024},
      eprint={2407.19054},
      archivePrefix={arXiv},
      primaryClass={stat.ML},
      url={https://arxiv.org/abs/2407.19054}, 
}

### **Licensing** 


MIT License


Copyright (c) 2023 The Reich Lab at UMass-Amherst


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:


The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.


THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

### **Data Cleaning**


These data were obtained via the Reich Lab GitHub. Daily reported hospitalization dates were modified by 1 following precedent in the code provided by the Reich Lab. Daily reports were aggregated to weekly reports also following code provided by the Reich Lab. If a week had no reported influenza hospitalizations for a given location, we report 0 hospitalizations, ensuring that the resulting time series has regular weekly reporting without gaps. We require there to be at least 10 observations per time series for inclusion. 

### **Date of Raw Data Access:**


2/23/2026
