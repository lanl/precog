# Johns Hopkins CSSEGIS and Our World in Data COVID-19 Data

### **Download Link**


https://github.com/CSSEGISandData/COVID-19_Unified-Dataset
https://github.com/owid/covid-19-data/tree/master/public/data 

### **Citation**


Badr, H.S., Zaitchik, B.F., Kerr, G.H. et al. Unified real-time environmental-epidemiological data for multiscale modeling of the COVID-19 pandemic. Sci Data 10, 367 (2023). https://doi.org/10.1038/s41597-023-02276-y

Edouard Mathieu, Hannah Ritchie, Lucas Rod\u00e9s-Guirao, Cameron Appel, Charlie Giattino, Joe Hasell, Bobbie Macdonald, Saloni Dattani, Diana Beltekian, Esteban Ortiz-Ospina and Max Roser (2020) - "Coronavirus Pandemic (COVID-19)". Published online at OurWorldInData.org. Retrieved from: 'https://ourworldindata.org/coronavirus'


### **Licensing (CSSEGIS)**


Copyright (c) 2023 Johns Hopkins University (JHU)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

### **Licensing (OWID)**


CC BY 4.0
All visualizations, data, and code produced by Our World in Data are completely open access under the Creative Commons BY license. You have the permission to use, distribute, and reproduce these in any medium, provided the source and authors are credited.
The data produced by third parties and made available by Our World in Data is subject to the license terms from the original third-party authors. We will always indicate the original source of the data in our database, and you should always check the license of any such third-party data before use.

### **Data Cleaning**


Johns Hopkins University (JHU) CSSEGIS COVID-19 data were subsetted to Admin0 (i.e., country) and Admin1 (i.e., region/state within country) observations. Several locations were excluded due to a comparatively small number of reported cases (e.g., Antarctica). The JHU COVID-19 dataset includes case and death reports for many different sources, often with duplicate and conflicting entries for each date and location. For each date and location, we included the cases/deaths counts recorded from the perceived highest quality source in the following order of prioritization: JHU, JRC, NYT, SES, RKI, DPC, CTP, NYC. See the Johns Hopkins CSSEGIS repo documentation (linked above) for more details. Additional data were available from JHU that were not included in the output data, e.g., vaccinations, Government policy, etc. 
We also sourced hospitalization data from Our World in Health (OWID). Locations were similarly subsampled to exclude locations with a comparatively small number of reported cases (e.g., the Vatican). This dataset also contained many more variables than were utilized in this repository. Users interested in digging into the full source dataset can find these data in the ./raw_data directory. 
JHU cases/deaths data were available between 1/2/2020 and 3/31/2023. OWID hospitalization data was available between 2/12/2020 and 9/20/2023. For a given outcome (cases, deaths, hospitalizations) and location, we required at least 10 weeks of data to be available for weekly data and at least 10 days of data to be available for daily data. Negative cases, deaths, hospitalizations reports were replaced with zeros.  

### **Date of Raw Data Access:**


9/20/2023 (JHU)

9/21/2023 (OWID)
