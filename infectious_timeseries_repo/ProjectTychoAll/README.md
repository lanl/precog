# Project Tycho 

### **Download Link**


https://www.tycho.pitt.edu/data/

### **Citation**


Willem G van Panhuis, Anne Cross, Donald S Burke, Project Tycho 2.0: a repository to improve the integration and reuse of data for global population health, Journal of the American Medical Informatics Association, Volume 25, Issue 12, December 2018, Pages 1608–1617, https://doi.org/10.1093/jamia/ocy123

### **Licensing**


Creative Commons Attribution 4.0 International Public License.

Subject to the terms and conditions of this Public License, the Licensor hereby grants You a worldwide, royalty-free, non-sublicensable, non-exclusive, irrevocable license to exercise the Licensed Rights in the Licensed Material to:

	a	reproduce and Share the Licensed Material, in whole or in part; and
  
	b	produce, reproduce, and Share Adapted Material.

Your exercise of the Licensed Rights is expressly made subject to the following conditions.

	1	If You Share the Licensed Material (including in modified form), You must:
	  a	retain the following if it is supplied by the Licensor with the Licensed Material:
	    i	identification of the creator(s) of the Licensed Material and any others designated to receive attribution, in any reasonable manner requested by the Licensor (including by pseudonym if designated);
	    ii	a copyright notice;
	    iii	a notice that refers to this Public License;
	    iv	a notice that refers to the disclaimer of warranties;
	    v	a URI or hyperlink to the Licensed Material to the extent reasonably practicable;
	  b	indicate if You modified the Licensed Material and retain an indication of any previous modifications; and
	  c	indicate the Licensed Material is licensed under this Public License, and include the text of, or the URI or hyperlink to, this Public License.
  
	2	You may satisfy the conditions in Section 3(a)(1) in any reasonable manner based on the medium, means, and context in which You Share the Licensed Material. For example, it may be reasonable to satisfy the conditions by providing a URI or hyperlink to a resource that includes the required information.
  
	3	If requested by the Licensor, You must remove any of the information required by Section 3(a)(1)(A) to the extent reasonably practicable.
  
	4	If You Share Adapted Material You produce, the Adapter's License You apply must not prevent recipients of the Adapted Material from complying with this Public License.

Section 5 – Disclaimer of Warranties and Limitation of Liability.

	a	Unless otherwise separately undertaken by the Licensor, to the extent possible, the Licensor offers the Licensed Material as-is and as-available, and makes no representations or warranties of any kind concerning the Licensed Material, whether express, implied, statutory, or other. This includes, without limitation, warranties of title, merchantability, fitness for a particular purpose, non-infringement, absence of latent or other defects, accuracy, or the presence or absence of errors, whether or not known or discoverable. Where disclaimers of warranties are not allowed in full or in part, this disclaimer may not apply to You.
  
	b	To the extent possible, in no event will the Licensor be liable to You on any legal theory (including, without limitation, negligence) or otherwise for any direct, special, indirect, incidental, consequential, punitive, exemplary, or other losses, costs, expenses, or damages arising out of this Public License or use of the Licensed Material, even if the Licensor has been advised of the possibility of such losses, costs, expenses, or damages. Where a limitation of liability is not allowed in full or in part, this limitation may not apply to You.
  
	c	The disclaimer of warranties and limitation of liability provided above shall be interpreted in a manner that, to the extent possible, most closely approximates an absolute disclaimer and waiver of all liability.

### **Data Cleaning**


Data were obtained using the Project Tycho API. Substantial data cleaning was implemented for one diseases. We will summarize these choices here, and we refer uses to the code used for data cleaning for details. For a given disease, locations were only included if they had at least 50 dates for which nonzero case counts were reported. 

For all diseases considered here, cases were reported cumulatively, and these cumulative case counts were occasionally restarted at zero in the middle of a given time series. These cumulative time series were converted to weekly incident cases. Occasionally, the transformed incident case reports would span multiple weeks (i.e., correspond to cases reported for several weeks combined). These cases were naively assigned to the first week in the corresponding time interval, and subsequent weeks were assigned zero cases. Users may also consider implementing their own code to instead allocate these cases evenly between the corresponding weeks. All weekly time series were mapped to a regular weekly reporting interval by determining the regular reporting week nearest in time to the reported start date in Project Tycho. This choice was made to facilitate generic forecast model development, and we recorded the original interval start and end dates in the _ts_dates_start_actual_ and _ts_dates_end_actual_ fields within the time series data list. The canonical reporting cadence (weekly/daily) for a given time series was determined as the median difference between the identified original reporting interval start and end dates. 

We emphasize that the released, cleaned time series reflect processed version of the original Project Tycho data in terms of the original reporting dates. Users interested in a more accurate-in-calendar-time time series for epidemiological research or for merging with other time-varying data streams may be better served by using our data cleaning/processing code as a starting place for developing their own data processing pipeline. 

### **Date of Raw Data Access:** 


11/19/2024
 

