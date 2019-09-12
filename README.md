# SpykeMapper
Tools for analyzing voltage imaging data.

This set of scripts takes an excel spreadsheet of raw voltage data as an input and outputs a variety of analyses, including separation of subthreshold and spiking behavior, pairwise cross-correlation between each trace, and a spike-triggered average using the subthreshold activity as the "input."


# To do:

Add Factor Analysis to the main script. 

Create a GUI for selecting what analyses you'd like to perform (performing all of them can be CPU intensive as the output files get rather large).

Import TIF and draw ROIs rather than Excel files.
  
* Using imported TIF, use first frame OR supplied brightfield image for drawing ROIs
* Save ROIs, load them later
* Compact data from ROIs into "average value" arrays - export this raw data if user wants

Batch Processing

User-inputted adjustments to the baseline fit (less important for neuronal activity)

* Four sliders (or two, if looking at neurons). 
* If user wants, apply correction from one trace to all traces if batch processing

Improved Spike Detection

* Ideas?
