# SPykemapper
Tools for analyzing voltage imaging data.

By running main.py, SPykemapper will ask for a TIFF input, allow you to draw an ROI around the cell and background, and then copy to clipboard the deltaF/F values for the trace after background subtraction and bleach correction. 

Runs best from command line (navigate terminal to folder, type "python main.py" to start script").


# To do:
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
