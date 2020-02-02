# SpykeMapper
This is a Python module to analyze voltage imaging data in a variety of ways, including but not limited to subthreshold waveform extraction, spike-triggered averaging, pairwise cross-correlation, and factor analysis (courtesy of the FactorAnalyzer package). It includes a class for importing and processing voltage imaging data, which is assumed to be in the form of an Excel workbook where the first column is fime and the remaining columns are raw data from individual cells. 

# Description
Voltage imaging enables collection of "intracellular" electrophysiology data from several cells at once. Separating this data into subthreshold waveforms and spike trains allows a researcher to infer several properties of the measured neurons that are inaccessible or challenging to study with classical, electrode-based techniques. 

# To do:

Import TIF and draw ROIs rather than Excel files.
  
* Using imported TIF, use first frame OR supplied brightfield image for drawing ROIs
* Save ROIs, load them later
* Compact data from ROIs into "average value" arrays - export this raw data if user wants

* Ideas?

# Requirements:
* Python 3.4 or higher
* numpy
* pandas
* scipy
* scikit-learn
* FactorAnalyzer
* easygui (will try to remove this dependency at some point)
