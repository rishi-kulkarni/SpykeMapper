import ROIextraction
import voltagestepfunctions
import pylab as pl
import pandas as pd

targetTiff = ROIextraction.openTiff()
TiffArray = ROIextraction.roi_to_array(targetTiff)
backgroundsub = voltagestepfunctions.backgroundSubtract(TiffArray)
bleachcorrection = voltagestepfunctions.baselineFlatten(backgroundsub)
deltaF = voltagestepfunctions.stepFinder(bleachcorrection)
pl.plot(deltaF[:,0],deltaF[:,1])
forexcel = pd.DataFrame(deltaF[:,1], index = deltaF[:,0], columns=["dF/F"])
forexcel.index.name = "mV"
forexcel.to_clipboard(excel=True)
print("Voltage Sensitivity Data copied to clipboard")
