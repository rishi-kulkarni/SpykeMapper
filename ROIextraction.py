# -*- coding: utf-8 -*-
"""
Functions to Import a TIFF Stack and Select an ROI
"""
import matplotlib
import matplotlib.pyplot as plt
plt.ioff() 
from roipoly import roipoly
import pylab as pl
from skimage import io
import numpy as np
import pandas as pd

plt.ioff() 
    
   
image = io.imread("isoTMR.tiff")
plt.imshow(image[1,:,:])

plt.title("Foreground: left click: line segment         right click: close region")
fgROI = roipoly(roicolor='r')
plt.show(block=True)
plt.ioff() 

plt.imshow(image[1,:,:], interpolation='nearest', cmap="Greys")
plt.title("Background: left click: line segment         right click: close region")
plt.colorbar()
fgROI.displayROI()
bgROI = roipoly(roicolor='b')

fgmask = fgROI.getMask(image[1,:,:])
bgmask = bgROI.getMask(image[1,:,:])
output = np.zeros((len(image[:,1,1]),3))
timepoints = np.arange(len(image[:,1,1]))
output[:,0] = timepoints
for key in range(len(output[:,1])):
    output[key][1] = pl.mean(image[key,fgmask])
    output[key][2] = pl.mean(image[key,bgmask])
 
forexcel = pd.DataFrame(output[:,[1,2]], index = output[:,0], columns=["cell","bg"])
forexcel.index.name = "time"
forexcel.to_clipboard(excel=True)