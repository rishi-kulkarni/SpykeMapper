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

plt.ioff() 
    
   
image = io.imread("isoTMR.tiff")
plt.imshow(image[1,:,:])

plt.title("left click: line segment         right click: close region")
ROI = roipoly(roicolor='r')
plt.show(block=True)

plt.imshow(image[1,:,:], interpolation='nearest', cmap="Greys")
plt.colorbar()
ROI.displayROI()

mask = ROI.getMask(image[1,:,:])
output = np.zeros((len(image[:,1,1]),2))
timepoints = np.arange(len(image[:,1,1]))
output[:,0] = timepoints
for key in range(len(output[:,1])):
    output[key][1] = pl.mean(image[key,mask])

plt.show(block=True)