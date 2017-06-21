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
import easygui

plt.ioff() 

def openTiff():
    #Find the path to a .tiff via a GUI fileopen box
    print("What TIFF stack do you want to open?")
    filename = ""
    while not filename.endswith(".tiff"):
        filename = easygui.fileopenbox(filetypes=["*.tiff"])
    #print(filename)
    return filename

def roi_to_array(tiffstack):
    tiffstack = io.imread(tiffstack)
    plt.imshow(tiffstack[1,:,:])
    plt.title("Foreground: left click: line segment         right click: close region")
    fgROI = roipoly(roicolor='r')
    plt.show(block=True)
    plt.ioff()
    plt.imshow(tiffstack[1,:,:], interpolation='nearest', cmap="Greys")
    plt.title("Background: left click: line segment         right click: close region")
    plt.colorbar()
    fgROI.displayROI()
    bgROI = roipoly(roicolor='b')
    plt.show(block=True)
    fgmask = fgROI.getMask(tiffstack[1, :, :])
    bgmask = bgROI.getMask(tiffstack[1, :, :])
    output = np.zeros((len(tiffstack[:,1,1]),3))
    timepoints = np.arange(len(tiffstack[:,1,1]))
    output[:,0] = timepoints
    for key in range(len(output[:,1])):
        output[key][1] = pl.mean(tiffstack[key,fgmask])
        output[key][2] = pl.mean(tiffstack[key,bgmask])
    return output
 
# forexcel = pd.DataFrame(output[:,[1,2]], index = output[:,0], columns=["cell","bg"])
# forexcel.index.name = "time"
# forexcel.to_clipboard(excel=True)