# -*- coding: utf-8 -*-
"""
Functions for Analyzing Voltage Dye Excel Data

Created on Tue Jun 13 16:24:12 2017

@author: Rishi
"""
import easygui
import numpy as np
import pandas as pd
import scipy.sparse.linalg as linalg
import scipy.sparse as sparse
from peakdet import peakdet
import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array
from scipy.signal import butter, lfilter
from scipy.signal import freqs
import scipy.signal as signal
import matplotlib.pyplot as plt

def workbooktoarray():
    #Find the path to a .xlsx via a GUI fileopen box
    print("What XLSX workbook do you want to open?")
    filename = ""
    while not filename.endswith(".xlsx"):
        filename = easygui.fileopenbox(filetypes=["*.xlsx"])        
    #print(filename)
    #Pick what sheet number you want to open
    print("What sheet number do you want to open? (The first sheet is 0)")
    while True:
        try:
            sheetnum = int(input())
        except ValueError:
            print("Sheet number must be an integer. What sheet number do you want to open? (The first sheet is 0)")
            continue
        break
    #Pick what columns to bring into the numpy array (A:C or A,B,C works fine)
    print("Which columns contain data? Input letters separated with commas.")
    parsecol = input()
    df = pd.read_excel(filename, sheet_name=sheetnum,header=0, usecols=parsecol)
    return df

def digitizeSpikes(spikearray):
    #this function will convert the filtered spikes to a 0-1 spike train, effectively digitizing the data
    digitize_spikes = np.zeros_like(spikearray) #this makes an array of 0s that will look like the input array
    digitize_spikes[:,0] = spikearray[:,0] #this sets the first column equal to time indices
    for cellnum in range(1,np.ma.size(spikearray,1)):
        threshold = 4.0*np.average(np.std(spikearray[:,1:],axis=0)) + np.median(spikearray[:,cellnum])
        above_thresh = np.where(spikearray[:,cellnum] > threshold)
        maxima = signal.argrelextrema(spikearray[:,cellnum], np.greater, order=2)
        spike_index = np.intersect1d(maxima, above_thresh)
        for index in spike_index:
            digitize_spikes[index,cellnum] = 1
    return digitize_spikes


def RasterSpikes(spikearray):
    #this function will convert the filtered spikes to a NaN-1 spike train, enabling easy plotting as a raster plot. If you plot this in matplotlib and place a line marker "|", you will get a traditional looking raster plot
    for i in range(1,np.ma.size(spikearray,1)):
        threshold = np.std(spikearray[:,i])
        rearm = 0
        for j in range(0,np.ma.size(spikearray,0)):
            if (spikearray[j,i] > 4.5*threshold):
                if (rearm == 0):
                    spikearray[j,i] = 1
                    rearm = 1
                else:
                    spikearray[j,i] = np.nan
            else:
                spikearray[j,i] = np.nan
                rearm = 0
    return spikearray

def plotRaster(RasterSpikes):
    #RasterSpikes, in this case, refers to an array that was produced via running the RasterSpikes() function on a filtered dataset containing noise and spikes. This function sets the x and y axes based on the size of the RasterSpikes array
    toplot = RasterSpikes[RasterSpikes > 0.1]
    for i in range(1,np.ma.size(RasterSpikes,1)):
        plt.plot(RasterSpikes[:,0], toplot[:,i] + i, label='baseline', marker='|', markersize=12)
    plt.xlim(0,np.ma.size(RasterSpikes[:,1]))
    plt.ylim(0,np.ma.size(RasterSpikes[1,:]) + 1)

def plotBaselineRaster(dSpikes, BaseNormal):
    #This function graphs the baseline values of all the cells in the BaseNormal array with an offset and then marks all the points that have an action potential with an upward line marker, effectively generating a raster plot that shows what the relative resting membrane potential was prior to that spike. The resting membrane potential is normalized to some value between 0 and 1 for ease of visualization. 
    for i in range(1,np.ma.size(BaseNormal,1)):
        markers = np.where(dSpikes[:,i] == 1)
        markers = np.asarray(markers)
        markers = markers.tolist()
        plt.plot(BaseNormal[:,0], i + ((BaseNormal[:,i]-np.min(BaseNormal[:,i]))/(np.max(BaseNormal[:,i]) - np.min(BaseNormal[:,i]))), label='baseline', markevery=markers, marker=2, markersize=12)
    plt.xlim(0,np.ma.size(BaseNormal[:,1]))
    
def plotBaseline(BaseNormal):
    #This function graphs the baseline values of all the cells in the BaseNormal array with an offset. The resting membrane potential is normalized to some value between 0 and 1 for ease of visualization.
    for i in range(1,np.ma.size(BaseNormal,1)):
        plt.plot(BaseNormal[:,0], i + ((BaseNormal[:,i]-np.min(BaseNormal[:,i]))/(np.max(BaseNormal[:,i]) - np.min(BaseNormal[:,i]))), label='baseline',)
    plt.xlim(0,np.ma.size(BaseNormal[:,1]))
    
def deltaF(baseline):
    for i in range(1,np.ma.size(baseline,1)):
        MaxFluor = np.mean(baseline[:,i])
        baseline[:,i] = np.divide(baseline[:,i],MaxFluor)
    return baseline
    
def backgroundSubtract(fluorarray):
    #generate a linear regression of the background, which is assumed to be column 3 of the array. column 1 is assumed to be timepoints.
    bgarray = fluorarray[:,[0,2]]
    y = bgarray[:,1]
    x = bgarray[:,0]
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A,y)[0]
    print(m,c)
    for keys in range(len(bgarray)):
        bgarray[keys][1] = m*keys + c
    sub=np.subtract(fluorarray[:,1],bgarray[:,1])
    fluorarray[:,1] = sub
    return fluorarray

def baseline_als(y, lam, p, niter=10):
    #performs a asymmetric least squares regression to figure out the baselinebas
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = linalg.spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

def baselineLowPass(rawdata):
    fs = 500  # Sampling frequency
    t = np.arange(2000) / fs
    rawx = rawdata
    fc = 3  # Cut-off frequency of the filter
    w = fc / (fs / 2) # Normalize the frequency
    b, a = signal.butter(5, w, 'low')
    output = signal.filtfilt(b, a, rawx)
    return output

def spikesHighPass(rawdata):
    fs = 500  # Sampling frequency
    t = np.arange(2000) / fs
    rawx = rawdata
    fc = 3  # Cut-off frequency of the filter
    w = fc / (fs / 2) # Normalize the frequency
    b, a = signal.butter(5, w, 'high')
    output = signal.filtfilt(b, a, rawx)
    return output

def baselineFlatten(bgcorrected):
    #divides out the baseline to produce dF/F plot. 
    ydata = bgcorrected[:,1]
    #you can change the fitting parameters for better fitting. 
    fitted = baseline_als(ydata, 10000, 0.0001)
#    fitted = baseline_als(fitted, 1000000, 1, niter=17)
    deltaf = np.divide(bgcorrected[:,1],fitted)
    bgcorrected[:,1] = deltaf
    return bgcorrected

def stepFinder(deltaFarray):
    #finds the steps in your trace via a "first derivative" and then uses them to calculate mean dF/F values
    zeroder = deltaFarray[:,1]
    smoothzd = smooth(zeroder,window_len=5,window="flat")
    firstder = smoothzd[1:] - smoothzd[:-1]
#    zeroder[1:]-zeroder[:-1]
#    
    maxpeaks, minpeaks = peakdet(firstder, (np.amax(firstder)/2))
    stepwidth = int(minpeaks[0][0] - maxpeaks[0][0])
    stepspace = (maxpeaks[2][0] - maxpeaks[0][0])/2
    firstpeak = int(maxpeaks[0][0])
    print("I think your step width is: " + str(stepwidth))
    print("Your first step is: " + str(firstpeak)) 
    print("Your spacing between steps is: "+str(stepspace))
    print("If I'm incorrect, enter new values, respectively. If I'm correct, enter blanks.")
    newstepwidth = input()
    newcurrentpeak = input()
    newstepspace = input()
    if newstepwidth != "":
        stepwidth = int(newstepwidth)
    if newcurrentpeak != "":
        firstpeak = int(newcurrentpeak)
    if newstepspace != "":
        stepspace = int(newstepspace)    
    currentpeak = firstpeak
    currentpeakend = currentpeak + stepwidth
    averagevalues = np.zeros((11,2))
    averagevalues[:,0] = [100, 80, 60, 40, 20, 0, -20, -40, -60, -80, -100]
    iterator = 0
    for step in range(len(averagevalues[:,0])):
        averagevalues[step][1] = ((np.mean(reject_outliers(zeroder[currentpeak:currentpeakend])) - 1)*100)
        print(str(currentpeak))
        iterator = iterator + 1
        currentpeak = firstpeak + int(iterator * stepspace)
        currentpeakend = currentpeak + stepwidth
    return averagevalues

def reject_outliers(data, m=4):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


    
#if __name__=="__main__":
#    from matplotlib.pyplot import plot, scatter, show
#    rawdata = workbooktoarray()
#    bgsub = backgroundSubtract(rawdata)
#    bleachcorr = baselineFlatten(bgsub)
#    dff = stepFinder(bleachcorr)
#    plot(dff[:,0],dff[:,1])
#    forexcel = pd.DataFrame(dff[:,1], index = dff[:,0], columns=["dF/F"])
#    forexcel.index.name = "mV"
#    forexcel.to_clipboard(excel=True)
#    print("Voltage Sensitivity Data copied to clipboard")