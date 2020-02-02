import numpy as np
import pandas as pd
import os
import easygui
import numpy as np
import pandas as pd
import scipy.sparse.linalg as linalg
import scipy.sparse as sparse
import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array
import scipy.signal as signal
import scipy
import itertools
import sklearn
import copy
import multiprocessing
import threading
import time 
import sklearn.preprocessing
from factor_analyzer import FactorAnalyzer


def threaded(fn):
    def wrapper(*args, **kwargs):
        thr = threading.Thread(target=fn, args=args, kwargs=kwargs)
        thr.start()
        return thr
    return wrapper

def workbooktoDF():
    #Find the path to a .xlsx via a GUI fileopen box 
    print("What XLSX workbook do you want to open?")
    global filename
    filename = ""
    while not filename.endswith(".xlsx"):
        filename = easygui.fileopenbox(filetypes=["*.xlsx"])        
    t0 = time.time()
    df = pd.concat(pd.read_excel(filename, sheet_name=None,index_col=0),axis=1)
    t1 = time.time()
    print(t1-t0)
    return df

def baseline_als(y, lam, p, niter=10):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w) # Do not create a new matrix, just update diagonal values
        Z = W + D
        z = linalg.spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

def digitizeSpikes(spikearray):
    #this function will convert the filtered spikes to a 0-1 spike train, effectively digitizing the data
    digitize_spikes = np.zeros_like(spikearray) #this makes an array of 0s that will look like the input array
    for cellnum in range(0,np.ma.size(spikearray,1)):
        threshold = 4.5*np.average(np.min(np.std(rolling_window(spikearray[:,cellnum],250),axis=-1))) + np.median(np.median(rolling_window(spikearray[:,cellnum],500),axis=-1))
        above_thresh = np.where(spikearray[:,cellnum] > threshold)
        maxima = signal.argrelextrema(spikearray[:,cellnum], np.greater, order=2)
        spike_index = np.intersect1d(maxima, above_thresh)
        for index in spike_index:
            digitize_spikes[index,cellnum] = 1
    return digitize_spikes

def deltaF(baseline):
    for i in range(0,np.ma.size(baseline,1)):
        MaxFluor = np.mean(baseline[:,i])
        baseline[:,i] = np.divide(baseline[:,i],MaxFluor)
    return baseline
    
def spike_triggered_baseline(normalized_baseline, spiketrain):
    #this function takes a continuous baseline function and a digitized spike train and outputs a matrix consisting of the time points +/- t on the baseline from each of the spikes in the spike train. It normalizes the added values to be dF/F values.
    t = 150 #this is the number of frames away from the spike the function will average
    spike_trigger_baseline = np.zeros((2*t,len(spiketrain[0,:])))
    markers = np.asarray(np.where(spiketrain == 1))
    markers = markers[:,markers[0]>t]
    markers = markers[:,markers[0]<(np.ma.size(spiketrain[:,0])-t)]
    for cellnum in range(0, np.max(markers[1])+1):
        indices = np.where(markers[1] == cellnum)
        for spiketime in markers[:,indices[0]][0]:
            toadd = (normalized_baseline[(spiketime-t):(spiketime+t),cellnum])
            spike_trigger_baseline[:,cellnum] = spike_trigger_baseline[:,cellnum] + toadd
        np.divide(spike_trigger_baseline[:,cellnum],(np.count_nonzero((markers[1,:] == cellnum))), out=spike_trigger_baseline[:,cellnum], where=(np.count_nonzero((markers[1,:] == cellnum)))!=0)
    return spike_trigger_baseline

def rolling_window(a, window):
    pad = np.ones(len(a.shape), dtype=np.int32)
    pad[-1] = window-1
    pad = list(zip(pad, np.zeros(len(a.shape), dtype=np.int32)))
    a = np.pad(a, pad,mode='reflect')
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

class VoltageTraceData:
    
        def __init__(self):
            self.manager = multiprocessing.Manager()
        
        
        def import_data(self, datatype="Excel"):
            print(datatype)
            if datatype == "Excel":
                self.RawData = workbooktoDF()
            self.SubThreshold = copy.deepcopy(self.RawData)
            self.Spikes = copy.deepcopy(self.RawData)
            
        @threaded  
        def preprocess(self, lam=5000, weight=0.01):
            for i in range(0,np.ma.size(self.SubThreshold.values,1)):
                fitted = baseline_als(self.SubThreshold.values[:,i], lam, weight)
                self.SubThreshold[self.SubThreshold.columns.get_level_values(0)[i],self.SubThreshold.columns.get_level_values(1)[i]] = fitted
                self.Spikes[self.Spikes.columns.get_level_values(0)[i],self.Spikes.columns.get_level_values(1)[i]] = self.Spikes[self.Spikes.columns.get_level_values(0)[i],self.Spikes.columns.get_level_values(1)[i]].values / fitted
        
        @threaded  
        def detrend(self, lam=1000000000000, weight=0.5):
            self.SubThresh_Stationary = copy.deepcopy(self.NormSubThresh)
            for i in range(0,np.ma.size(self.SubThresh_Stationary.values,1)):
                fitted = baseline_als(self.SubThresh_Stationary.values[:,i], lam, weight)
                self.SubThresh_Stationary[self.SubThresh_Stationary.columns.get_level_values(0)[i],self.SubThresh_Stationary.columns.get_level_values(1)[i]] = self.SubThresh_Stationary[self.SubThresh_Stationary.columns.get_level_values(0)[i],self.SubThresh_Stationary.columns.get_level_values(1)[i]] - fitted + 1      
        
        
        @threaded            
        def make_spiketrain(self):
            self.SpikeTrain = copy.deepcopy(self.Spikes)
            for key in self.SpikeTrain.columns.unique(level=0):
                self.SpikeTrain[key] = digitizeSpikes(self.SpikeTrain[key].values)
                
                
        @threaded
        def normalize_subthreshold(self):
            self.NormSubThresh = copy.deepcopy(self.SubThreshold)
            for key in self.NormSubThresh.columns.unique(level=0):
                self.NormSubThresh[key] = deltaF(self.NormSubThresh[key].values)
             
        @threaded
        def gen_STA(self):
            self.STA = self.NormSubThresh.copy(deep=True)
            self.STA = self.STA.iloc[:300]
            self.STA[:] = spike_triggered_baseline(self.NormSubThresh.values, self.SpikeTrain.values)
            self.STA.values[self.STA.values == 0] = np.nan
        
        @threaded
        def FactorAnalyze(self, rotate='varimax'):
            self.SharedVariance = copy.deepcopy(self.SubThresh_Stationary)
            self.SharedVariance = self.SharedVariance.iloc[:1]
            self.SharedVariance.index.rename('Shared Variance', inplace=True)
            self.FactorLoadings = copy.deepcopy(self.SubThresh_Stationary)
            self.FactorLoadings = self.FactorLoadings.iloc[:3]
            self.FactorLoadings.index.rename('Factor Number', inplace=True)
            for key in self.SharedVariance.columns.unique(level=0):
                factor = FactorAnalyzer(n_factors=3,rotation=rotate)
                factor.fit(self.SubThresh_Stationary[key].values)
                self.SharedVariance[key] = np.atleast_2d(factor.get_communalities())
                self.FactorLoadings[key] = factor.loadings_.T