import voltagestepfunctions as VI
import seaborn as sns; sns.set()
import numpy as np
import pandas as pd
import os
import easygui
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
import scipy
import matplotlib.pyplot as plt
import itertools
import sklearn
import copy
import multiprocessing
import threading
import time 
import sklearn.preprocessing
from factor_analyzer import FactorAnalyzer
from numba import njit, jit


def baseline_als(y, lam, p, niter=5):
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

def fit_column(dfColumn, lam = 5000, weight = 0.01, niter=10):
    fitted = baseline_als(dfColumn, lam, weight, niter)
    return fitted

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

def digitizeSpikes(spikearray):
    #this function will convert the filtered spikes to a 0-1 spike train, effectively digitizing the data
    thresh = 6*(pd.Series(spikearray).rolling(250).std().min()) +  pd.Series(spikearray).rolling(500).median().median()
    above_thresh = np.where(spikearray > thresh)
    maxima = signal.argrelextrema(spikearray, np.greater, order = 2)
    spike_index = np.intersect1d(maxima, above_thresh)
    spiketrain = np.copy(spikearray)
    spiketrain[spike_index] = 1
    spiketrain = (spiketrain == 1).astype(int)
    return spiketrain

def deltaF(baseline):
    for i in range(0,np.ma.size(baseline,1)):
        MaxFluor = np.mean(baseline[:,i])
        baseline[:,i] = np.divide(baseline[:,i],MaxFluor)
    return baseline

def plotBaseline(BaseNormal):
    #This function graphs the baseline values of all the cells in the BaseNormal array with an offset. The resting membrane potential is normalized to some value between 0 and 1 for ease of visualization.
    for i in range(0,np.ma.size(BaseNormal,0)-1):
        plt.plot(BaseNormal.index, i + ((BaseNormal.values[:,i]-np.min(BaseNormal.values[:,i]))/(np.max(BaseNormal.values[:,i]) - np.min(BaseNormal.values[:,i]))), label='baseline',)
    plt.xlim(0,np.ma.size(BaseNormal.index))
    
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

def generate_ISI(spiketrain, subthresh):      
    ISI = []
    ISI_baseline = []
    for i in range(0,len(spiketrain.values[0,:])):
        for j in range(0,np.ma.size(np.where(spiketrain.iloc[:,i])[0])-1):
            ISI_baseline = np.append(ISI_baseline, np.mean(subthresh.iloc[np.where(spiketrain.iloc[:,i])[0][j]:np.where(spiketrain.iloc[:,i])[0][j+1],i]))
        ISI_add = np.diff(np.where(spiketrain.iloc[:,i]))
        #ISI_baseline = np.delete(ISI_baseline, -1, 0)
        ISI = np.append(ISI,ISI_add)
    return ISI, ISI_baseline

    
    
class VoltageTraceData:
    
        def __init__(self):
            self.manager = multiprocessing.Manager()
        
        
        def import_data(self, datatype="Excel"):
            print(datatype)
            if datatype == "Excel":
                self.RawData = workbooktoDF()

        def preprocess_data(self):
            self.process_raw()
            self.normalize_subthreshold()
            self.make_spiketrain()

        def subtract_bg(self):
            if 'BG' in self.RawData.columns.levels[1]:
                for key in self.RawData.columns.unique(level=0):
                    BGfit = baseline_als(self.RawData[key]['BG'], 1000000000,0.5)
                    self.RawData[key] = self.RawData[key].drop(columns='BG')
                    self.RawData[key] = self.RawData[key] - BGfit[:,None]
                self.RawData = self.RawData.dropna(axis=1)
            else:
                print("No BG trace!")

        def process_raw(self, smooth=5000, weight = 0.01, iterat=10):
            self.SubThreshold = self.RawData.transform(baseline_als, axis=0, raw=True, lam = smooth, p = weight, niter=iterat)
            self.Spikes = self.RawData.divide(self.SubThreshold, axis="columns")
   
        def make_spiketrain(self):
            self.SpikeTrain = self.Spikes.transform(digitizeSpikes, axis=0, raw=True)
                
        def normalize_subthreshold(self):
            self.NormSubThresh = self.SubThreshold.divide(self.SubThreshold.mean(axis=0))
             
        def gen_STA(self):
            self.STA = self.NormSubThresh.copy(deep=True)
            self.STA = self.STA.iloc[:300]
            self.STA[:] = spike_triggered_baseline(self.NormSubThresh.values, self.SpikeTrain.values)
            self.STA.values[self.STA.values == 0] = np.nan
        
        def FactorAnalyze(self, rotate='varimax'):
            self.SharedVariance = copy.deepcopy(self.NormSubThresh)
            self.SharedVariance = self.SharedVariance.iloc[:1]
            self.SharedVariance.index.rename('Shared Variance', inplace=True)
            self.FactorLoadings = copy.deepcopy(self.NormSubThresh)
            self.FactorLoadings = self.FactorLoadings.iloc[:3]
            self.FactorLoadings.index.rename('Factor Number', inplace=True)
            for key in self.SharedVariance.columns.unique(level=0):
                factor = FactorAnalyzer(n_factors=3,rotation=rotate)
                factor.fit(sklearn.preprocessing.StandardScaler().fit_transform(self.NormSubThresh[key].values))
                self.SharedVariance[key] = np.atleast_2d(factor.get_communalities())
                self.FactorLoadings[key] = factor.loadings_.T