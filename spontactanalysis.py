"""
Script for importing an XLSX worksheet of fluorescence traces and creating a new workbook with the following format:

Worksheet 1: Raw Data
Worksheet 2: Baseline, which contains slow depolarization info
Worksheet 3: Baseline-corrected spiking data

Created by Rishi Kulkarni
"""

import voltagestepfunctions as VI
import numpy as np
import pandas as pd
import os
import easygui
import itertools
import sklearn
import sklearn.preprocessing
from sklearn.decomposition import FactorAnalysis
from factor_analyzer import FactorAnalyzer


def workbooktoDF():
    #Find the path to a .xlsx via a GUI fileopen box 
    print("What XLSX workbook do you want to open?")
    global filename
    filename = ""
    while not filename.endswith(".xlsx"):
        filename = easygui.fileopenbox(filetypes=["*.xlsx"])        
    #print(filename)
    df = pd.read_excel(filename, sheet_name=None)
    return df

def spike_triggered_baseline(normalized_baseline, spiketrain):
    #this function takes a continuous baseline function and a digitized spike train and outputs a matrix consisting of the time points +/- t on the baseline from each of the spikes in the spike train. It normalizes the added values to be dF/F values.
    t = 150 #this is the number of frames away from the spike the function will average
    spike_trigger_baseline = np.zeros((2*t,1))
    spike_trigger_baseline[:,0] = np.arange(-t,t)
    markers = np.asarray(np.where(spiketrain == 1))
    markers = markers[:,markers[0]>t]
    markers = markers[:,markers[0]<(np.ma.size(spiketrain[:,0])-t)]
    for cellnum in range(1, np.max(markers[1])+1):
        indices = np.where(markers[1] == cellnum)
        for spiketime in markers[:,indices[0]][0]:
            toadd = (normalized_baseline[(spiketime-t):(spiketime+t),cellnum])/np.mean(normalized_baseline[:,cellnum])
            spike_trigger_baseline = np.hstack((spike_trigger_baseline, toadd.reshape(-1,1)))
    return spike_trigger_baseline
            
def spike_triggered_average(sta_matrix):
    sta_average = np.zeros((300,2))
    sta_stdev = np.zeros((300,2))
    sta_average[:,0] = sta_matrix[:,0]
    sta_stdev[:,0] = sta_matrix[:,0]
    for i in range(np.ma.size(sta_matrix[:,0])):
        sta_average[i,1] = np.mean(sta_matrix[i,1:])
    for i in range(np.ma.size(sta_matrix[:,0])):
        sta_stdev[i,1] = np.std(sta_matrix[i,1:])
    return sta_average, sta_stdev

def pairwise_cc(baseline):
    baselineDF = pd.DataFrame(data=baseline)
    cross_corr_matrix = np.zeros((len(baselineDF.columns[1:]),len(baselineDF.columns[1:])))
    cross_corr_lag = np.zeros((len(baselineDF.columns[1:]),len(baselineDF.columns[1:])))
    for cell1, cell2 in itertools.combinations(baselineDF.columns[1:], 2):
        cell1base = baselineDF[cell1].values
        cell2base = baselineDF[cell2].values
        time = len(cell1base)
        cell1base = (cell1base - np.mean(cell1base)) / (np.std(cell1base) * len(cell1base))
        cell2base = (cell2base - np.mean(cell2base)) / (np.std(cell2base))
        fullcorr = np.correlate(cell1base, cell2base,'full')[time-50:time+50]
        corr_value = np.max(fullcorr)
        corr_lag = np.argmax(fullcorr) - 51
        cross_corr_matrix[cell1-1,cell2-1] = corr_value
        cross_corr_lag[cell1-1,cell2-1] = corr_lag
    return cross_corr_matrix, cross_corr_lag

def factor_analysis(baseline):
    norm_base = np.copy(baseline)
    norm_base[:,1:] = sklearn.preprocessing.StandardScaler().fit_transform(norm_base[:,1:])
    #these commands copy the baseline array and scales the data that isn't part of the x-axis
    factor = FactorAnalyzer(n_factors=3,rotation=None)
    factor.fit(norm_base[:,1:])
    eigenvalues = factor.get_eigenvalues()
    variance = factor.get_factor_variance()
    communalities = factor.get_communalities()
    loadings = factor.loadings_
    transform = factor.transform(norm_base[:,1:])
    outputFA = pd.DataFrame(data=transform,columns=["Factor 1","Factor 2","Factor 3"])
    outputFA["Eigenvalues (Raw) (index is factor #)"] = pd.Series(eigenvalues[0])
    outputFA["Eigenvalues (Common Factor)"] = pd.Series(eigenvalues[1])
    outputFA["Factor 1 Loadings (index is cell #)"] = pd.Series(loadings[:,0])
    outputFA["Factor 2 Loadings"] = pd.Series(loadings[:,1])
    outputFA["Factor 3 Loadings"] = pd.Series(loadings[:,2])
    outputFA["Communalities"] = pd.Series(communalities)
    outputFA["Proportional Explained Variance (index is factor #)"] = pd.Series(variance[1])
    return outputFA

    

rawDF = workbooktoDF()
for key in rawDF:
    folder, outputfile = os.path.split(filename)
    rawdata = rawDF[key].values
    for i in range(1,np.ma.size(rawdata,1)):
        rawdata[:,i] = rawdata[:,i] - VI.baseline_als(rawdata[:,i], 100000000,0.001)    
    #this corrects out any linear photobleaching, which is important for some dyes
    baseline = np.copy(rawdata)
    spikes = np.copy(rawdata)
    
    for i in range(1,np.ma.size(baseline,1)):
        actionpotentials = VI.spikesHighPass(spikes[:,i])
        fitted = VI.baseline_als(baseline[:,i], 10000,0.01)
        baseline[:,i] = fitted
        spikes[:,i] = spikes[:,i] - fitted

    dSpikes = np.copy(spikes)
    dSpikes = VI.digitizeSpikes(dSpikes)
    BaseNormal = np.copy(baseline)
    BaseNormal = VI.deltaF(BaseNormal)
    RasterSpikes = np.copy(spikes)
    RasterSpikes = VI.RasterSpikes(RasterSpikes)
    FactorAnalysisDF = factor_analysis(baseline)
    if key != "Concatenated":
        cross_corr_matrix, cross_corr_lag = pairwise_cc(baseline)
    
    baselineDF = pd.DataFrame(data=baseline, index = rawdata[:,0], columns = rawDF[key].columns)
    BaseNormalDF = pd.DataFrame(data=BaseNormal, index = rawdata[:,0], columns = rawDF[key].columns)
    spikesDF = pd.DataFrame(data=spikes, index = rawdata[:,0], columns = rawDF[key].columns)
    dSpikesDF = pd.DataFrame(data=dSpikes, index = rawdata[:,0], columns = rawDF[key].columns)
    if key != "Concatenated":
        cross_corr_matrix = pd.DataFrame(data=cross_corr_matrix, index = rawDF[key].columns[1:], columns = rawDF[key].columns[1:])
        cross_corr_lag = pd.DataFrame(data=cross_corr_lag, index = rawDF[key].columns[1:], columns = rawDF[key].columns[1:])
    
    if not os.path.exists(folder+"\\"+os.path.splitext(outputfile)[0]):
        os.mkdir(folder+"\\"+os.path.splitext(outputfile)[0])

    with pd.ExcelWriter(folder+"\\"+os.path.splitext(outputfile)[0]+"\\"+key+".xlsx") as writer:
        rawDF[key].to_excel(writer, sheet_name="Raw Data")
        baselineDF.to_excel(writer, sheet_name="Baseline_Filtered")
        BaseNormalDF.to_excel(writer, sheet_name="Baseline_Normalized")
        spikesDF.to_excel(writer, sheet_name="Spikes_Filtered")
        dSpikesDF.to_excel(writer,sheet_name="Spike Train")
        FactorAnalysisDF.to_excel(writer,sheet_name="Factor Analysis")
        if key != "Concatenated":
            cross_corr_matrix.to_excel(writer,sheet_name="Baseline_CC_R_Values")
            cross_corr_lag.to_excel(writer,sheet_name="Baseline_CC_Time_Lags")


print("Done!")