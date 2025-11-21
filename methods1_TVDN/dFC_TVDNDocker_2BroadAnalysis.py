## libraries
from pyTVDN import TVDNDetect
from pathlib import Path
import numpy as np
import scipy.io
import os
import glob
import pandas as pd

## settings and directories
outputdir = '/out'
datadir = '/data'
AllTasks = ['BreathNVHA', 'BreathNVLA', 'Meditation']

## Get list of subjects within the datadir
AllSubj = glob.glob(datadir + '/*/')
AllSubj.sort()
for iDir in range(len(AllSubj)):
    AllSubj[iDir] = AllSubj[iDir][:-1] # remove the final "/" from each directory
    AllSubj[iDir] = AllSubj[iDir].replace("/data/","")
AllSubj.remove('sub-2021') # remove sub-2021. incomplete dataset


# create empty data frame
rvalues = np.linspace(0.4,0.9,6)
kappavalues = np.linspace(1,4,31)

for iRvalue in range(len(rvalues)):
    temprvalue = rvalues[iRvalue]
    temprvalue = np.round(temprvalue,1)
    # if the wider r-folder is not created, make directory
    temprpath = outputdir + '/r_' + str(temprvalue)
    if not os.path.exists(temprpath):
        os.mkdir(temprpath)  
    
    for ikappavalue in range(len(kappavalues)):
        tempkappavalue = kappavalues[ikappavalue]
        tempkappavalue = np.round(tempkappavalue,1)
        
        #if the wider kappa folder is not created, make directory
        tempkappapath = temprpath + '/k_' + str(tempkappavalue)
        if not os.path.exists(tempkappapath):
            os.mkdir(tempkappapath)
        
        #check if df1 is already made
        tempdfsavename = tempkappapath + '/df1.csv'
        if os.path.exists(tempdfsavename):
            tempdf = pd.read_csv(tempdfsavename)
        else:
            # create empty df1 dataframe
            tempdf = pd.DataFrame(columns = ['ID', 'Task', 'r', 'kappa', 'ECPTS', 'NumECPTS', 'CurMSE', 'ParasR'])
        
        for iSubj in range(len(AllSubj)):
            tempSubj = AllSubj[iSubj]
            # create subject directory
            tempsubjpath = tempkappapath + '/' + tempSubj
            if not os.path.exists(tempsubjpath):
                os.mkdir(tempsubjpath)
            for iTask in range(len(AllTasks)):
                tempTask = AllTasks[iTask]
                temptaskpath = tempsubjpath + '/' + tempTask
                if not os.path.exists(temptaskpath):
                    os.mkdir(temptaskpath)
                # check if final object is created or not; if created, skip
                if os.path.exists(temptaskpath + '/EigenCurve.jpg'):
                    continue
                #now read in the time series (TS) file
                tempTSfile = datadir + '/' + tempSubj + '/func/' + tempSubj + '_task-' + tempTask + '_space-MNI152NLin2009cAsym_atlas-4S156Parcels_timeseries.tsv'
                tempTS = pd.read_csv(tempTSfile, sep = '\t')
                #transpose the time series file into d x n matrix (d = ROI, n = time)
                tempTS = tempTS.transpose()
                # now drop all rows containing NaNs
                tempTS = tempTS.dropna()
                # begin TVDN detection
                Detection = TVDNDetect(Ymat = tempTS, saveDir = None, dataType = 'fMRI', fName = None, r = temprvalue, kappa = tempkappavalue, decimateRate = 1, freq = 1/0.85, lamb = 1e-4, Lmin = 2, downRate = 2, MaxM = 40, fct = 0.5, T = 2)
                Detection()
                # read outputs
                tempECPTS = Detection.ecpts
                tempNumECPTS = len(tempECPTS)
                tempCurMSE = Detection.GetCurMSE()
                tempParasR = Detection.paras['r']
                # place outputs into dataframe
                tempdf.loc[len(tempdf.index)] = [tempSubj, tempTask, temprvalue, tempkappavalue, tempECPTS, tempNumECPTS, tempCurMSE, tempParasR]
                tempdf.to_csv(tempdfsavename, index=False, mode='w')
                # plot outputs
                Detection.PlotEcpts(saveFigPath = temptaskpath + '/DetectionResults.jpg')
                Detection.PlotRecCurve(saveFigPath = temptaskpath + '/RecCurve.jpg')
                Detection.PlotEigenCurve(saveFigPath = temptaskpath + '/EigenCurve.jpg')
                print('r val:' + str(temprvalue) + ', of k val:' + str(tempkappavalue) + ' of ' + tempSubj + ' of ' + tempTask + ' completed.')
        


