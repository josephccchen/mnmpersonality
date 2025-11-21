## libraries
from pyTVDN import TVDNDetect
from pathlib import Path
import numpy as np
import scipy.io
import os
import glob
import pandas as pd
from numpy.linalg import inv

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
AllSubj.remove('sub-2021') # remove sub-2021 (incomplete dataset)

# create empty data frame
rvalues = np.linspace(0.4,0.9,6)
kappavalues = np.linspace(1,4,31)

masterdf = pd.DataFrame(columns = ['rvalue', 'kappavalue', 'TotalCurMSE', 'NumZeroECPTS', 'AverageECPTS', 'NumRows'])

# now read through all df1.csv and create masterdf
for irvalues in range(len(rvalues)):
    temprvalue = rvalues[irvalues]
    temprvalue = temprvalue.round(1)
    temprpath = outputdir + '/r_' + str(temprvalue)
    for ikappavalue in range(len(kappavalues)):
        tempkappavalue = kappavalues[ikappavalue]
        tempkappavalue = np.round(tempkappavalue,1)
        tempkappapath = temprpath + '/k_' + str(tempkappavalue)
        tempdffile = tempkappapath + '/df1.csv'
        #read df1
        df1 = pd.read_csv(tempdffile)
        
        #calculate total CurMSE
        tempTotalCurMSE = sum(df1.CurMSE)
        
        # calculate Number of Zero ECPTS
        tempNumZeros = df1['NumECPTS'].value_counts().get(0, 0)
        
        # calculate Average ECPTS
        tempAvgECPTS = df1['NumECPTS'].mean()
        
        #get number of rows of df1
        tempshape = df1.shape
        tempRows = tempshape[0]
        
        # append to masterdf
        masterdf.loc[len(masterdf.index)] = [temprvalue, tempkappavalue, tempTotalCurMSE, tempNumZeros, tempAvgECPTS, tempRows]

masterdfsavename = outputdir + '/masterdf.csv'
masterdf.to_csv(masterdfsavename)

