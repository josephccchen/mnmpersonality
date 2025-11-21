## directories
#output = "/Users/josephchen/data/MnM/MnM-dFC-TVDN"
#ratingsdir = "/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/data/MnMscripts/Behavioural"
#datadir = "/Users/josephchen/data/MnM/MnM-XCPD-NoGSRNoCensor"
#scriptdir = "/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/data/MnMscripts/fMRI_dFC"

from pyTVDN import TVDNDetect
from pathlib import Path
import numpy as np
import scipy.io
import os
import glob
import pandas as pd
from numpy.linalg import inv
## settings
outputdir = '/out'
datadir = '/data'
#ratingsdir = '/ratings'
AllTasks = ['BreathNVHA', 'BreathNVLA', 'Meditation']


## set a function to obtain FC matrix
def obt_absFC(Amat):
    """get abs FC
    (1) transform to correlation
    (2) z transform
    (3) take absolute value
    (remove the diaginal term)
    """
    FC = np.corrcoef(Amat) # Pearson's r
    FC = FC - np.diag(np.diag(FC)) # remove diagonal
    FC = np.arctanh(FC) # Fisher's transform
    abs_FC = np.abs(FC)
    return abs_FC

# function to extract the upper triangular part of the matrix
def extract_upper_triangle(matrix):
    upper_tri_indices = np.triu_indices(matrix.shape[0], k=1)
    upper_tri_values = matrix[upper_tri_indices]
    return upper_tri_values, upper_tri_indices


# Reconstruct the symmetric matrix from the 1xN array 
def reconstruct_matrix_from_upper_tri(modified_values, upper_tri_indices, matrix_size):
    reconstructed_matrix = np.zeros((matrix_size, matrix_size))
    reconstructed_matrix[upper_tri_indices] = modified_values
    reconstructed_matrix.T[upper_tri_indices] = modified_values  # Symmetrize the matrix
    return reconstructed_matrix

## Get list of subjects within the datadir
AllSubj = glob.glob(datadir + '/*/')
AllSubj.sort()
for iDir in range(len(AllSubj)):
    AllSubj[iDir] = AllSubj[iDir][:-1] # remove the final "/" from each directory
    AllSubj[iDir] = AllSubj[iDir].replace("/data/","")


# remove sub-2021
AllSubj.remove('sub-2021')

# create empty data frame
finalvalues = pd.DataFrame(np.array([[0.5, 1.5],
                                     [0.5, 1.6],
                                     [0.6, 1.6],
                                     [0.6, 1.7],
                                     [0.7, 1.4],
                                     [0.7, 1.5], 
                                     [0.7, 1.6], 
                                     [0.8, 1.2],
                                     [0.8, 1.3],
                                     [0.9, 1.3],
                                     [0.9, 1.4]]),
                           columns = ['r', 'kappa'])
NumRows = np.shape(finalvalues)[0]


finaldir = outputdir + '/final_analyses'
# check if finaldir is made; if not, make it
if not os.path.exists(finaldir):
    os.mkdir(finaldir)

for iRow in range(NumRows):
    #final r + k values
    tempr = finalvalues['r'][iRow]
    tempkappa = finalvalues['kappa'][iRow]
    # working directory
    newdir = finaldir + '/final_r' + str(tempr) + '_kappa' + str(tempkappa)
    if not os.path.exists(newdir):
        os.mkdir(newdir)
    #check if df1 is already made
    df1savename = newdir + '/df1.csv'
    
    # create empty df1 dataframe
    df1 = pd.DataFrame(columns = ['ID', 'Task', 'r', 'kappa', 'ECPTS', 'NumECPTS', 'CurMSE', 'ParasR', 'NumMatrices'])
    
    for iSubj in range(len(AllSubj)):
        tempSubj = AllSubj[iSubj]
        # create subject directory
        tempsubjpath = newdir + '/' + tempSubj
        if not os.path.exists(tempsubjpath):
            os.mkdir(tempsubjpath)
        for iTask in range(len(AllTasks)):
            tempTask = AllTasks[iTask]
            temptaskpath = tempsubjpath + '/' + tempTask
            if not os.path.exists(temptaskpath):
                os.mkdir(temptaskpath)
            
            # check if final object is created or not; if created, skip
            #if os.path.exists(temptaskpath + '/EigenCurve.jpg'):
            #    continue
            #now read in the time series (TS) file
            tempTSfile = datadir + '/' + tempSubj + '/func/' + tempSubj + '_task-' + tempTask + '_space-MNI152NLin2009cAsym_atlas-4S156Parcels_timeseries.tsv'
            tempTS = pd.read_csv(tempTSfile, sep = '\t')
            #transpose the time series file into d x n matrix (d = ROI, n = time)
            tempTS = tempTS.transpose()
            # now drop all rows containing NaNs
            tempTS = tempTS.dropna()
            # begin TVDN detection
            #Detection = TVDNDetect(Ymat=Ymat, saveDir=None, dataType="fMRI", fName=None, r=20, kappa=1.5, decimateRate = 1, freq=0.85, lamb = 1e-6, Lmin = 2, downRate = 2, MaxM=30, fct = 0.5, T = 2)
            #Detection = TVDNDetect(Ymat = tempTS, saveDir = None, dataType = 'fMRI', fName = None, r = 0.6, kappa = 2, decimateRate = 1, freq = 1/0.85, lamb = 1e-4, Lmin = 2, downRate = 2, MaxM = 30, fct = 0.5, T = 2)
            Detection = TVDNDetect(Ymat = tempTS, saveDir = None, dataType = 'fMRI', fName = None, r = tempr, kappa = tempkappa, decimateRate = 1, freq = 1/0.85, lamb = 1e-4, Lmin = 2, downRate = 2, MaxM = 40, fct = 0.5, T = 2)
            Detection()
            
            #get A(t) for each segment
            Detection.GetFeatures()
            U = Detection.midRes.eigVecs
            r = Detection.paras.r
            Urinv = inv(U)[:r, :]
            Ur = U[:, :r]
            lamMs = Detection.curEigVals
            As = []
            for lamM in lamMs:
                At = Ur @ np.diag(lamM) @ Urinv
                As.append(At)
            
            #for each A(t), print out the conn matrix
            for iFC in range(len(As)):
                tempAt = As[iFC]
                tempFC = obt_absFC(tempAt)
                # now print out each 
                tempFCnum = str(iFC+1).zfill(2)
                tempFCsavename = temptaskpath + "/FC_matrix-" + tempFCnum + '.csv'
                np.savetxt(tempFCsavename, tempFC, delimiter=",")
            
            # read outputs
            tempECPTS = Detection.ecpts
            tempNumECPTS = len(tempECPTS)
            tempCurMSE = Detection.GetCurMSE()
            tempParasR = Detection.paras['r']
            tempNumAs = len(As)
            # place outputs into dataframe
            df1.loc[len(df1.index)] = [tempSubj, tempTask, tempr, tempkappa, tempECPTS, tempNumECPTS, tempCurMSE, tempParasR, tempNumAs]
    df1.to_csv(df1savename, index=False, mode='w')
    # plot outputs
    Detection.PlotEcpts(saveFigPath = temptaskpath + '/DetectionResults.jpg')
    Detection.PlotRecCurve(saveFigPath = temptaskpath + '/RecCurve.jpg')
    Detection.PlotEigenCurve(saveFigPath = temptaskpath + '/EigenCurve.jpg')
    print('r val:' + str(tempr) + ', of k val:' + str(tempkappa) + ' of ' + tempSubj + ' of ' + tempTask + ' completed.')



