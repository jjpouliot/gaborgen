# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas
import os
import numpy as np


sub = "120"
os.chdir("/Users/CSEA1/Documents/GaborGen24/Data/GABORGEN24_" + sub + "/DAT")


habCSplus = []
habGS1 = []
habGS2 = []
habGS3 = []

acqCSplusShock = []
acqCSplusNoShock = []
acqCSplus = []
acqGS1 = []
acqGS2 = []
acqGS3 = []

extCSplus = []
extGS1 = []
extGS2 = []
extGS3 = []

datStr = []

with open ("gaborgen24_fMRI_Day1_" + sub + "_logfile.dat") as datFile:
    data = pandas.read_csv(datFile, sep=",")

os.system("mkdir /Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli")

# %% Separate and store time series based on condition

#Habituation
hab = data[data['phase'] == 1]

habCSplus = hab[hab['stim'] == 1]
habCSplus['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/habCSplus_times.1D", index=False, header=False)

habGS1 = hab[hab['stim'] == 2]
habGS1['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/habGS1_times.1D", index=False, header=False)

habGS2 = hab[hab['stim'] == 3]
habGS2['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/habGS2_times.1D", index=False, header=False)

habGS3 = hab[hab['stim'] == 4]
habGS3['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/habGS3_times.1D", index=False, header=False)


#Acquisition
acq = data[data['phase'] == 2]

acqCSplus = acq[acq['stim'] == 1]
acqCSplus['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/acqCSplus_times.1D", index=False, header=False)

acqCSplusWShock = acqCSplus[acqCSplus['paired'] == 1]
acqCSplusWShock['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/acqCSplusWShock_times.1D", index=False, header=False)

acqCSplusNoShock = acqCSplus[acqCSplus['paired'] == 0]
acqCSplusNoShock['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/acqCSplusNoShock_times.1D", index=False, header=False)


acqGS1 = acq[acq['stim'] == 2]
acqGS1['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/acqGS1_times.1D", index=False, header=False)

acqGS2 = acq[acq['stim'] == 3]
acqGS2['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/acqGS2_times.1D", index=False, header=False)

acqGS3 = acq[acq['stim'] == 4]
acqGS3['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/acqGS3_times.1D", index=False, header=False)


#Extinction
ext = data[data['phase'] == 3]

extCSplus = ext[ext['stim'] == 1]
extCSplus['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/extCSplus_times.1D", index=False, header=False)

extGS1 = ext[ext['stim'] == 2]
extGS1['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/extGS1_times.1D", index=False, header=False)

extGS2 = ext[ext['stim'] == 3]
extGS2['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/extGS2_times.1D", index=False, header=False)

extGS3 = ext[ext['stim'] == 4]
extGS3['timeSinceFirstTR'].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/stimuli/extGS3_times.1D", index=False, header=False)


# %% Extract nuisance regressors

with open ("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/sub-" + sub + "_task-gaborgen_desc-confounds_timeseries.tsv") as confFile:
    confounds = pandas.read_csv(confFile, sep="\t")

headMots = ['trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'] 

for i in range(6):
    confounds[headMots[i]].to_csv("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/" + headMots[i] + "_run01.txt", index=False, header=False)
    
numTRS = confounds.shape[0]

goodTRs = np.ones((1070,1))

filter_col = [col for col in confounds if col.startswith('motion_outlier')]

motOutlierTRs = confounds[filter_col]
motOutlierTRsCol = motOutlierTRs.sum(axis=1)

badTRs_FWDbool = confounds['framewise_displacement'] > 0.3
badTRs_FWD = badTRs_FWDbool.replace({True: 1, False: 0})

badTRs_FWDnp = badTRs_FWD.to_numpy() / 2
motOutlierTRsColnp = motOutlierTRsCol.to_numpy() / 2
badTRs = np.round(badTRs_FWDnp + motOutlierTRsColnp, 0)

goodTRs = np.subtract(goodTRs, badTRs[:,None])
np.savetxt("/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/derivatives/sub-" + sub + "/func/motion_" + sub + "_censor.1D",goodTRs,delimiter=',',fmt='%i')
