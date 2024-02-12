#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 22:53:37 2024

@author: CSEA
"""


## Convert Data to BIDS format
from bids import BIDSLayout
from bids.tests import get_test_data_path
import os

bidsDir = "/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen";
bidsConfig = "/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/code/BIDS_config.json"
subSourceDir = "/Users/CSEA1/Documents/GaborGen24/Data/GABORGEN24_118/fMRI"

bidsSubID = "118"


os.chdir(bidsDir)

os.system("dcm2bids_scaffold --force")

os.system("dcm2bids_helper -d " + subSourceDir + " --force")

os.system("dcm2bids -d " + subSourceDir + " -p " + bidsSubID + " -c " + bidsConfig + " --auto_extract_entities")



## Inspect images in nibabel
#from nilearn import plotting
#from matplotlib import pyplot as plt
#import numpy as np
#import nibabel as nb

#funcDir = "/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/sub-119/func"
#anatDir = "/Users/CSEA1/Documents/GaborGen24/Data/bidsGaborGen/sub-119/anat"

#t1w = anatDir + "/sub-119_T1w.nii.gz"

#plotting.plot_glass_brain(t1w, threshold=3, colorbar=True,
#                          title='plot_glass_brain with display_mode="lyrz"',
#                          plot_abs=False, display_mode='lyrz', cmap='rainbow')