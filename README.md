# Gaborgen

Currently in this repository we are building up an fMRI preprocessing pipeline for the Gaborgen study. Eventually this may evolve into analyses for the EEG and fMRI data from this study. Base R is used so that there are no added dependencies besides AFNI and dcm2niix. It performs an initial regression with all stimuli to get nice quality control features which can be used for subsequent regressions and analyses.

To recreate the analyses, data should be taken from the dropbox in the same format as it is kept there for example:

GABORGEN24_101

- fMRI

- EEG

- dat

## 001fMRI_prepro.R
This script does the initial preprocessing per participant. It imposes restrictions to force the user to have the data in the correct format and to prevent incorrect analyses. It also forces the user to put their data and local git directory paths at the top of the script. The script forces the use of the tcsh terminal which is recommended for AFNI analyses. The system2 R function is used to force tcsh terminal calls.

