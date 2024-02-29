# Gaborgen

Currently in this repository we are building up an fMRI preprocessing pipeline for the Gaborgen study. Eventually this may evolve into analyses for the EEG and fMRI data from this study.

To recreate the analyses, data should be taken from the dropbox in the same format as it is kept there for example:
GABORGEN24_101
- fMRI
- EEG
- dat



The master script imposes restrictions by making sure that the data is in the correct format and that the tcsh terminal is available (which is recommended for AFNI analyses). The system2 R function is used to force tcsh terminal calls. It also forces the user to put their data and local git directory paths at the top of the script.
