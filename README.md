# Gaborgen

Currently in this repository we are building up an fMRI preprocessing pipeline for the Gaborgen study. Eventually this may evolve into analyses for the EEG and fMRI data from this study. Base R is used so that there are no added dependencies besides AFNI and dcm2niix. It performs an initial regression with all stimuli to get nice quality control features which can be used for subsequent regressions and analyses.

To install dcm2niix on a mac (https://formulae.brew.sh/formula/dcm2niix), open a terminal and run:
brew install dcm2niix

If you don't have homebrew, look at the instructions for installing it (https://brew.sh/)

Installing AFNI is done through the terminal as well. There are various options for doing so depending on your operating system (https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/background_install/install_instructs/index.html). Double check this website to make sure terminal commands below are still up to date.

For a Mac, first check if your Mac as an ARM or Intel CPU. 

If Intel, open terminal and run:
cd
curl -O https://raw.githubusercontent.com/afni/afni/master/src/other_builds/OS_notes.macos_12_intel_a_admin_pt1.zsh
curl -O https://raw.githubusercontent.com/afni/afni/master/src/other_builds/OS_notes.macos_12_intel_a_admin_pt2.zsh
curl -O https://raw.githubusercontent.com/afni/afni/master/src/other_builds/OS_notes.macos_12_intel_b_user.tcsh

zsh OS_notes.macos_12_intel_a_admin_pt1.zsh

zsh OS_notes.macos_12_intel_a_admin_pt2.zsh

REBOOT YOUR COMPUTER NOW

tcsh OS_notes.macos_12_intel_b_user.tcsh


If ARM (M1 CPU), open terminal and run:
cd
curl -O https://raw.githubusercontent.com/afni/afni/master/src/other_builds/OS_notes.macos_12_ARM_a_admin_pt1.zsh
curl -O https://raw.githubusercontent.com/afni/afni/master/src/other_builds/OS_notes.macos_12_ARM_a_admin_pt2.zsh
curl -O https://raw.githubusercontent.com/afni/afni/master/src/other_builds/OS_notes.macos_12_ARM_b_user.tcsh

zsh OS_notes.macos_12_ARM_a_admin_pt1.zsh

zsh OS_notes.macos_12_ARM_a_admin_pt2.zsh

REBOOT YOUR COMPUTER NOW

tcsh OS_notes.macos_12_ARM_b_user.tcsh

Everything should now should now run if everything is installed correctly and dat ais in the correct structure as described below. If this doesn't run well and you get x11 errors or no qualty checks outputs, then you should look at the home_terminal_configuration_files/ folder. There are examples of .cshrc & .tcshrc configuration files that should be in your home directory. If it looks like an empty folder, that is because you need to reveal "hidden" files as files that start with a period are usually not visible by default. 

To recreate the analyses, data should be taken from the dropbox in the same format as it is kept there for example:

GABORGEN24_101

- fMRI

- EEG

- DAT

Participants 123 and up have this format because there are two potential days of recording

GABORGEN24_DAY1_123

- fMRI

- EEG

- DAT

GABORGEN24_DAY2_123

- fMRI

- EEG

- DAT

## a001fMRI_prepro.R
This script does the initial preprocessing per participant. It imposes restrictions to force the user to have the data in the correct format and to prevent incorrect analyses. It also forces the user to put their data and local git directory paths at the top of the script. The script forces the use of the tcsh terminal which is recommended for AFNI analyses. The system2 R function is used to force tcsh terminal calls.

### Steps for fMRI preprocessing

- Make sure the raw data are in two locations: 1) connected passport hard-drive (/Volumes/MY_PASSPORT/gaborgen24_eeg_fmri_raw) and dropbox (/Users/andrewfarkas/UFL Dropbox/Andrew Farkas/TwoDayFearConditioning/Data). If not, get the data from the share drive by opening a finder window and click : go > connect to server > smb://ICT-BA1-FS01.AD.UFL.EDU/AMRIS-AMRIS3TS$ > connect > KEil. There you will find the raw fMRI data. EEG and DAT files need to be taken from the EEG presentation and data recording computers.

- After double-checking that the participant is backed up in two locations, you can move that participant folder to the local raw_data folder for processing (/Users/andrewfarkas/research_data/gaborgen/raw_data). The folder structure needs to be how it is described above.

- Open the a001fMRI_prepro.R in RStudio, adjust the top lines as needed to process specific participants or Day1/Day2. The script should run with no problems, so if it fails, double check the structure of folders, paths, necessary tools (dcm2niix, afni, tcsh).

  - The script will save the results to the folder you specify. Each subject will take a few hours and many gigabytes of space. Only a subset of the participants will completely fill your hard drive and you'll run out of space.  

- The quickest way to check the the quality of the data after preprocessing is to open the html file located at GABORGEN24_XXX.results/QC_GABORGEN24_122/index.html. Double-click it to see if everything is lined up correctly, how much movement there was, and if there was a basic visual effect from the cues. Update the excel Mastersheet in the fMRI_preprocessing tab based on if the participant looks like they have useable data or not.

- Now back up the results. Take the massive preprocessed results folder, and the logs (the files that start with output.proc and proc) and put those in the AFNI_results folder (/Users/andrewfarkas/UFL Dropbox/Andrew Farkas/TwoDayFearConditioning/AFNI_results_andrewf). The folder is so large that we cannot back it up locally, but this is okay. The most important thing is that there are multiple backups of the raw data because the preprocessed results can always be regenerated. So that is why it is so important to have the raw data backed up on the dropbox and local hard drive. Andrew also has a copy of the raw data. 

- Lastly, move the nifti, json, BRIK, and HEAD files into the raw data folder for each participant. This will help with data reporting later, because since the software update, we no longer get individual folder for each scan.



## PREPRO_gaborgen_FinalVersion
This script does the initial preprocessing per participant.
It imposes restrictions to force the user to have the data in the correct format and to prevent incorrect analyses.
This script save one output for each step in the participant's folder

1. Format transformation
For the vdhr to set in first step we should be sure that we have the three files of EEG:
vdhr
eeg
vmrk
2. Chan Locations, Renaming T1, removing MRI artifacts, resampling

You need to write the path of the chanel locations files: '/Users/jcedielescobar/Documents/GitHub/gaborgen1/standard_1005.elc'  (31 channels)
Resampling: 500

3. Cardibalistic artifacts

4. Filtering
 pop_eegfiltnew, 3 to 40

5. Bad Channels (For the first time)
This is step opens each file and you need to check manually the bad channels with the scroll and with the spectra maps, you write the numbers and save the files. 
It will create a structure called badChannelsStruct.mat with the name of the subject and the bad channels numbers for be used in the interpolation step.

6. ICA
Type: Sobi

7. Eye components.

The first time you need to select the channels with the ICLabel tool but you need to check manually. You will need the eyecomponent.mat structure with the name of 
the subject and the components that need to be removed. Be sure the struct is in your main path. 

8. Interpolation 
Be sure you have the badChannelsStruct.mat struct in your main path

9. Reference 

10. CSD Transformation 
