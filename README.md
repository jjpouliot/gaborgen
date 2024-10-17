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


If ARM, open terminal and run:
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

## 001fMRI_prepro.R
This script does the initial preprocessing per participant. It imposes restrictions to force the user to have the data in the correct format and to prevent incorrect analyses. It also forces the user to put their data and local git directory paths at the top of the script. The script forces the use of the tcsh terminal which is recommended for AFNI analyses. The system2 R function is used to force tcsh terminal calls.

