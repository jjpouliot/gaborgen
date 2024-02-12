#!/bin/tcsh
# For the fMRIPrep tutorial, copy and paste this into the "func" directory of ${subj} in the "derivatives/fmriprep" folder
# and type: "tcsh doDecon.sh ${subj}"

if ( $#argv > 0 ) then
    set subj = $argv[1]
else
    set subj = 118_CSPlusVsGSs
endif

3dDeconvolve -force_TR 2.0                         \
    -input r*_scale.nii                            \
    -censor motion_118_censor.1D                                         \
    -mask sub-118_task-gaborgen_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz						     \
    -polort 15                                                                \
    -num_stimts 18                                                           \
    -stim_times 1 stimuli/acqCSplus_times.1D 'CSPLINzero(0,12,7)'                          \
    -stim_label 1 acqCSplus                                                  \
    -stim_times 2 stimuli/acqGS1_times.1D 'CSPLINzero(0,12,7)'                          \
    -stim_label 2 acqGS1                                                  \
    -stim_times 3 stimuli/acqGS2_times.1D 'CSPLINzero(0,12,7)'                          \
    -stim_label 3 acqGS2                                                  \
    -stim_times 4 stimuli/acqGS3_times.1D 'CSPLINzero(0,12,7)'                          \
    -stim_label 4 acqGS3                                                  \
    -stim_times 5 stimuli/habCSplus_times.1D 'CSPLINzero(0,12,7)'           \
    -stim_label 5 habCSplus                                                      \
    -stim_times 6 stimuli/habGS1_times.1D 'CSPLINzero(0,12,7)'                   \
    -stim_label 6 habGS1                                                         \
    -stim_times 7 stimuli/habGS2_times.1D 'CSPLINzero(0,12,7)'                   \
    -stim_label 7 habGS2                                                         \
    -stim_times 8 stimuli/habGS3_times.1D 'CSPLINzero(0,12,7)'                   \
    -stim_label 8 habGS3                                                         \
    -stim_times 9 stimuli/extCSplus_times.1D 'CSPLINzero(0,12,7)'                   \
    -stim_label 9 extCSplus                                                         \
    -stim_times 10 stimuli/extGS1_times.1D 'CSPLINzero(0,12,7)'                   \
    -stim_label 10 extGS1                                                         \
    -stim_times 11 stimuli/extGS2_times.1D 'CSPLINzero(0,12,7)'                   \
    -stim_label 11 extGS2                                                         \
    -stim_times 12 stimuli/extGS3_times.1D 'CSPLINzero(0,12,7)'                   \
    -stim_label 12 extGS3                                                         \
    -stim_file 13 trans_x_run01.txt'[0]' -stim_base 13 -stim_label 13 trans_x_01   \
    -stim_file 14 trans_y_run01.txt'[0]' -stim_base 14 -stim_label 14 trans_y_01  \
    -stim_file 15 trans_z_run01.txt'[0]' -stim_base 15 -stim_label 15 trans_z_01    \
    -stim_file 16 rot_x_run01.txt'[0]' -stim_base 16 -stim_label 16 rot_x_01     \
    -stim_file 17 rot_y_run01.txt'[0]' -stim_base 17 -stim_label 17 rot_y_01     \
    -stim_file 18 rot_z_run01.txt'[0]' -stim_base 18 -stim_label 18 rot_z_01     \
    -jobs 10                                                                  \
    -gltsym 'SYM: 3*acqCSplus -acqGS1 -acqGS2 -acqGS3'                                      \
    -glt_label 1 acqCSplus-GSs                                                   \
    -gltsym 'SYM: 3*habCSplus -habGS1 -habGS2 -habGS3'                                      \
    -glt_label 2 habCSplus-GSs                                                   \
    -gltsym 'SYM: 3*extCSplus -extGS1 -extGS2 -extGS3'                                      \
    -glt_label 3 extCSplus-GSs                                                   \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                                  \
    -x1D_uncensored X.nocensor.xmat.1D                                       \
    -fitts fitts.$subj                                                       \
    -errts errts.${subj}                                                     \
    -bucket stats.$subj

