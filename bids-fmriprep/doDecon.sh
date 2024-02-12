#!/bin/tcsh
# For the fMRIPrep tutorial, copy and paste this into the "func" directory of ${subj} in the "derivatives/fmriprep" folder
# and type: "tcsh doDecon.sh ${subj}"

if ( $#argv > 0 ) then
    set subj = $argv[1]
else
    set subj = 120_motCens
endif

3dDeconvolve -force_TR 2.0     \
    -input r*_scale.nii                            \
    -censor motion_120_censor.1D                                         \
    -mask sub-120_task-gaborgen_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz	     \
    -polort 15                                                                \
    -num_stimts 7                                                           \
    -stim_times 1 stimuli/stim_times.1D 'CSPLINzero(0,8,5)'                          \
    -stim_label 1 stimulus                                                  \
    -stim_file 2 trans_x_run01.txt'[0]' -stim_base 2 -stim_label 2 trans_x_01   \
    -stim_file 3 trans_y_run01.txt'[0]' -stim_base 3 -stim_label 3 trans_y_01  \
    -stim_file 4 trans_z_run01.txt'[0]' -stim_base 4 -stim_label 4 trans_z_01    \
    -stim_file 5 rot_x_run01.txt'[0]' -stim_base 5 -stim_label 5 rot_x_01     \
    -stim_file 6 rot_y_run01.txt'[0]' -stim_base 6 -stim_label 6 rot_y_01     \
    -stim_file 7 rot_z_run01.txt'[0]' -stim_base 7 -stim_label 7 rot_z_01     \
    -jobs 7                                                                  \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                                  \
    -xjpeg X.jpg                                                             \
    -fitts fitts.$subj                                                       \
    -errts errts.${subj}                                                     \
    -bucket stats.$subj
