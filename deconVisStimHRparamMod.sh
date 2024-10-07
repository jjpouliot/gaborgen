#!/bin/tcsh

set subj = $argv[1]

3dDeconvolve -input pb05.GABORGEN24_${subj}.r01.scale+tlrc.HEAD          \
    -polort A           \
    -censor censor_GABORGEN24_${subj}_combined_2.1D            \
    -mask full_mask.GABORGEN24_${subj}+tlrc.HEAD                    \
    -num_stimts 7                   \
    -stim_times_AM2 1 stimuli/acqAllStims-HRdc.1D 'BLOCK(2,1)'               \
    -stim_label 1 visStim+HRdc            \
    -stim_file 2 motion_demean.1D'[0]' -stim_base 2 -stim_label 2 roll   \
    -stim_file 3 motion_demean.1D'[1]' -stim_base 3 -stim_label 3 pitch  \
    -stim_file 4 motion_demean.1D'[2]' -stim_base 4 -stim_label 4 yaw    \
    -stim_file 5 motion_demean.1D'[3]' -stim_base 5 -stim_label 5 dS     \
    -stim_file 6 motion_demean.1D'[4]' -stim_base 6 -stim_label 6 dL     \
    -stim_file 7 motion_demean.1D'[5]' -stim_base 7 -stim_label 7 dP  \
    -iresp 1 GABORGEN24_${subj}.BLOCK.all.IRF                        \
    -jobs 10            \
    -fout                            \
    -tout                             \
    -bucket GABORGEN24_${subj}.BLOCK.all.stats                       \
    -xsave                                                   \
    -xout                                                   \
    -xjpeg X.jpg                                           \
    -x1D X.xmat.1D                                     