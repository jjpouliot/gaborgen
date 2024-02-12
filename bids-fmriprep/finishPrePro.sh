#!/bin/tcsh
#Finishes the prepro pipeline by smoothing and scaling the images

3dMerge -1blur_fwhm 4.0 -doall -prefix r01_blur.nii sub-118_task-gaborgen_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz

3dTstat -prefix rm.mean_r01.nii r01_blur.nii

3dcalc -a r01_blur.nii -b rm.mean_r01.nii -c sub-118_task-gaborgen_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz \
-expr 'c * min(200, a/b*100)*step(a)*step(b)' -prefix r01_scale.nii

