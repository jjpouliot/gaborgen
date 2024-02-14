#!/bin/tcsh

afni_proc.py                                                             \
        -subj_id                  GABORGEN24_120_B0corr                                \
        -dsets                    /blue/akeil/jourdan.pouliot/GaborGen24/bidsGaborGen/sub-120/func/sub-120_task-gaborgen_bold.nii.gz  \
        -copy_anat                /blue/akeil/jourdan.pouliot/GaborGen24/bidsGaborGen/derivatives/afni_23.3.09/sub-120/ANAT_DEOB+orig.HEAD             \
        -blocks                   tshift align tlrc volreg mask blur scale regress     \
        -tcat_remove_first_trs    0                                      \
        -radial_correlate_blocks  tcat volreg                            \
	-blip_forward_dset	  /blue/akeil/jourdan.pouliot/GaborGen24/bidsGaborGen/sub-120/fmap/sub-120_dir-PA_epi.nii.gz                          \
	-blip_reverse_dset	  /blue/akeil/jourdan.pouliot/GaborGen24/bidsGaborGen/sub-120/fmap/sub-120_dir-AP_epi.nii.gz                         \
        -align_unifize_epi        yes                                  \
        -align_opts_aea           -cost lpc+ZZ                           \
                                  -giant_move                            \
                                  -check_flip                            \
        -tlrc_base                MNI152_2009_template.nii.gz            \
        -tlrc_NL_warp                                                   \
        -volreg_align_to          MIN_OUTLIER                            \
        -volreg_align_e2a                                                \
        -volreg_tlrc_warp                                                \
        -volreg_compute_tsnr	  yes                                    \
        -mask_epi_anat            yes                                    \
        -blur_size                4.0                                    \
        -regress_stim_times	  /blue/akeil/jourdan.pouliot/GaborGen24/bidsGaborGen/derivatives/afni_23.3.09/sub-120/stim_times.1D         \
        -regress_stim_labels	  stimulus                              \
        -regress_basis            'CSPLINzero(0,8,5)'                  \
        -regress_motion_per_run                                           \
        -regress_censor_motion    0.3                                     \
        -regress_censor_outliers  0.05                                     \
        -regress_reml_exec                                              \
        -regress_compute_fitts                                          \
        -regress_make_ideal_sum sum_ideal.1D                           \
        -regress_est_blur_epits                                           \
        -regress_est_blur_errts                                         \
        -execute


