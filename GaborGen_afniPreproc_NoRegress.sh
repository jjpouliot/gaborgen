#!/bin/tcsh


afni_proc.py                                                             \
        -subj_id                  afni103                                \
        -copy_anat                afni103/sub-03_acq-t1mpragesagp2iso_run-1_T1w.nii.gz                        \
        -dsets                    afni103/sub-03_task-BOLDEPIcmrr2s_run-1_echo-1_bold.nii.gz                \
        -blocks                   tshift align tlrc volreg mask blur     \
                                  scale                                  \
        -radial_correlate_blocks  tcat volreg                            \
        -tcat_remove_first_trs    0                                      \
        -align_unifize_epi        yes                                  \
        -align_opts_aea           -cost lpc+ZZ                           \
                                  -giant_move                            \
                                  -check_flip                            \
        -tlrc_base                MNI152_2009_template.nii.gz            \
        -volreg_align_to          MIN_OUTLIER                            \
        -volreg_align_e2a                                                \
        -volreg_tlrc_warp                                                \
        -volreg_compute_tsnr      yes                                    \
        -mask_epi_anat            yes                                    \
        -blur_size                4.0                                    \
        -execute
