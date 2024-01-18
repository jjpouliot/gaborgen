participant_index <- 1

current_structural_path <- paste0(participant_directories[participant_index],
                                  '/fMRI/T1_MPRAGE_SAG_P2_ISO_0005')

current_functional_path <- paste0(participant_directories[participant_index],
                                  '/fMRI/BOLD-EPI-CMRR-2S_0006')

current_structural_NIFTI_path <- list.files(current_structural_path, 
                                            pattern = '.nii$', 
                                            full.names = T)

current_functional_NIFTI_path <- list.files(current_functional_path, 
                                            pattern = '.nii$', 
                                            full.names = T)

paste0('afni_proc.py
        -subj_id ', paste("GABORGEN24_", participants_to_preprocess[participant_index]),
       '-dsets ', current_functional_NIFTI_path,
       '-copy_anat ', current_functional_path,
       '-blocks tshift align tlrc volreg mask blur scale regress
        -tcat_remove_first_trs 0
        -radial_correlate_blocks  tcat volreg
        -align_unifize_epi        yes
        -align_opts_aea           -cost lpc+ZZ   
                                  -giant_move  
                                  -check_flip
        -tlrc_base                MNI152_2009_template.nii.gz
        -tlrc_NL_warp
        -volreg_align_to          MIN_OUTLIER
        -volreg_align_e2a
        -volreg_tlrc_warp
        -volreg_compute_tsnr      yes
        -mask_epi_anat            yes
        -blur_size                4.0
        -regress_stim_times       /Users/andrewfarkas/research_data/gaborgen/raw_data/GABORGEN24_120/stim_times.1D         \
        -regress_stim_labels      stimulus                              \
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



"