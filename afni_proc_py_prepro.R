setwd(where_results_should_be_saved)

for(participant_index in 1:length(participant_directories)) {
  
  current_structural_path <- paste0(participant_directories[participant_index],
                                    '/fMRI/T1_MPRAGE_SAG_P2_ISO_0005')
  
  current_functional_path <- paste0(participant_directories[participant_index],
                                    '/fMRI/BOLD-EPI-CMRR-2S_0006')
  
  current_stim_onsets_path <- paste0(participant_directories[participant_index],
                                     '/stim_times.1D')
  
  current_structural_NIFTI_path <- list.files(current_structural_path, 
                                              pattern = '.nii$', 
                                              full.names = T)
  
  current_functional_NIFTI_path <- list.files(current_functional_path, 
                                              pattern = '.nii$', 
                                              full.names = T)
  
  afni_proc_py_script <- paste0('afni_proc.py
        -subj_id ', paste0("GABORGEN24_", participants_to_preprocess[participant_index], " "),
                                '-dsets ', paste0(current_functional_NIFTI_path," "),
                                '-copy_anat ', paste0(current_structural_NIFTI_path, " "), 
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
        
        -regress_stim_times ', paste0(current_stim_onsets_path, " "),
                                '-regress_stim_labels      stimulus                              
        -regress_basis      "CSPLINzero(0,8,5)"                 
        -regress_motion_per_run                                           
        -regress_censor_motion    0.3                                     
        -regress_censor_outliers  0.05                                     
        -regress_reml_exec                                              
        -regress_compute_fitts                                          
        -regress_make_ideal_sum sum_ideal.1D                           
        -regress_est_blur_epits                                           
        -regress_est_blur_errts                                         
        -execute')
  
  afni_proc_py_script <- gsub(pattern = "\n", replacement = "", x = afni_proc_py_script)
  
  system2('tcsh', 
          args = c('-c', 
                   shQuote(afni_proc_py_script)))
  
}

setwd(local_git_directory)
