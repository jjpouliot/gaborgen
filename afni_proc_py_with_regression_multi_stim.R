setwd(where_results_should_be_saved)

for(directory_index in 1:length(participant_directories)) {
  
  current_participant_found_structural_directories <- list.files(
    paste0(participant_directories[directory_index],'/fMRI'),
    pattern = 'T1_MPRAGE_SAG_P2_ISO', full.names = T)
  
  current_structural_directory <- 
    current_participant_found_structural_directories[
      which.max(as.numeric(
        sub(".*?(\\d+)$", "\\1", 
            current_participant_found_structural_directories)))]
  
  current_participant_found_functional_directories <- list.files(
    paste0(participant_directories[directory_index],'/fMRI'),
    pattern = 'BOLD-EPI-CMRR-2S', full.names = T)
  
  current_functional_directory <- current_participant_found_functional_directories[
    which.max(as.numeric(
      sub(".*?(\\d+)$", "\\1", 
          current_participant_found_functional_directories)))]
  
  current_participant_found_blip_forward <- list.files(
    paste0(participant_directories[directory_index],'/fMRI'),
    pattern = 'CMMR-DISTMAP_AP', full.names = T)
  
  current_blip_forward_directory <- current_participant_found_blip_forward[
    which.max(as.numeric(
      sub(".*?(\\d+)$", "\\1", 
          current_participant_found_blip_forward)))]
  
  current_participant_found_blip_reverse <- list.files(
    paste0(participant_directories[directory_index],'/fMRI'),
    pattern = 'CMMR-DISTMAP_PA', full.names = T)
  
  current_blip_reverse_directory <- current_participant_found_blip_reverse[
    which.max(as.numeric(
      sub(".*?(\\d+)$", "\\1", 
          current_participant_found_blip_reverse)))]
  
  current_structural_path <- list.files(current_structural_directory, 
                                        pattern = '.HEAD$',
                                        full.names = T)
  
  current_functional_path <- list.files(current_functional_directory, 
                                        pattern = '.HEAD$',
                                        full.names = T)
  
  current_blip_forward_path <- list.files(current_blip_forward_directory, 
                                          pattern = '.HEAD$', 
                                          full.names = T)
  
  current_blip_reverse_path <- list.files(current_blip_reverse_directory, 
                                          pattern = '.HEAD$',
                                          full.names = T)
  
  current_stim_onsets_path <- paste(list.files(participant_directories[directory_index], 
                                               pattern = '_stim_times.1D$',
                                               full.names = T), 
                                    collapse = " ")
  
  stim_names <- gsub(pattern = "_stim_times\\.1D", 
                     replacement = "", 
                     paste(list.files(participant_directories[directory_index], 
                                      pattern = '_stim_times.1D$'), 
                           collapse = " "))
  
  afni_proc_py_script <- paste0('python3 ~/abin/afni_proc.py
        -subj_id ', paste0(basename(participant_directories[directory_index]), " "),
       '-dsets ', paste0(current_functional_path," "),
       '-copy_anat ', paste0(current_structural_path, " "), 
       '-blocks tshift align tlrc volreg mask blur scale regress
        -tcat_remove_first_trs 0
        -radial_correlate_blocks  tcat volreg
        -blip_forward_dset ', paste0(current_blip_forward_path, " "),
       '-blip_reverse_dset ', paste0(current_blip_reverse_path, " "),
        '-align_unifize_epi        yes
        -align_opts_aea           -cost lpc+ZZ  
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
        '-regress_stim_labels ', paste0(stim_names, " "),                              
        '-regress_basis      "BLOCK(2,1)"                 
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
